"""Reference cases for genotype encoding.

Every expected mask is derived from the pair numbering: copy 0 holding
base ``a`` and copy 1 holding ``b`` is bit ``5*a + b``, with A=0, C=1, G=2, T=3
and slot 4 for ABSENT.
"""

from __future__ import annotations

import copy
import itertools
import pickle
import re
from collections.abc import Mapping

import numpy as np
import pandas as pd
import pytest

from lineage.reconstruction.encoding import _Calls, align, encode, haplotypes
from lineage.types import ABSENT, UNKNOWN, Haplotypes, Loci, Sex

# --------------------------------------------------------------------------------------
# Masks and kinds written out as literals
# --------------------------------------------------------------------------------------

A, C, G, T = 0, 1, 2, 3
NO_COPY = 4  # the ABSENT slot of a pair


def pair_bit(a: int, b: int) -> int:
    """Bit for copy 0 holding ``a`` and copy 1 holding ``b``."""
    return 1 << (5 * a + b)


# Two copies, nothing known: the 16 base/base pairs, bits 0-3, 5-8, 10-13 and 15-18.
# 15 + (15 << 5) + (15 << 10) + (15 << 15) = 15 + 480 + 15_360 + 491_520.
TWO_ANY = 507_375
# One copy, nothing known: (base, ABSENT), bits 4, 9, 14 and 19.
# 16 + 512 + 16_384 + 524_288.
ONE_ANY = 541_200
AG = 1_028  # unordered AG: bit 2 (A,G) = 4 and bit 10 (G,A) = 1_024
AG_ORDERED = 4  # bit 2 (A,G): copy 0 holds A
GA_ORDERED = 1_024  # bit 10 (G,A): copy 0 holds G
# A on one copy, the other unknown, no order: bits 0-3 (A,x) and 5, 10, 15 (x,A).
# 15 + 32 + 1_024 + 32_768; seven pairs.
A_ON_ONE_COPY = 33_839
CT = 65_792  # bit 8 (C,T) = 256 and bit 16 (T,C) = 65_536
CT_ORDERED = 256  # bit 8 (C,T)
TC_ORDERED = 65_536  # bit 16 (T,C)
AT = 32_776  # bit 3 (A,T) = 8 and bit 15 (T,A) = 32_768
TA_ORDERED = 32_768  # bit 15 (T,A)
CA_ORDERED = 32  # bit 5 (C,A)
GT = 139_264  # bit 13 (G,T) = 8_192 and bit 17 (T,G) = 131_072
AA = 1  # bit 0 (A,A)
CC = 64  # bit 6 (C,C)
GG = 4_096  # bit 12 (G,G)
TT = 262_144  # bit 18 (T,T)
A_ONE = 16  # bit 4 (A, ABSENT)
C_ONE = 512  # bit 9 (C, ABSENT)
G_ONE = 16_384  # bit 14 (G, ABSENT)
T_ONE = 524_288  # bit 19 (T, ABSENT)
COPY0_A = 15  # bits 0-3: (A,A), (A,C), (A,G), (A,T); copy 0 is A, copy 1 unknown
COPY0_G = 15_360  # bits 10-13: (G,A), (G,C), (G,G), (G,T); copy 0 is G
COPY1_A = 33_825  # bits 0, 5, 10, 15: (A,A), (C,A), (G,A), (T,A); copy 1 is A

# Kinds: (genotype, phase, origin) with 0 missing, 1 observed, 2 forced,
# 3 estimated and 4 asserted.
OBSERVED = (1, 0, 0)
ASSERTED = (1, 4, 0)
MISSING = (0, 0, 0)

X_POS = 3_000_000  # on X, between the two pseudoautosomal regions


def test_literal_masks_match_their_bits():
    bases = range(4)
    assert TWO_ANY == sum(pair_bit(a, b) for a in bases for b in bases)
    assert ONE_ANY == sum(pair_bit(a, NO_COPY) for a in bases)
    assert AG == pair_bit(A, G) | pair_bit(G, A)
    assert AG_ORDERED == pair_bit(A, G)
    assert GA_ORDERED == pair_bit(G, A)
    pairs_with_a = {(A, b) for b in bases} | {(b, A) for b in bases}
    assert len(pairs_with_a) == 7
    assert A_ON_ONE_COPY == sum(pair_bit(a, b) for a, b in pairs_with_a)
    assert CT == pair_bit(C, T) | pair_bit(T, C)
    assert CT_ORDERED == pair_bit(C, T)
    assert TC_ORDERED == pair_bit(T, C)
    assert AT == pair_bit(A, T) | pair_bit(T, A)
    assert TA_ORDERED == pair_bit(T, A)
    assert CA_ORDERED == pair_bit(C, A)
    assert GT == pair_bit(G, T) | pair_bit(T, G)
    assert (AA, CC, GG, TT) == tuple(pair_bit(b, b) for b in bases)
    assert (A_ONE, C_ONE, G_ONE, T_ONE) == tuple(pair_bit(b, NO_COPY) for b in bases)
    assert COPY0_A == sum(pair_bit(A, b) for b in bases)
    assert COPY0_G == sum(pair_bit(G, b) for b in bases)
    assert COPY1_A == sum(pair_bit(b, A) for b in bases)
    assert (UNKNOWN, ABSENT) == (-1, -2)


# --------------------------------------------------------------------------------------
# Input builders
# --------------------------------------------------------------------------------------


def row(
    rsid: str,
    chrom: str | None,
    pos: int,
    genotype: str | None,
    *,
    phased: bool = False,
    phase_set: int | None = None,
    conflicted: bool = False,
) -> tuple:
    return (rsid, chrom, pos, genotype, phased, phase_set, conflicted)


def make_frame(rows: list[tuple]) -> pd.DataFrame:
    """A normalized snps frame with object ``chrom`` and ``genotype`` columns."""
    frame = pd.DataFrame(
        {
            "chrom": pd.Series([r[1] for r in rows], dtype=object),
            "pos": pd.Series(np.array([r[2] for r in rows], dtype=np.uint32)),
            "genotype": pd.Series([r[3] for r in rows], dtype=object),
        }
    )
    frame.index = pd.Index([r[0] for r in rows], dtype="str", name="rsid")
    return frame


def bool_series(frame: pd.DataFrame, values: list[bool]) -> pd.Series:
    return pd.Series(np.array(values, dtype=bool), index=frame.index.copy())


def set_series(frame: pd.DataFrame, values: list[int | None]) -> pd.Series:
    return pd.Series(pd.array(values, dtype="Int64"), index=frame.index.copy())


def build_inputs(rows: list[tuple]):
    frame = make_frame(rows)
    phased = bool_series(frame, [r[4] for r in rows])
    sets = set_series(frame, [r[5] for r in rows])
    conflicted = bool_series(frame, [r[6] for r in rows])
    return frame, phased, sets, conflicted


def run_encode(rows: list[tuple], *, sex: Sex | None, omit: tuple[str, ...] = ()):
    frame, phased, sets, conflicted = build_inputs(rows)
    series = {"phased": phased, "phase_set": sets, "genotype_conflicted": conflicted}
    for name in omit:
        series[name] = None
    return encode(frame, sex=sex, **series)


def make_calls(loci, allowed, phase_set, kinds) -> _Calls:
    def as_array(values):
        if isinstance(values, np.ndarray):
            return values
        return np.array(values, dtype=np.int64)

    return _Calls(
        loci,
        as_array(allowed),
        as_array(phase_set),
        as_array(kinds),
    )


def block(chrom, pos, rsid, allowed, kinds, *, phase_set=None) -> _Calls:
    n = len(pos)
    return make_calls(
        Loci(chrom, pos, rsid),
        allowed,
        [-1] * n if phase_set is None else phase_set,
        [list(k) for k in kinds],
    )


# --------------------------------------------------------------------------------------
# Result readers
# --------------------------------------------------------------------------------------

DIAGNOSTIC_COLUMNS = ["code", "chrom", "pos", "start", "end", "detail"]
AGGREGATED_CODES = {"non_base_genotype", "ploidy_mismatch", "excluded_region"}
ROW_COUNT = re.compile(r"(\d+) (rows?)(?![A-Za-z0-9_])")


def summary(calls: _Calls) -> dict[str, object]:
    return {
        "chrom": calls.loci.chrom,
        "pos": calls.loci.pos.tolist(),
        "rsid": calls.loci.rsid.tolist(),
        "allowed": calls.allowed.tolist(),
        "phase_set": calls.phase_set.tolist(),
        "kinds": calls.kinds.tolist(),
    }


def summaries(blocks: Mapping[str, _Calls]) -> list[tuple[str, dict[str, object]]]:
    return [(chrom, summary(calls)) for chrom, calls in blocks.items()]


def assert_block(calls, *, chrom, pos, rsid, allowed, kinds, phase_set=None) -> None:
    n = len(pos)
    assert isinstance(calls, _Calls)
    assert isinstance(calls.loci, Loci)
    assert calls.loci.chrom == chrom
    assert calls.loci.pos.dtype == np.uint32
    assert calls.allowed.dtype == np.uint32
    assert calls.phase_set.dtype == np.int32
    assert calls.kinds.dtype == np.uint8
    assert calls.allowed.shape == calls.phase_set.shape == (n,)
    assert calls.kinds.shape == (n, 3)
    assert calls.loci.pos.tolist() == list(pos)
    assert calls.loci.rsid.tolist() == list(rsid)
    assert calls.allowed.tolist() == list(allowed)
    assert calls.kinds.tolist() == [list(k) for k in kinds]
    assert calls.phase_set.tolist() == (
        [-1] * n if phase_set is None else list(phase_set)
    )


def leading_row_count(detail: object) -> int:
    assert isinstance(detail, str), detail
    match = ROW_COUNT.match(detail)
    assert match is not None, detail
    count = int(match.group(1))
    assert match.group(2) == ("row" if count == 1 else "rows"), detail
    return count


def assert_diagnostics_schema(diagnostics: pd.DataFrame) -> None:
    assert list(diagnostics.columns) == DIAGNOSTIC_COLUMNS
    assert isinstance(diagnostics.index, pd.RangeIndex)
    assert diagnostics.index.equals(pd.RangeIndex(len(diagnostics)))
    for column in ("code", "chrom", "detail"):
        assert diagnostics[column].dtype == "str", column
    for column in ("pos", "start", "end"):
        assert diagnostics[column].dtype == "Int64", column


def _optional_int(value: object) -> int | None:
    return None if pd.isna(value) else int(value)


def diagnostic_rows(diagnostics: pd.DataFrame) -> list[tuple]:
    """(code, chrom, pos, start, end, leading row count) per row, in order."""
    assert_diagnostics_schema(diagnostics)
    columns = {name: diagnostics[name].tolist() for name in DIAGNOSTIC_COLUMNS}
    rows = []
    for i in range(len(diagnostics)):
        code = columns["code"][i]
        detail = columns["detail"][i]
        count = leading_row_count(detail) if code in AGGREGATED_CODES else None
        rows.append(
            (
                code,
                columns["chrom"][i],
                _optional_int(columns["pos"][i]),
                _optional_int(columns["start"][i]),
                _optional_int(columns["end"][i]),
                count,
            )
        )
    return rows


def locus_diag(code: str, chrom: str, pos: int) -> tuple:
    return (code, chrom, pos, None, None, None)


def region_diag(code: str, chrom: str, start: int, end: int, rows: int) -> tuple:
    return (code, chrom, None, start, end, rows)


def hap_columns(haps: Haplotypes) -> list[tuple[int, int]]:
    return [(int(haps[0, i]), int(haps[1, i])) for i in range(haps.shape[1])]


def assert_read_only(array: np.ndarray) -> None:
    assert isinstance(array, np.ndarray)
    assert not array.flags.writeable
    with pytest.raises(ValueError):
        array[...] = array.copy()
    base = array.base
    while base is not None:
        if isinstance(base, np.ndarray):
            assert not base.flags.writeable
            base = base.base
        else:
            assert memoryview(base).readonly
            base = None


def assert_read_only_mapping(mapping: Mapping, key: str) -> None:
    assert isinstance(mapping, Mapping)
    with pytest.raises(TypeError):
        mapping[key] = mapping[key]  # type: ignore[index]
    with pytest.raises(TypeError):
        del mapping[key]  # type: ignore[attr-defined]
    with pytest.raises(TypeError):
        mapping["new"] = mapping[key]  # type: ignore[index]
    assert key in mapping


# --------------------------------------------------------------------------------------
# 1. The allowed pairs of each kind of call
# --------------------------------------------------------------------------------------


def _encoded_one_row(chrom, pos, genotype, *, sex, phased=False) -> _Calls:
    blocks, _ = run_encode([row("rs1", chrom, pos, genotype, phased=phased)], sex=sex)
    return blocks[chrom]


def _a_known_on_one_copy() -> _Calls:
    # encode cannot produce a partial two-copy call; build it with a forced genotype.
    return make_calls(Loci("1", [1_000]), [A_ON_ONE_COPY], [-1], [[2, 0, 0]])


ALLOWED_TABLE = [
    pytest.param(
        lambda: _encoded_one_row("1", 1_000, None, sex=Sex.FEMALE),
        TWO_ANY,
        16,
        MISSING,
        -1,
        (UNKNOWN, UNKNOWN),
        id="two-copies-nothing-known",
    ),
    pytest.param(
        lambda: _encoded_one_row("1", 1_000, "AG", sex=Sex.FEMALE),
        AG,
        2,
        OBSERVED,
        -1,
        (UNKNOWN, UNKNOWN),
        id="unordered-AG",
    ),
    pytest.param(
        lambda: _encoded_one_row("1", 1_000, "AG", sex=Sex.FEMALE, phased=True),
        AG_ORDERED,
        1,
        ASSERTED,
        0,
        (A, G),
        id="AG-known-order",
    ),
    pytest.param(
        _a_known_on_one_copy,
        A_ON_ONE_COPY,
        7,
        (2, 0, 0),
        -1,
        (UNKNOWN, UNKNOWN),
        id="A-on-one-copy-other-unknown",
    ),
    pytest.param(
        lambda: _encoded_one_row("X", X_POS, "A", sex=Sex.MALE),
        A_ONE,
        1,
        OBSERVED,
        -1,
        (A, ABSENT),
        id="one-copy-A",
    ),
    pytest.param(
        lambda: _encoded_one_row("X", X_POS, None, sex=Sex.MALE),
        ONE_ANY,
        4,
        MISSING,
        -1,
        (UNKNOWN, ABSENT),
        id="one-copy-nothing-known",
    ),
]


@pytest.mark.parametrize(
    ("build", "mask", "n_pairs", "kinds", "phase_set", "copies"), ALLOWED_TABLE
)
def test_allowed_table_row(build, mask, n_pairs, kinds, phase_set, copies):
    calls = build()
    assert calls.allowed.dtype == np.uint32
    assert calls.allowed.tolist() == [mask]
    assert mask.bit_count() == n_pairs
    assert calls.kinds.tolist() == [list(kinds)]
    assert calls.phase_set.tolist() == [phase_set]
    haps = haplotypes(calls)
    assert isinstance(haps, np.ndarray)
    assert haps.dtype == np.int8
    assert haps.shape == (2, 1)
    assert hap_columns(haps) == [copies]


def test_haplotypes_give_a_copy_only_where_every_pair_agrees():
    calls = make_calls(
        Loci("1", [100, 200, 300, 400, 500]),
        [
            COPY0_A,  # every pair has A on copy 0; copy 1 is A, C, G or T
            COPY1_A,  # every pair has A on copy 1
            TT,  # (T,T)
            pair_bit(A, C) | pair_bit(A, G),  # bits 1, 2: copy 0 A, copy 1 C or G
            AG,  # (A,G) or (G,A): neither copy is settled
        ],
        [0, 0, -1, 0, -1],
        [[2, 2, 0], [2, 2, 0], [3, 0, 0], [2, 2, 0], [1, 0, 0]],
    )
    haps = haplotypes(calls)
    assert haps.dtype == np.int8
    assert haps.shape == (2, 5)
    assert haps[0].tolist() == [A, UNKNOWN, T, A, UNKNOWN]
    assert haps[1].tolist() == [UNKNOWN, A, T, UNKNOWN, UNKNOWN]
    again = haplotypes(calls)
    assert again is not haps
    assert not np.shares_memory(again, haps)


def test_haplotypes_of_one_copy_rows_mark_the_second_copy_absent():
    calls = make_calls(
        Loci("X", [3_000_000, 4_000_000, 5_000_000]),
        [G_ONE, ONE_ANY, A_ONE | C_ONE],  # bit 14; bits 4, 9, 14, 19; bits 4 and 9
        [-1, -1, -1],
        [[1, 0, 0], [0, 0, 0], [3, 0, 0]],
    )
    haps = haplotypes(calls)
    assert haps.dtype == np.int8
    assert haps[0].tolist() == [G, UNKNOWN, UNKNOWN]
    assert haps[1].tolist() == [ABSENT, ABSENT, ABSENT]


# --------------------------------------------------------------------------------------
# 2. encode on single rows
# --------------------------------------------------------------------------------------

SINGLE_ROWS = [
    # chrom, pos, sex, genotype, phased, allowed, kinds, phase_set, diagnostics
    pytest.param(
        "1",
        1_000,
        Sex.FEMALE,
        "AG",
        False,
        AG,
        OBSERVED,
        -1,
        [],
        id="autosome-unordered-AG",
    ),
    pytest.param(
        "1", 1_000, None, "AG", False, AG, OBSERVED, -1, [], id="autosome-unknown-sex"
    ),
    pytest.param(
        "1",
        1_000,
        Sex.MALE,
        "CC",
        False,
        CC,
        OBSERVED,
        -1,
        [],
        id="autosome-homozygous-CC",
    ),
    # Phased: copy 0 holds the first character; one implicit set, numbered 0.
    pytest.param(
        "1",
        1_000,
        Sex.FEMALE,
        "AG",
        True,
        AG_ORDERED,
        ASSERTED,
        0,
        [],
        id="autosome-phased-AG",
    ),
    pytest.param(
        "1",
        1_000,
        Sex.MALE,
        "TC",
        True,
        TC_ORDERED,
        ASSERTED,
        0,
        [],
        id="autosome-phased-TC",
    ),
    pytest.param(
        "1",
        1_000,
        Sex.MALE,
        "A",
        False,
        TWO_ANY,
        MISSING,
        -1,
        [region_diag("ploidy_mismatch", "1", 1_000, 1_000, 1)],
        id="autosome-one-character",
    ),
    pytest.param(
        "X",
        X_POS,
        Sex.FEMALE,
        "A",
        False,
        TWO_ANY,
        MISSING,
        -1,
        [region_diag("ploidy_mismatch", "X", X_POS, X_POS, 1)],
        id="female-x-one-character",
    ),
    pytest.param(
        "X",
        X_POS,
        Sex.FEMALE,
        "CT",
        False,
        CT,
        OBSERVED,
        -1,
        [],
        id="female-x-two-bases",
    ),
    pytest.param(
        "X", X_POS, Sex.MALE, "A", False, A_ONE, OBSERVED, -1, [], id="male-x-A"
    ),
    pytest.param(
        "X", X_POS, Sex.MALE, "AA", False, A_ONE, OBSERVED, -1, [], id="male-x-AA"
    ),
    pytest.param(
        "X", X_POS, Sex.MALE, "TT", False, T_ONE, OBSERVED, -1, [], id="male-x-TT"
    ),
    pytest.param(
        "X",
        X_POS,
        Sex.MALE,
        "AG",
        False,
        ONE_ANY,
        MISSING,
        -1,
        [region_diag("ploidy_mismatch", "X", X_POS, X_POS, 1)],
        id="male-x-AG",
    ),
    *[
        pytest.param(
            "1",
            1_000,
            Sex.FEMALE,
            value,
            False,
            TWO_ANY,
            MISSING,
            -1,
            [region_diag("non_base_genotype", "1", 1_000, 1_000, 1)],
            id=f"autosome-non-base-{value or 'empty'}",
        )
        for value in ["II", "DI", "D", "ag", "AAA", "0", ""]
    ],
    *[
        pytest.param(
            "X",
            X_POS,
            Sex.MALE,
            value,
            False,
            ONE_ANY,
            MISSING,
            -1,
            [region_diag("non_base_genotype", "X", X_POS, X_POS, 1)],
            id=f"male-x-non-base-{value}",
        )
        for value in ["D", "a", "AAA"]
    ],
    pytest.param(
        "1",
        1_000,
        Sex.FEMALE,
        None,
        False,
        TWO_ANY,
        MISSING,
        -1,
        [],
        id="autosome-null",
    ),
    pytest.param(
        "X", X_POS, Sex.MALE, None, False, ONE_ANY, MISSING, -1, [], id="male-x-null"
    ),
]


@pytest.mark.parametrize(
    (
        "chrom",
        "pos",
        "sex",
        "genotype",
        "phased",
        "allowed",
        "kinds",
        "phase_set",
        "expected_diagnostics",
    ),
    SINGLE_ROWS,
)
def test_encode_single_row(
    chrom, pos, sex, genotype, phased, allowed, kinds, phase_set, expected_diagnostics
):
    blocks, diagnostics = run_encode(
        [row("rs1", chrom, pos, genotype, phased=phased)], sex=sex
    )
    assert list(blocks) == [chrom]
    assert_block(
        blocks[chrom],
        chrom=chrom,
        pos=[pos],
        rsid=["rs1"],
        allowed=[allowed],
        kinds=[kinds],
        phase_set=[phase_set],
    )
    assert diagnostic_rows(diagnostics) == expected_diagnostics


def test_encode_returns_blocks_in_natural_order():
    rows = [
        row("rs1", "X", X_POS, "AG"),
        row("rs2", "10", 5_000, "CC"),
        row("rs3", "2", 5_000, "AA"),
        row("rs4", "1", 5_000, "GG"),
    ]
    female, _ = run_encode(rows, sex=Sex.FEMALE)
    assert list(female) == ["1", "2", "10", "X"]
    assert [female[c].allowed.tolist() for c in female] == [[GG], [AA], [CC], [AG]]
    unknown, diagnostics = run_encode(rows, sex=None)
    assert list(unknown) == ["1", "2", "10"]
    assert diagnostic_rows(diagnostics) == [
        region_diag("excluded_region", "X", X_POS, X_POS, 1)
    ]


def test_encode_without_series_treats_every_row_as_unphased_and_unflagged():
    frame = make_frame(
        [
            row("rs1", "1", 100, "AG"),
            row("rs2", "1", 200, None),
            row("rs3", "1", 200, "CT"),
        ]
    )
    blocks, diagnostics = encode(frame, sex=Sex.FEMALE)
    assert_block(
        blocks["1"],
        chrom="1",
        pos=[100, 200],
        rsid=["rs1", "rs2;rs3"],
        allowed=[AG, CT],
        kinds=[OBSERVED, OBSERVED],
    )
    assert diagnostic_rows(diagnostics) == []


# --------------------------------------------------------------------------------------
# 3. Rows at one position
# --------------------------------------------------------------------------------------

# The control at 2000 is (C,C), bit 6 = 64, observed, and must not be touched.
CONTROL = row("rs99", "1", 2_000, "CC")

SAME_POSITION = [
    pytest.param(
        [row("rs1", "1", 1_000, "AG"), row("rs2", "1", 1_000, "GA")],
        AG,
        OBSERVED,
        "rs1;rs2",
        [],
        id="agreeing-aliases",
    ),
    pytest.param(
        [row("rs1", "1", 1_000, "AG"), row("rs1", "1", 1_000, "AG")],
        AG,
        OBSERVED,
        "rs1",
        [],
        id="repeated-row",
    ),
    pytest.param(
        [row("rs3", "1", 1_000, "CT"), row("rs4", "1", 1_000, None)],
        CT,
        OBSERVED,
        "rs3;rs4",
        [],
        id="call-with-no-call-alias",
    ),
    pytest.param(
        # Unflagged null alias permits the call; "rs10" sorts before "rs9" as text.
        [row("rs10", "1", 1_000, None), row("rs9", "1", 1_000, "GG")],
        GG,
        OBSERVED,
        "rs10;rs9",
        [],
        id="unflagged-null-alias-permits-call",
    ),
    pytest.param(
        [row("rs1", "1", 1_000, "AG"), row("rs2", "1", 1_000, "II")],
        AG,
        OBSERVED,
        "rs1;rs2",
        [region_diag("non_base_genotype", "1", 1_000, 1_000, 1)],
        id="call-with-non-base-alias",
    ),
    pytest.param(
        [row("rs1", "1", 1_000, "AG"), row("rs2", "1", 1_000, "A")],
        AG,
        OBSERVED,
        "rs1;rs2",
        [region_diag("ploidy_mismatch", "1", 1_000, 1_000, 1)],
        id="call-with-one-character-alias",
    ),
    pytest.param(
        # {(A,G),(G,A)} and {(A,A)} share no pair.
        [row("rs5", "1", 1_000, "AG"), row("rs6", "1", 1_000, "AA")],
        TWO_ANY,
        MISSING,
        "rs5;rs6",
        [locus_diag("position_conflict", "1", 1_000)],
        id="conflicting-calls",
    ),
    pytest.param(
        [
            row("rs1", "1", 1_000, "AG"),
            row("rs2", "1", 1_000, "GA"),
            row("rs3", "1", 1_000, "AA"),
        ],
        TWO_ANY,
        MISSING,
        "rs1;rs2;rs3",
        [locus_diag("position_conflict", "1", 1_000)],
        id="two-agree-one-conflicts",
    ),
    pytest.param(
        [row("rs7", "1", 1_000, None, conflicted=True), row("rs8", "1", 1_000, "AG")],
        TWO_ANY,
        MISSING,
        "rs7;rs8",
        [locus_diag("position_conflict", "1", 1_000)],
        id="flagged-with-called-alias",
    ),
    pytest.param(
        [
            row("rs7", "1", 1_000, None, conflicted=True),
            row("rs8", "1", 1_000, "AG"),
            row("rs9", "1", 1_000, "AA"),
        ],
        TWO_ANY,
        MISSING,
        "rs7;rs8;rs9",
        [locus_diag("position_conflict", "1", 1_000)],
        id="flagged-and-conflicting-one-diagnostic",
    ),
    pytest.param(
        [row("rs7", "1", 1_000, None, conflicted=True)],
        TWO_ANY,
        MISSING,
        "rs7",
        [locus_diag("position_conflict", "1", 1_000)],
        id="flagged-alone",
    ),
]


@pytest.mark.parametrize(
    ("rows", "allowed", "kinds", "rsid", "expected_diagnostics"), SAME_POSITION
)
def test_rows_at_one_position(rows, allowed, kinds, rsid, expected_diagnostics):
    blocks, diagnostics = run_encode([*rows, CONTROL], sex=Sex.FEMALE)
    assert list(blocks) == ["1"]
    assert_block(
        blocks["1"],
        chrom="1",
        pos=[1_000, 2_000],
        rsid=[rsid, "rs99"],
        allowed=[allowed, CC],
        kinds=[kinds, OBSERVED],
    )
    assert diagnostic_rows(diagnostics) == expected_diagnostics


@pytest.mark.parametrize(
    ("chrom", "pos", "names", "expected"),
    [
        pytest.param("1", 1_000, ["1:1000"], None, id="own-placeholder-is-not-a-name"),
        pytest.param("1", 1_000, ["1:999"], "1:999", id="other-position-is-a-name"),
        pytest.param(
            "1",
            1_000,
            ["rs5", "1:1000", "1:999"],
            "1:999;rs5",
            id="placeholder-dropped-names-sorted",
        ),
        pytest.param("1", 1_000, ["rs9", "rs10", "rs9"], "rs10;rs9", id="text-order"),
        pytest.param("2", 1_000, ["1:1000"], "1:1000", id="other-chromosome-is-a-name"),
        pytest.param("X", X_POS, ["X:3000000", "rs4"], "rs4", id="x-placeholder"),
    ],
)
def test_locus_names(chrom, pos, names, expected):
    rows = [row(name, chrom, pos, "AG") for name in names]
    rows.append(row(f"{chrom}:{pos + 1}", chrom, pos + 1, "CC"))
    blocks, diagnostics = run_encode(rows, sex=Sex.FEMALE)
    assert blocks[chrom].loci.pos.tolist() == [pos, pos + 1]
    assert blocks[chrom].loci.rsid.tolist() == [expected, None]
    assert diagnostic_rows(diagnostics) == []


# --------------------------------------------------------------------------------------
# 4. Phase
# --------------------------------------------------------------------------------------


def test_one_explicit_set_is_numbered_zero():
    rows = [
        row("rs1", "1", 100, "AG", phased=True, phase_set=7),
        row("rs2", "1", 200, "CT", phased=True, phase_set=7),
        row("rs3", "1", 300, "GG"),
        row("rs4", "1", 400, "TA", phased=True, phase_set=7),
    ]
    blocks, diagnostics = run_encode(rows, sex=Sex.FEMALE)
    assert_block(
        blocks["1"],
        chrom="1",
        pos=[100, 200, 300, 400],
        rsid=["rs1", "rs2", "rs3", "rs4"],
        # bit 2 (A,G), bit 8 (C,T), bit 12 (G,G), bit 15 (T,A)
        allowed=[AG_ORDERED, CT_ORDERED, GG, TA_ORDERED],
        kinds=[ASSERTED, ASSERTED, OBSERVED, ASSERTED],
        phase_set=[0, 0, -1, 0],
    )
    assert diagnostic_rows(diagnostics) == []


def test_phase_set_none_links_the_phased_rows_of_each_chromosome():
    rows = [
        row("rs1", "1", 100, "AG", phased=True),
        row("rs2", "1", 200, "CT", phased=True),
        row("rs3", "2", 100, "CA", phased=True),
        row("rs4", "2", 200, "GT"),
    ]
    blocks, diagnostics = run_encode(rows, sex=Sex.MALE, omit=("phase_set",))
    assert list(blocks) == ["1", "2"]
    assert_block(
        blocks["1"],
        chrom="1",
        pos=[100, 200],
        rsid=["rs1", "rs2"],
        allowed=[AG_ORDERED, CT_ORDERED],
        kinds=[ASSERTED, ASSERTED],
        phase_set=[0, 0],
    )
    assert_block(
        blocks["2"],
        chrom="2",
        pos=[100, 200],
        rsid=["rs3", "rs4"],
        allowed=[CA_ORDERED, GT],  # bit 5 (C,A); bits 13 and 17
        kinds=[ASSERTED, OBSERVED],
        phase_set=[0, -1],
    )
    assert diagnostic_rows(diagnostics) == []


def test_phased_none_phases_nothing():
    rows = [
        row("rs1", "1", 100, "AG", phased=True, phase_set=3),
        row("rs2", "1", 200, "CT", phased=True, phase_set=3),
    ]
    blocks, diagnostics = run_encode(rows, sex=Sex.FEMALE, omit=("phased",))
    assert_block(
        blocks["1"],
        chrom="1",
        pos=[100, 200],
        rsid=["rs1", "rs2"],
        allowed=[AG, CT],
        kinds=[OBSERVED, OBSERVED],
    )
    assert diagnostic_rows(diagnostics) == []


def test_two_explicit_sets_are_renumbered_by_first_position():
    # Set 9 starts at 100 and becomes 0; set 2 starts at 200 and becomes 1.
    rows = [
        row("rs1", "1", 100, "AG", phased=True, phase_set=9),
        row("rs2", "1", 200, "CT", phased=True, phase_set=2),
        row("rs3", "1", 300, "GA", phased=True, phase_set=9),
        row("rs4", "1", 400, "TC", phased=True, phase_set=2),
    ]
    blocks, diagnostics = run_encode(rows, sex=Sex.FEMALE)
    assert_block(
        blocks["1"],
        chrom="1",
        pos=[100, 200, 300, 400],
        rsid=["rs1", "rs2", "rs3", "rs4"],
        # bits 2 (A,G), 8 (C,T), 10 (G,A), 16 (T,C)
        allowed=[AG_ORDERED, CT_ORDERED, GA_ORDERED, TC_ORDERED],
        kinds=[ASSERTED] * 4,
        phase_set=[0, 1, 0, 1],
    )
    assert diagnostic_rows(diagnostics) == []


def test_na_set_on_phased_row_is_an_implicit_set_apart_from_explicit_zero():
    # The implicit set starts at 100 and becomes 0; explicit set 0 starts at 200 and
    # becomes 1. An unphased row with a set number is not phased.
    rows = [
        row("rs1", "1", 100, "AG", phased=True),
        row("rs2", "1", 200, "CT", phased=True, phase_set=0),
        row("rs3", "1", 300, "GA", phased=True),
        row("rs4", "1", 400, "TC", phased=True, phase_set=0),
        row("rs5", "1", 500, "AT", phase_set=0),
    ]
    blocks, diagnostics = run_encode(rows, sex=Sex.FEMALE)
    assert_block(
        blocks["1"],
        chrom="1",
        pos=[100, 200, 300, 400, 500],
        rsid=["rs1", "rs2", "rs3", "rs4", "rs5"],
        allowed=[AG_ORDERED, CT_ORDERED, GA_ORDERED, TC_ORDERED, AT],
        kinds=[ASSERTED, ASSERTED, ASSERTED, ASSERTED, OBSERVED],
        phase_set=[0, 1, 0, 1, -1],
    )
    assert diagnostic_rows(diagnostics) == []


@pytest.mark.parametrize(
    "aliases",
    [
        pytest.param(
            [
                row("rs1", "1", 1_000, "AG", phased=True, phase_set=3),
                row("rs2", "1", 1_000, "AG", phased=True, phase_set=3),
            ],
            id="same-order-same-set",
        ),
        pytest.param(
            [
                row("rs1", "1", 1_000, "AG", phased=True),
                row("rs2", "1", 1_000, "AG", phased=True),
            ],
            id="same-order-implicit-set",
        ),
        pytest.param(
            # Only phased rows take part; {(A,G)} meets the unordered {(A,G),(G,A)}.
            [
                row("rs1", "1", 1_000, "AG", phased=True, phase_set=3),
                row("rs2", "1", 1_000, "GA"),
            ],
            id="phased-with-unphased-alias",
        ),
        pytest.param(
            [
                row("rs1", "1", 1_000, "AG", phased=True, phase_set=3),
                row("rs2", "1", 1_000, None),
            ],
            id="phased-with-no-call-alias",
        ),
    ],
)
def test_agreeing_phased_rows_keep_their_order(aliases):
    blocks, diagnostics = run_encode([*aliases, CONTROL], sex=Sex.FEMALE)
    assert_block(
        blocks["1"],
        chrom="1",
        pos=[1_000, 2_000],
        rsid=["rs1;rs2", "rs99"],
        allowed=[AG_ORDERED, CC],
        kinds=[ASSERTED, OBSERVED],
        phase_set=[0, -1],
    )
    assert diagnostic_rows(diagnostics) == []


@pytest.mark.parametrize(
    "aliases",
    [
        pytest.param(
            # {(A,G)} and {(G,A)} leave no ordered pair.
            [
                row("rs1", "1", 1_000, "AG", phased=True, phase_set=3),
                row("rs2", "1", 1_000, "GA", phased=True, phase_set=3),
            ],
            id="opposite-order",
        ),
        pytest.param(
            [
                row("rs1", "1", 1_000, "AG", phased=True, phase_set=3),
                row("rs2", "1", 1_000, "AG", phased=True, phase_set=4),
            ],
            id="different-sets",
        ),
        pytest.param(
            [
                row("rs1", "1", 1_000, "AG", phased=True, phase_set=3),
                row("rs2", "1", 1_000, "AG", phased=True),
            ],
            id="explicit-and-implicit-set",
        ),
    ],
)
def test_disagreeing_phased_rows_keep_the_unordered_call(aliases):
    # The locus at 1000 keeps unordered AG and is not in any set; set 3 is first
    # retained at 2000, where it becomes set 0.
    control = row("rs99", "1", 2_000, "CT", phased=True, phase_set=3)
    blocks, diagnostics = run_encode([*aliases, control], sex=Sex.FEMALE)
    assert_block(
        blocks["1"],
        chrom="1",
        pos=[1_000, 2_000],
        rsid=["rs1;rs2", "rs99"],
        allowed=[AG, CT_ORDERED],
        kinds=[OBSERVED, ASSERTED],
        phase_set=[-1, 0],
    )
    assert diagnostic_rows(diagnostics) == [locus_diag("phase_conflict", "1", 1_000)]


def test_set_numbering_skips_a_locus_whose_phase_was_dropped():
    # Set 4 conflicts at 100, so it is first retained at 300 and comes after set 8,
    # first retained at 200.
    rows = [
        row("rs1", "1", 100, "AG", phased=True, phase_set=4),
        row("rs2", "1", 100, "GA", phased=True, phase_set=4),
        row("rs3", "1", 200, "CT", phased=True, phase_set=8),
        row("rs4", "1", 300, "TC", phased=True, phase_set=4),
    ]
    blocks, diagnostics = run_encode(rows, sex=Sex.FEMALE)
    assert_block(
        blocks["1"],
        chrom="1",
        pos=[100, 200, 300],
        rsid=["rs1;rs2", "rs3", "rs4"],
        allowed=[AG, CT_ORDERED, TC_ORDERED],
        kinds=[OBSERVED, ASSERTED, ASSERTED],
        phase_set=[-1, 0, 1],
    )
    assert diagnostic_rows(diagnostics) == [locus_diag("phase_conflict", "1", 100)]


@pytest.mark.parametrize(
    "aliases",
    [
        pytest.param(
            [
                row("rs1", "1", 1_000, None, conflicted=True),
                row("rs2", "1", 1_000, "AG", phased=True, phase_set=1),
            ],
            id="flagged-with-phased-alias",
        ),
        pytest.param(
            [
                row("rs1", "1", 1_000, "AG", phased=True, phase_set=1),
                row("rs2", "1", 1_000, "AA", phased=True, phase_set=1),
            ],
            id="conflicting-phased-calls",
        ),
    ],
)
def test_position_conflict_is_settled_before_phase(aliases):
    blocks, diagnostics = run_encode([*aliases, CONTROL], sex=Sex.FEMALE)
    assert_block(
        blocks["1"],
        chrom="1",
        pos=[1_000, 2_000],
        rsid=["rs1;rs2", "rs99"],
        allowed=[TWO_ANY, CC],
        kinds=[MISSING, OBSERVED],
    )
    assert diagnostic_rows(diagnostics) == [locus_diag("position_conflict", "1", 1_000)]


@pytest.mark.parametrize(
    ("chrom", "pos", "sex", "genotype", "phased", "allowed", "expected_diagnostics"),
    [
        pytest.param("X", X_POS, Sex.MALE, "AA", True, A_ONE, [], id="one-copy-AA"),
        pytest.param("X", X_POS, Sex.MALE, "G", True, G_ONE, [], id="one-copy-G"),
        pytest.param(
            "X",
            X_POS,
            Sex.MALE,
            "AG",
            True,
            ONE_ANY,
            [region_diag("ploidy_mismatch", "X", X_POS, X_POS, 1)],
            id="one-copy-two-bases",
        ),
        pytest.param(
            "1",
            1_000,
            Sex.FEMALE,
            "A",
            True,
            TWO_ANY,
            [region_diag("ploidy_mismatch", "1", 1_000, 1_000, 1)],
            id="two-copy-one-base",
        ),
        pytest.param(
            "1", 1_000, Sex.FEMALE, None, True, TWO_ANY, [], id="two-copy-no-call"
        ),
        pytest.param(
            "1",
            1_000,
            Sex.FEMALE,
            "DI",
            True,
            TWO_ANY,
            [region_diag("non_base_genotype", "1", 1_000, 1_000, 1)],
            id="two-copy-non-base",
        ),
    ],
)
def test_phase_flags_are_ignored_off_two_base_two_copy_rows(
    chrom, pos, sex, genotype, phased, allowed, expected_diagnostics
):
    rows = [row("rs1", chrom, pos, genotype, phased=phased, phase_set=5)]
    blocks, diagnostics = run_encode(rows, sex=sex)
    called = allowed not in (TWO_ANY, ONE_ANY)
    assert_block(
        blocks[chrom],
        chrom=chrom,
        pos=[pos],
        rsid=["rs1"],
        allowed=[allowed],
        kinds=[OBSERVED if called else MISSING],
        phase_set=[-1],
    )
    assert diagnostic_rows(diagnostics) == expected_diagnostics


def test_phased_homozygous_call_is_an_asserted_two_base_call():
    # "GG" phased is the ordered domain {(G,G)}, bit 12, kept with kind asserted.
    rows = [row("rs1", "1", 100, "GG", phased=True, phase_set=5)]
    blocks, diagnostics = run_encode(rows, sex=Sex.FEMALE)
    assert_block(
        blocks["1"],
        chrom="1",
        pos=[100],
        rsid=["rs1"],
        allowed=[GG],
        kinds=[ASSERTED],
        phase_set=[0],
    )
    assert diagnostic_rows(diagnostics) == []


# --------------------------------------------------------------------------------------
# 5. X boundaries, excluded chromosomes and invalid positions
# --------------------------------------------------------------------------------------

# GRCh37 PAR1 is 60,001-2,699,520 and PAR2 is 154,931,044-155,260,560, inclusive.
X_BOUNDARIES = [
    60_000,  # before PAR1
    60_001,  # PAR1 start
    2_699_520,  # PAR1 end
    2_699_521,  # after PAR1
    154_931_043,  # before PAR2
    154_931_044,  # PAR2 start
    155_260_560,  # PAR2 end
    155_260_561,  # after PAR2
]
X_KEPT = [60_000, 2_699_521, 154_931_043, 155_260_561]
PAR_DIAGNOSTICS = [
    region_diag("excluded_region", "X", 60_001, 2_699_520, 2),
    region_diag("excluded_region", "X", 154_931_044, 155_260_560, 2),
]


@pytest.mark.parametrize(
    ("sex", "genotypes", "allowed"),
    [
        pytest.param(
            Sex.MALE,
            # PAR rows hold two different bases: excluded rows get no ploidy check.
            ["A", "AG", "AG", "C", "GG", "AG", "CT", "T"],
            [A_ONE, C_ONE, G_ONE, T_ONE],  # bits 4, 9, 14, 19
            id="male",
        ),
        pytest.param(
            Sex.FEMALE,
            # PAR rows hold one character: excluded rows get no ploidy check.
            ["AG", "A", "A", "CT", "GG", "T", "C", "TA"],
            [AG, CT, GG, AT],  # bits 2+10, 8+16, 12, 3+15
            id="female",
        ),
    ],
)
def test_x_pseudoautosomal_boundaries(sex, genotypes, allowed):
    rows = [
        row(f"rs{i + 1}", "X", pos, genotype)
        for i, (pos, genotype) in enumerate(zip(X_BOUNDARIES, genotypes, strict=True))
    ]
    blocks, diagnostics = run_encode(rows, sex=sex)
    assert list(blocks) == ["X"]
    assert_block(
        blocks["X"],
        chrom="X",
        pos=X_KEPT,
        rsid=["rs1", "rs4", "rs5", "rs8"],
        allowed=allowed,
        kinds=[OBSERVED] * 4,
    )
    assert diagnostic_rows(diagnostics) == PAR_DIAGNOSTICS


def test_x_for_unknown_sex_is_excluded_in_one_row():
    rows = [
        row("rs1", "1", 100, "AG"),
        row("rs2", "X", 60_000, "AG"),
        row("rs3", "X", 100_000, "A"),  # in PAR1
        row("rs4", "X", 5_000_000, "A"),
    ]
    blocks, diagnostics = run_encode(rows, sex=None)
    assert list(blocks) == ["1"]
    assert diagnostic_rows(diagnostics) == [
        region_diag("excluded_region", "X", 60_000, 5_000_000, 3)
    ]


def test_y_mt_and_other_names_are_excluded():
    rows = [
        row("rs1", "1", 100, "AG"),
        row("rs2", "Y", 2_000_000, "A"),
        row("rs3", "Y", 3_000_000, "AG"),
        row("rs4", "MT", 100, "II"),  # excluded rows get no genotype check
        row("rs5", "0", 7, "AA"),
        row("rs6", "PAR", 500, "AG"),
    ]
    for sex in (Sex.MALE, Sex.FEMALE, None):
        blocks, diagnostics = run_encode(rows, sex=sex)
        assert list(blocks) == ["1"]
        assert_block(
            blocks["1"],
            chrom="1",
            pos=[100],
            rsid=["rs1"],
            allowed=[AG],
            kinds=[OBSERVED],
        )
        # Y, MT, then the other names in sorted order: "0" before "PAR".
        assert diagnostic_rows(diagnostics) == [
            region_diag("excluded_region", "Y", 2_000_000, 3_000_000, 2),
            region_diag("excluded_region", "MT", 100, 100, 1),
            region_diag("excluded_region", "0", 7, 7, 1),
            region_diag("excluded_region", "PAR", 500, 500, 1),
        ]


@pytest.mark.parametrize(
    ("rows", "sex", "kept", "expected_diagnostics"),
    [
        pytest.param(
            # chr1 is 249,250,621 long; two rows at one invalid position, one report.
            [
                row("rs1", "1", 249_250_621, "AG"),
                row("rs2", "1", 249_250_622, "CT"),
                row("rs3", "1", 249_250_622, "CT"),
            ],
            Sex.FEMALE,
            {"1": ([249_250_621], ["rs1"], [AG])},
            [locus_diag("invalid_position", "1", 249_250_622)],
            id="chr1-end",
        ),
        pytest.param(
            # X is 155,270,560 long.
            [row("rs1", "X", 155_270_560, "A"), row("rs2", "X", 155_270_561, "C")],
            Sex.MALE,
            {"X": ([155_270_560], ["rs1"], [A_ONE])},
            [locus_diag("invalid_position", "X", 155_270_561)],
            id="x-end",
        ),
        pytest.param(
            # A block is never empty: a chromosome with only invalid rows has none.
            [row("rs1", "2", 1_000, "CC"), row("rs2", "1", 249_250_622, "CT")],
            Sex.FEMALE,
            {"2": ([1_000], ["rs1"], [CC])},
            [locus_diag("invalid_position", "1", 249_250_622)],
            id="only-invalid-rows",
        ),
    ],
)
def test_invalid_positions_stay_off_the_axis(rows, sex, kept, expected_diagnostics):
    blocks, diagnostics = run_encode(rows, sex=sex)
    assert list(blocks) == list(kept)
    for chrom, (pos, rsid, allowed) in kept.items():
        assert_block(
            blocks[chrom],
            chrom=chrom,
            pos=pos,
            rsid=rsid,
            allowed=allowed,
            kinds=[OBSERVED] * len(pos),
        )
    assert diagnostic_rows(diagnostics) == expected_diagnostics


# --------------------------------------------------------------------------------------
# 6. The diagnostics frame
# --------------------------------------------------------------------------------------


def test_diagnostics_frame_schema_and_order():
    rows = [
        row("rs1", "1", 100, "AG"),
        row("rs2", "1", 150, "AG"),
        row("rs3", "1", 150, "AA"),  # conflicts with rs2
        row("rs4", "1", 200, "II"),
        row("rs5", "1", 300, "A"),
        row("rs6", "1", 400, "DI"),
        row("rs7", "1", 500, "AG", phased=True, phase_set=1),
        row("rs8", "1", 500, "GA", phased=True, phase_set=1),  # opposite order
        row("rs9", "1", 249_250_622, "AG"),
        row("rs10", "2", 50, "A"),
        row("rs11", "2", 50, "I"),
        row("rs12", "2", 60, "CC"),
        row("rs13", "X", 100_000, "AG"),  # PAR1
        row("rs14", "X", X_POS, "AG"),
        row("rs15", "Y", 10, "A"),
        row("rs16", "Y", 20, "C"),
        row("rs17", "MT", 5, "A"),
        row("rs18", "0", 7, "AA"),
        row("rs19", "PAR", 1, "AG"),
    ]
    blocks, diagnostics = run_encode(rows, sex=Sex.FEMALE)

    assert list(diagnostics.columns) == [
        "code",
        "chrom",
        "pos",
        "start",
        "end",
        "detail",
    ]
    assert diagnostics["code"].dtype == "str"
    assert diagnostics["chrom"].dtype == "str"
    assert diagnostics["detail"].dtype == "str"
    assert diagnostics["pos"].dtype == "Int64"
    assert diagnostics["start"].dtype == "Int64"
    assert diagnostics["end"].dtype == "Int64"
    assert isinstance(diagnostics.index, pd.RangeIndex)

    # Ordered by chromosome, then pos or start, then code. Aggregated codes carry
    # start/end and a leading row count; per-locus codes carry pos only.
    assert diagnostic_rows(diagnostics) == [
        locus_diag("position_conflict", "1", 150),
        region_diag("non_base_genotype", "1", 200, 400, 2),
        region_diag("ploidy_mismatch", "1", 300, 300, 1),
        locus_diag("phase_conflict", "1", 500),
        locus_diag("invalid_position", "1", 249_250_622),
        region_diag("non_base_genotype", "2", 50, 50, 1),
        region_diag("ploidy_mismatch", "2", 50, 50, 1),
        region_diag("excluded_region", "X", 100_000, 100_000, 1),
        region_diag("excluded_region", "Y", 10, 20, 2),
        region_diag("excluded_region", "MT", 5, 5, 1),
        region_diag("excluded_region", "0", 7, 7, 1),
        region_diag("excluded_region", "PAR", 1, 1, 1),
    ]
    assert diagnostics["detail"].notna().all()

    assert list(blocks) == ["1", "2", "X"]
    assert_block(
        blocks["1"],
        chrom="1",
        pos=[100, 150, 200, 300, 400, 500],
        rsid=["rs1", "rs2;rs3", "rs4", "rs5", "rs6", "rs7;rs8"],
        allowed=[AG, TWO_ANY, TWO_ANY, TWO_ANY, TWO_ANY, AG],
        kinds=[OBSERVED, MISSING, MISSING, MISSING, MISSING, OBSERVED],
    )
    assert_block(
        blocks["2"],
        chrom="2",
        pos=[50, 60],
        rsid=["rs10;rs11", "rs12"],
        allowed=[TWO_ANY, CC],
        kinds=[MISSING, OBSERVED],
    )
    assert_block(
        blocks["X"],
        chrom="X",
        pos=[X_POS],
        rsid=["rs14"],
        allowed=[AG],
        kinds=[OBSERVED],
    )


@pytest.mark.parametrize("sex", [Sex.MALE, Sex.FEMALE, None])
def test_empty_frame_gives_empty_results_with_full_schema(sex):
    blocks, diagnostics = encode(make_frame([]), sex=sex)
    assert isinstance(blocks, Mapping)
    assert len(blocks) == 0
    assert len(diagnostics) == 0
    assert_diagnostics_schema(diagnostics)


def test_no_diagnostics_still_has_full_schema():
    _, diagnostics = run_encode([row("rs1", "1", 100, "AG")], sex=Sex.FEMALE)
    assert len(diagnostics) == 0
    assert_diagnostics_schema(diagnostics)


# --------------------------------------------------------------------------------------
# encode input rejections
# --------------------------------------------------------------------------------------


def _one_row_frame() -> pd.DataFrame:
    return make_frame([row("rs1", "1", 100, "AG")])


def _with_pos(values, dtype):
    def build():
        frame = make_frame([row("rs1", "1", 100, "AG")] * len(values))
        frame["pos"] = np.array(values, dtype=dtype)
        return frame, {}

    return build


def _integer_phased():
    frame = _one_row_frame()
    return frame, {
        "phased": pd.Series(np.array([1], dtype=np.int64), index=frame.index.copy())
    }


def _text_conflicted():
    frame = _one_row_frame()
    flags = pd.Series(["yes"], dtype=object, index=frame.index.copy())
    return frame, {"genotype_conflicted": flags}


ENCODE_REJECTIONS = [
    pytest.param(
        lambda: (_one_row_frame().drop(columns="genotype"), {}), id="no-genotype"
    ),
    pytest.param(lambda: (_one_row_frame().drop(columns="pos"), {}), id="no-pos"),
    pytest.param(lambda: (_one_row_frame().drop(columns="chrom"), {}), id="no-chrom"),
    pytest.param(
        lambda: (make_frame([row("rs1", None, 100, "AG")]), {}), id="na-chrom"
    ),
    pytest.param(lambda: (make_frame([row("rs1", "1", 0, "AG")]), {}), id="zero-pos"),
    pytest.param(
        lambda: (make_frame([row("rs1", "1", 100, "AG"), row("rs2", "Y", 0, "A")]), {}),
        id="zero-pos-on-excluded-row",
    ),
    pytest.param(_with_pos([-5], np.int64), id="negative-pos"),
    pytest.param(_with_pos([2**32], np.int64), id="pos-beyond-uint32"),
    pytest.param(_with_pos([100.5], np.float64), id="fractional-pos"),
    pytest.param(_with_pos([np.nan], np.float64), id="nan-pos"),
    pytest.param(_with_pos([np.inf], np.float64), id="infinite-pos"),
    pytest.param(_integer_phased, id="integer-phased"),
    pytest.param(_text_conflicted, id="text-conflicted"),
]


@pytest.mark.parametrize("build", ENCODE_REJECTIONS)
def test_encode_rejects_malformed_input(build):
    frame, series = build()
    with pytest.raises(ValueError):
        encode(frame, sex=Sex.FEMALE, **series)


# --------------------------------------------------------------------------------------
# 7. _Calls construction
# --------------------------------------------------------------------------------------


def autosome_args() -> dict:
    return {
        "loci": Loci("1", [100, 200]),
        "allowed": [AG, TWO_ANY],
        "phase_set": [-1, -1],
        "kinds": [[1, 0, 0], [0, 0, 0]],
    }


def male_x_args() -> dict:
    return {
        "loci": Loci("X", [3_000_000, 4_000_000]),
        "allowed": [A_ONE, ONE_ANY],
        "phase_set": [-1, -1],
        "kinds": [[1, 0, 0], [0, 0, 0]],
    }


def female_x_args() -> dict:
    return {**autosome_args(), "loci": Loci("X", [3_000_000, 4_000_000])}


VALID_CALLS = [
    pytest.param(autosome_args, {}, id="autosome"),
    pytest.param(male_x_args, {}, id="male-x"),
    pytest.param(female_x_args, {}, id="female-x"),
    pytest.param(
        male_x_args,
        {"kinds": [[1, 0, 2], [0, 0, 0]]},
        id="male-x-forced-origin",
    ),
    pytest.param(
        autosome_args,
        {"allowed": [A_ON_ONE_COPY, TWO_ANY], "kinds": [[2, 0, 0], [0, 0, 0]]},
        id="forced-partial-genotype",
    ),
    pytest.param(
        autosome_args,
        {
            "allowed": [AG_ORDERED, CT_ORDERED],
            "phase_set": [0, 0],
            "kinds": [[1, 2, 2], [1, 2, 2]],
        },
        id="forced-phase-and-origin",
    ),
    pytest.param(
        autosome_args,
        {
            "allowed": [AG_ORDERED, CT_ORDERED],
            "phase_set": [0, 3],
            "kinds": [[1, 4, 0], [1, 4, 0]],
        },
        id="asserted-phase",
    ),
    pytest.param(
        autosome_args,
        {
            "allowed": [AG_ORDERED, CT_ORDERED],
            "phase_set": [2, 2],
            "kinds": [[3, 3, 3], [3, 3, 3]],
        },
        id="estimated",
    ),
    pytest.param(
        autosome_args,
        {
            "allowed": np.array([AG, TWO_ANY], dtype=np.uint64),
            "phase_set": np.array([-1, -1], dtype=np.int8),
            "kinds": np.array([[1, 0, 0], [0, 0, 0]], dtype=np.uint8),
        },
        id="other-integer-dtypes",
    ),
]


@pytest.mark.parametrize(("base", "changes"), VALID_CALLS)
def test_calls_accepts_valid_blocks(base, changes):
    args = {**base(), **changes}
    calls = make_calls(**args)
    assert calls.loci is args["loci"]
    assert calls.allowed.dtype == np.uint32
    assert calls.phase_set.dtype == np.int32
    assert calls.kinds.dtype == np.uint8
    for name in ("allowed", "phase_set", "kinds"):
        assert getattr(calls, name).tolist() == np.asarray(args[name]).tolist(), name


REJECTED_CALLS = [
    # TypeError: not a Loci, or an array without an integer dtype.
    pytest.param(
        autosome_args, {"loci": ("1", [100, 200])}, TypeError, id="loci-not-loci"
    ),
    pytest.param(
        autosome_args,
        {"allowed": np.array([AG, TWO_ANY], dtype=np.float64)},
        TypeError,
        id="allowed-float",
    ),
    pytest.param(
        autosome_args, {"allowed": np.array([True, True])}, TypeError, id="allowed-bool"
    ),
    pytest.param(
        autosome_args,
        {"allowed": np.array([AG, TWO_ANY], dtype=object)},
        TypeError,
        id="allowed-object",
    ),
    pytest.param(
        autosome_args,
        {"phase_set": np.array([-1.0, -1.0])},
        TypeError,
        id="phase-set-float",
    ),
    pytest.param(
        autosome_args,
        {"kinds": np.array([[1.0, 0, 0], [0, 0, 0]])},
        TypeError,
        id="kinds-float",
    ),
    pytest.param(
        autosome_args,
        {"kinds": np.array([[True, False, False], [False] * 3])},
        TypeError,
        id="kinds-bool",
    ),
    # ValueError: shapes.
    pytest.param(
        autosome_args, {"allowed": [[AG], [TWO_ANY]]}, ValueError, id="allowed-2d"
    ),
    pytest.param(
        autosome_args,
        {"allowed": [AG, TWO_ANY, TWO_ANY]},
        ValueError,
        id="allowed-too-long",
    ),
    pytest.param(
        autosome_args, {"phase_set": [-1]}, ValueError, id="phase-set-too-short"
    ),
    pytest.param(
        autosome_args, {"kinds": [[1, 0], [0, 0]]}, ValueError, id="kinds-n-by-2"
    ),
    pytest.param(
        autosome_args, {"kinds": [1, 0, 0, 0, 0, 0]}, ValueError, id="kinds-flat"
    ),
    pytest.param(
        autosome_args,
        {"kinds": [[1, 0, 0], [0, 0, 0], [0, 0, 0]]},
        ValueError,
        id="kinds-too-long",
    ),
    # ValueError: masks.
    pytest.param(autosome_args, {"allowed": [0, TWO_ANY]}, ValueError, id="zero-mask"),
    pytest.param(
        autosome_args, {"allowed": [AG | 1 << 25, TWO_ANY]}, ValueError, id="bit-25"
    ),
    pytest.param(
        autosome_args, {"allowed": [AG | 1 << 31, TWO_ANY]}, ValueError, id="bit-31"
    ),
    # 1028 | 1 << 32 would narrow to a valid uint32 1028.
    pytest.param(
        autosome_args,
        {"allowed": [AG | 1 << 32, TWO_ANY]},
        ValueError,
        id="bit-32-before-narrowing",
    ),
    pytest.param(
        autosome_args, {"allowed": [-1, TWO_ANY]}, ValueError, id="negative-mask"
    ),
    pytest.param(
        male_x_args,
        {"allowed": [1 << 20, ONE_ANY]},
        ValueError,
        id="absent-copy-0-bit-20",
    ),
    pytest.param(
        male_x_args,
        {"allowed": [1 << 24, ONE_ANY]},
        ValueError,
        id="absent-both-bit-24",
    ),
    # bit 0 (A,A) with bit 4 (A, ABSENT).
    pytest.param(
        male_x_args,
        {"allowed": [AA | A_ONE, ONE_ANY]},
        ValueError,
        id="mixed-ploidy-mask-on-x",
    ),
    pytest.param(
        autosome_args,
        {"allowed": [AG | A_ONE, TWO_ANY]},
        ValueError,
        id="mixed-ploidy-mask-on-autosome",
    ),
    pytest.param(
        male_x_args,
        {"allowed": [AG, ONE_ANY]},
        ValueError,
        id="rows-of-different-ploidy",
    ),
    pytest.param(
        autosome_args,
        {"allowed": [A_ONE, ONE_ANY]},
        ValueError,
        id="one-copy-on-autosome",
    ),
    # ValueError: phase sets.
    pytest.param(
        autosome_args, {"phase_set": [-2, -1]}, ValueError, id="phase-set-below-minus-1"
    ),
    # 2**32 would narrow to int32 0.
    pytest.param(
        autosome_args,
        {"phase_set": [2**32, -1], "kinds": [[1, 4, 0], [0, 0, 0]]},
        ValueError,
        id="phase-set-before-narrowing",
    ),
    pytest.param(
        male_x_args,
        {"phase_set": [0, -1], "kinds": [[1, 2, 0], [0, 0, 0]]},
        ValueError,
        id="one-copy-with-phase",
    ),
    pytest.param(
        autosome_args,
        {
            "allowed": [AG_ORDERED, CT_ORDERED],
            "phase_set": [0, 0],
            "kinds": [[1, 2, 2], [1, 2, 0]],
        },
        ValueError,
        id="origin-known-at-part-of-a-set",
    ),
    # ValueError: kind codes outside their column's set.
    pytest.param(
        autosome_args,
        {"kinds": [[4, 0, 0], [0, 0, 0]]},
        ValueError,
        id="genotype-kind-asserted",
    ),
    # 257 would narrow to uint8 1.
    pytest.param(
        autosome_args,
        {"kinds": [[257, 0, 0], [0, 0, 0]]},
        ValueError,
        id="genotype-kind-before-narrowing",
    ),
    pytest.param(
        autosome_args,
        {"phase_set": [0, -1], "kinds": [[1, 1, 0], [0, 0, 0]]},
        ValueError,
        id="phase-kind-observed",
    ),
    pytest.param(
        autosome_args,
        {
            "phase_set": [0, -1],
            "kinds": [[1, 2, 1], [0, 0, 0]],
        },
        ValueError,
        id="origin-kind-observed",
    ),
    pytest.param(
        autosome_args,
        {
            "phase_set": [0, -1],
            "kinds": [[1, 2, 4], [0, 0, 0]],
        },
        ValueError,
        id="origin-kind-asserted",
    ),
    # ValueError: a kind is 0 exactly when its attribute is missing.
    pytest.param(
        autosome_args,
        {"kinds": [[0, 0, 0], [0, 0, 0]]},
        ValueError,
        id="called-genotype-kind-0",
    ),
    pytest.param(
        autosome_args,
        {"kinds": [[1, 0, 0], [1, 0, 0]]},
        ValueError,
        id="unrestricted-two-copy-genotype-kind-1",
    ),
    pytest.param(
        male_x_args,
        {"kinds": [[1, 0, 0], [2, 0, 0]]},
        ValueError,
        id="unrestricted-one-copy-genotype-kind-2",
    ),
    pytest.param(
        autosome_args, {"phase_set": [0, -1]}, ValueError, id="phase-set-with-kind-0"
    ),
    pytest.param(
        autosome_args,
        {"kinds": [[1, 4, 0], [0, 0, 0]]},
        ValueError,
        id="phase-kind-without-set",
    ),
]


@pytest.mark.parametrize(("base", "changes", "error"), REJECTED_CALLS)
def test_calls_rejects(base, changes, error):
    args = {**base(), **changes}
    with pytest.raises(error):
        make_calls(**args)


def test_calls_rejects_an_empty_block():
    with pytest.raises(ValueError):
        _Calls(
            Loci("1", []),
            np.zeros(0, dtype=np.int64),
            np.zeros(0, dtype=np.int64),
            np.zeros((0, 3), dtype=np.int64),
        )


# --------------------------------------------------------------------------------------
# 8. align
# --------------------------------------------------------------------------------------


def test_align_builds_one_union_axis_per_chromosome():
    p1 = {
        "10": block("10", [700], ["rs7"], [GG], [OBSERVED]),
        "1": block(
            "1",
            [100, 300, 500],
            ["rs1", "rs3;rs9", None],
            [AG, AG_ORDERED, TWO_ANY],
            [OBSERVED, ASSERTED, MISSING],
            phase_set=[-1, 0, -1],
        ),
    }
    p2 = {
        "1": block(
            "1",
            [200, 300, 600],
            ["rs2", "rs3;rs30", None],
            [CT, GG, CC],
            [OBSERVED] * 3,
        ),
        "2": block("2", [50], [None], [AA], [OBSERVED]),
    }
    out = align({"p1": p1, "p2": p2}, sexes={"p2": Sex.FEMALE, "p1": Sex.MALE})

    assert list(out) == ["p2", "p1"]
    for name in out:
        assert list(out[name]) == ["1", "2", "10"]
        for chrom in out[name]:
            assert out[name][chrom].loci is out["p1"][chrom].loci

    # Axis for chr1: 100, 200, 300, 500, 600. At 300 the names "rs3;rs9" and
    # "rs3;rs30" merge into the distinct names rs3, rs30, rs9 in text order.
    axis = ([100, 200, 300, 500, 600], ["rs1", "rs2", "rs3;rs30;rs9", None, None])
    assert_block(
        out["p1"]["1"],
        chrom="1",
        pos=axis[0],
        rsid=axis[1],
        allowed=[AG, TWO_ANY, AG_ORDERED, TWO_ANY, TWO_ANY],
        kinds=[OBSERVED, MISSING, ASSERTED, MISSING, MISSING],
        phase_set=[-1, -1, 0, -1, -1],
    )
    assert_block(
        out["p2"]["1"],
        chrom="1",
        pos=axis[0],
        rsid=axis[1],
        allowed=[TWO_ANY, CT, GG, TWO_ANY, CC],
        kinds=[MISSING, OBSERVED, OBSERVED, MISSING, OBSERVED],
    )
    assert_block(
        out["p1"]["2"],
        chrom="2",
        pos=[50],
        rsid=[None],
        allowed=[TWO_ANY],
        kinds=[MISSING],
    )
    assert_block(
        out["p2"]["2"], chrom="2", pos=[50], rsid=[None], allowed=[AA], kinds=[OBSERVED]
    )
    assert_block(
        out["p1"]["10"],
        chrom="10",
        pos=[700],
        rsid=["rs7"],
        allowed=[GG],
        kinds=[OBSERVED],
    )
    assert_block(
        out["p2"]["10"],
        chrom="10",
        pos=[700],
        rsid=["rs7"],
        allowed=[TWO_ANY],
        kinds=[MISSING],
    )


def test_align_gives_an_untested_participant_unrestricted_blocks():
    tested = {
        "1": block("1", [100], ["rs1"], [AG], [OBSERVED]),
        "X": block("X", [X_POS], ["rs2"], [CT], [OBSERVED]),
    }
    out = align(
        {"mother": tested},
        sexes={"mother": Sex.FEMALE, "son": Sex.MALE, "daughter": Sex.FEMALE},
    )
    assert list(out) == ["mother", "son", "daughter"]
    assert list(out["son"]) == list(out["daughter"]) == ["1", "X"]
    assert_block(
        out["son"]["1"],
        chrom="1",
        pos=[100],
        rsid=["rs1"],
        allowed=[TWO_ANY],
        kinds=[MISSING],
    )
    assert_block(
        out["son"]["X"],
        chrom="X",
        pos=[X_POS],
        rsid=["rs2"],
        allowed=[ONE_ANY],
        kinds=[MISSING],
    )
    assert_block(
        out["daughter"]["X"],
        chrom="X",
        pos=[X_POS],
        rsid=["rs2"],
        allowed=[TWO_ANY],
        kinds=[MISSING],
    )
    assert_block(
        out["mother"]["X"],
        chrom="X",
        pos=[X_POS],
        rsid=["rs2"],
        allowed=[CT],
        kinds=[OBSERVED],
    )
    assert out["son"]["X"].loci is out["mother"]["X"].loci


def test_align_gives_no_x_block_for_unknown_sex():
    p = {
        "1": block("1", [100], ["rs1"], [AG], [OBSERVED]),
        "X": block("X", [X_POS], ["rs2"], [CT], [OBSERVED]),
    }
    q = {"1": block("1", [200], ["rs5"], [CC], [OBSERVED])}
    out = align({"p": p, "q": q}, sexes={"p": Sex.FEMALE, "q": None, "r": None})
    assert list(out) == ["p", "q", "r"]
    assert list(out["p"]) == ["1", "X"]
    assert list(out["q"]) == ["1"]
    assert list(out["r"]) == ["1"]
    assert_block(
        out["q"]["1"],
        chrom="1",
        pos=[100, 200],
        rsid=["rs1", "rs5"],
        allowed=[TWO_ANY, CC],
        kinds=[MISSING, OBSERVED],
    )
    assert_block(
        out["r"]["1"],
        chrom="1",
        pos=[100, 200],
        rsid=["rs1", "rs5"],
        allowed=[TWO_ANY, TWO_ANY],
        kinds=[MISSING, MISSING],
    )
    assert_block(
        out["p"]["X"],
        chrom="X",
        pos=[X_POS],
        rsid=["rs2"],
        allowed=[CT],
        kinds=[OBSERVED],
    )


def test_align_fills_x_by_ploidy():
    male = {
        "X": block(
            "X",
            [3_000_000, 5_000_000],
            ["rs1", "rs3"],
            [A_ONE, G_ONE],
            [OBSERVED, OBSERVED],
        )
    }
    female = {
        "X": block(
            "X", [4_000_000, 5_000_000], ["rs2", "rs3"], [AG, GG], [OBSERVED, OBSERVED]
        )
    }
    out = align({"m": male, "f": female}, sexes={"m": Sex.MALE, "f": Sex.FEMALE})
    axis = [3_000_000, 4_000_000, 5_000_000]
    names = ["rs1", "rs2", "rs3"]
    assert_block(
        out["m"]["X"],
        chrom="X",
        pos=axis,
        rsid=names,
        allowed=[A_ONE, ONE_ANY, G_ONE],
        kinds=[OBSERVED, MISSING, OBSERVED],
    )
    assert_block(
        out["f"]["X"],
        chrom="X",
        pos=axis,
        rsid=names,
        allowed=[TWO_ANY, AG, GG],
        kinds=[MISSING, OBSERVED, OBSERVED],
    )
    assert out["m"]["X"].loci is out["f"]["X"].loci


def test_align_preserves_a_partial_inferred_block():
    child = {
        "1": block(
            "1",
            [100, 200, 300],
            ["rs1", "rs2", "rs3"],
            [COPY0_A, COPY0_G, TWO_ANY],
            [(2, 2, 2), (3, 2, 2), MISSING],
            phase_set=[0, 0, -1],
        )
    }
    son = {"X": block("X", [3_000_000], ["rs40"], [A_ONE], [(3, 0, 2)])}
    mother = {
        "1": block("1", [150, 300], ["rs15", "rs3"], [AG, CC], [OBSERVED, OBSERVED]),
        "X": block("X", [4_000_000], ["rs41"], [AG], [OBSERVED]),
    }
    out = align(
        {"child": child, "son": son, "mother": mother},
        sexes={"mother": Sex.FEMALE, "child": Sex.FEMALE, "son": Sex.MALE},
    )
    assert list(out) == ["mother", "child", "son"]
    assert_block(
        out["child"]["1"],
        chrom="1",
        pos=[100, 150, 200, 300],
        rsid=["rs1", "rs15", "rs2", "rs3"],
        allowed=[COPY0_A, TWO_ANY, COPY0_G, TWO_ANY],
        kinds=[(2, 2, 2), MISSING, (3, 2, 2), MISSING],
        phase_set=[0, -1, 0, -1],
    )
    assert_block(
        out["son"]["X"],
        chrom="X",
        pos=[3_000_000, 4_000_000],
        rsid=["rs40", "rs41"],
        allowed=[A_ONE, ONE_ANY],
        kinds=[(3, 0, 2), MISSING],
    )
    assert_block(
        out["child"]["X"],
        chrom="X",
        pos=[3_000_000, 4_000_000],
        rsid=["rs40", "rs41"],
        allowed=[TWO_ANY, TWO_ANY],
        kinds=[MISSING, MISSING],
    )
    assert_block(
        out["son"]["1"],
        chrom="1",
        pos=[100, 150, 200, 300],
        rsid=["rs1", "rs15", "rs2", "rs3"],
        allowed=[TWO_ANY] * 4,
        kinds=[MISSING] * 4,
    )


def _autosome_block() -> _Calls:
    return block("1", [100], ["rs1"], [AG], [OBSERVED])


def _female_x_block() -> _Calls:
    return block("X", [X_POS], ["rs2"], [AG], [OBSERVED])


def _male_x_block() -> _Calls:
    return block("X", [X_POS], ["rs2"], [A_ONE], [OBSERVED])


ALIGN_REJECTIONS = [
    pytest.param(
        lambda: {"p": {"1": _autosome_block()}, "q": {"1": _autosome_block()}},
        {"p": Sex.FEMALE},
        id="person-missing-from-sexes",
    ),
    pytest.param(
        lambda: {"p": {"2": _autosome_block()}},
        {"p": Sex.FEMALE},
        id="key-differs-from-chromosome",
    ),
    pytest.param(
        lambda: {"p": {"X": _female_x_block()}},
        {"p": None},
        id="x-block-for-unknown-sex",
    ),
    pytest.param(
        lambda: {"p": {"X": _female_x_block()}},
        {"p": Sex.MALE},
        id="male-with-two-copy-x",
    ),
    pytest.param(
        lambda: {"p": {"X": _male_x_block()}},
        {"p": Sex.FEMALE},
        id="female-with-one-copy-x",
    ),
]


@pytest.mark.parametrize(("people", "sexes"), ALIGN_REJECTIONS)
def test_align_rejects(people, sexes):
    with pytest.raises(ValueError):
        align(people(), sexes=sexes)


# --------------------------------------------------------------------------------------
# 9. Properties
# --------------------------------------------------------------------------------------

ALL_SHORT_GENOTYPES = [
    "".join(chars) for n in (1, 2) for chars in itertools.product("ACGTID", repeat=n)
]


def reference_domain(genotype: str, copies: int) -> tuple[int | None, str | None]:
    """The per-row genotype rule stated independently: (mask or None, diagnostic code)."""
    if not 1 <= len(genotype) <= 2 or any(ch not in "ACGT" for ch in genotype):
        return None, "non_base_genotype"
    codes = ["ACGT".index(ch) for ch in genotype]
    if copies == 1:
        if len(set(codes)) != 1:
            return None, "ploidy_mismatch"
        pairs = {(codes[0], NO_COPY)}
    else:
        if len(codes) != 2:
            return None, "ploidy_mismatch"
        pairs = {(codes[0], codes[1]), (codes[1], codes[0])}
    return sum(1 << (5 * a + b) for a, b in pairs), None


@pytest.mark.parametrize(
    ("chrom", "sex", "copies", "first"),
    [
        pytest.param("1", Sex.FEMALE, 2, 1_000, id="autosome-female"),
        pytest.param("1", Sex.MALE, 2, 1_000, id="autosome-male"),
        pytest.param("1", None, 2, 1_000, id="autosome-unknown-sex"),
        pytest.param("X", Sex.FEMALE, 2, X_POS, id="x-female"),
        pytest.param("X", Sex.MALE, 1, X_POS, id="x-male"),
    ],
)
def test_every_short_genotype_follows_the_reference_rule(chrom, sex, copies, first):
    positions = [first + 10 * i for i in range(len(ALL_SHORT_GENOTYPES))]
    rows = [
        row(f"rs{i + 1}", chrom, pos, genotype)
        for i, (pos, genotype) in enumerate(
            zip(positions, ALL_SHORT_GENOTYPES, strict=True)
        )
    ]
    blocks, diagnostics = run_encode(rows, sex=sex)

    unrestricted = TWO_ANY if copies == 2 else ONE_ANY
    allowed, kinds = [], []
    affected: dict[str, list[int]] = {"non_base_genotype": [], "ploidy_mismatch": []}
    for pos, genotype in zip(positions, ALL_SHORT_GENOTYPES, strict=True):
        mask, code = reference_domain(genotype, copies)
        allowed.append(unrestricted if mask is None else mask)
        kinds.append(MISSING if mask is None else OBSERVED)
        if code is not None:
            affected[code].append(pos)

    # 6 one-character and 36 two-character strings. Two copies: 16 calls, 4 single
    # bases, 22 holding I or D. One copy: 4 single bases and 4 doubled bases are
    # calls, 12 two-base strings are mismatched, 22 hold I or D.
    assert len(ALL_SHORT_GENOTYPES) == 42
    counts = (
        kinds.count(OBSERVED),
        len(affected["ploidy_mismatch"]),
        len(affected["non_base_genotype"]),
    )
    assert counts == ((16, 4, 22) if copies == 2 else (8, 12, 22))

    assert list(blocks) == [chrom]
    assert_block(
        blocks[chrom],
        chrom=chrom,
        pos=positions,
        rsid=[f"rs{i + 1}" for i in range(len(positions))],
        allowed=allowed,
        kinds=kinds,
    )
    expected = sorted(
        (
            region_diag(code, chrom, min(found), max(found), len(found))
            for code, found in affected.items()
        ),
        key=lambda diag: (diag[3], diag[0]),
    )
    assert diagnostic_rows(diagnostics) == expected


RANDOM_CHROMS = ["1", "1", "1", "2", "X", "X", "Y", "MT", "PAR"]
RANDOM_FIRST = {
    "1": 1_000,
    "2": 1_000,
    "X": 59_990,
    "Y": 3_000_000,
    "MT": 100,
    "PAR": 1,
}
RANDOM_GENOTYPES = [
    "AG",
    "GA",
    "AA",
    "CC",
    "CT",
    "TC",
    "GT",
    "TT",
    "A",
    "C",
    "II",
    "DI",
    "D",
    "ag",
    "AAA",
    None,
    None,
]


# Always present, so every property run sees a phase conflict, agreeing phased
# aliases and a flagged alias.
FIXED_ROWS = [
    row("rs9001", "2", 5_000, "AG", phased=True, phase_set=1),
    row("rs9002", "2", 5_000, "GA", phased=True, phase_set=1),
    row("rs9003", "2", 5_100, "CT", phased=True),
    row("rs9004", "2", 5_100, "CT", phased=True),
    row("rs9005", "2", 5_200, "TC", phased=True, phase_set=1),
    row("rs9006", "2", 5_300, None, conflicted=True),
    row("rs9007", "2", 5_300, "GG"),
]


def random_rows(rng: np.random.Generator, n: int) -> list[tuple]:
    """Synthetic rows with repeats, conflicts, flags, phase and excluded regions."""
    names: dict[tuple[str, int, int], str] = {}
    rows = list(FIXED_ROWS)
    for _ in range(n):
        chrom = RANDOM_CHROMS[int(rng.integers(len(RANDOM_CHROMS)))]
        # X positions straddle the start of PAR1 at 60,001.
        pos = RANDOM_FIRST[chrom] + int(rng.integers(0, 60))
        if chrom == "1" and rng.random() < 0.15:
            pos = 249_250_610 + int(
                rng.integers(0, 20)
            )  # straddles the end at 249,250,621
        conflicted = bool(rng.random() < 0.05)
        genotype = (
            None
            if conflicted
            else RANDOM_GENOTYPES[int(rng.integers(len(RANDOM_GENOTYPES)))]
        )
        phased = bool(
            genotype is not None and len(genotype) == 2 and rng.random() < 0.5
        )
        phase_set = (
            int(rng.choice([0, 1, 7])) if phased and rng.random() < 0.7 else None
        )
        if rng.random() < 0.1:
            rsid = f"{chrom}:{pos}"
        else:
            key = (chrom, pos, int(rng.integers(3)))
            rsid = names.setdefault(key, f"rs{len(names) + 1}")
        rows.append(
            row(
                rsid,
                chrom,
                pos,
                genotype,
                phased=phased,
                phase_set=phase_set,
                conflicted=conflicted,
            )
        )
    return rows


def assert_rows_are_varied(rows: list[tuple]) -> None:
    positions = [(r[1], r[2]) for r in rows]
    assert len(set(positions)) < len(positions)  # repeated positions
    assert any(r[6] for r in rows)  # flagged rows
    assert any(r[4] and r[5] is None for r in rows)  # implicit set
    assert any(r[4] and r[5] is not None for r in rows)  # explicit sets
    assert {r[1] for r in rows} == set(RANDOM_FIRST)


@pytest.mark.parametrize("seed", [11, 12, 13])
@pytest.mark.parametrize("sex", [Sex.MALE, Sex.FEMALE, None])
def test_shuffling_rows_with_their_series_changes_nothing(seed, sex):
    rng = np.random.default_rng(seed)
    rows = random_rows(rng, 300)
    assert_rows_are_varied(rows)
    frame, phased, sets, conflicted = build_inputs(rows)
    blocks, diagnostics = encode(
        frame, phased=phased, phase_set=sets, genotype_conflicted=conflicted, sex=sex
    )
    order = rng.permutation(len(frame))
    shuffled_blocks, shuffled_diagnostics = encode(
        frame.iloc[order],
        phased=phased.iloc[order],
        phase_set=sets.iloc[order],
        genotype_conflicted=conflicted.iloc[order],
        sex=sex,
    )
    assert summaries(shuffled_blocks) == summaries(blocks)
    assert diagnostic_rows(shuffled_diagnostics) == diagnostic_rows(diagnostics)

    # FIXED_ROWS on chr2. Set numbers there depend on the random rows before 5000,
    # so only their separation is checked: the implicit set at 5100, set 1 at 5200.
    chr2 = blocks["2"]
    at = {pos: i for i, pos in enumerate(chr2.loci.pos.tolist())}
    assert chr2.allowed[at[5_000]] == AG  # opposite orders: unordered call kept
    assert chr2.kinds[at[5_000]].tolist() == list(OBSERVED)
    assert chr2.allowed[at[5_100]] == CT_ORDERED  # bit 8 (C,T)
    assert chr2.kinds[at[5_100]].tolist() == list(ASSERTED)
    assert chr2.allowed[at[5_200]] == TC_ORDERED  # bit 16 (T,C)
    assert chr2.kinds[at[5_200]].tolist() == list(ASSERTED)
    assert chr2.phase_set[at[5_100]] != chr2.phase_set[at[5_200]]
    assert chr2.allowed[at[5_300]] == TWO_ANY  # flagged alias suppresses the call
    assert chr2.kinds[at[5_300]].tolist() == list(MISSING)
    found = diagnostic_rows(diagnostics)
    assert locus_diag("phase_conflict", "2", 5_000) in found
    assert locus_diag("position_conflict", "2", 5_300) in found


def relabel_mask(mask: int, sigma: list[int]) -> int:
    """Move pair (a, b) to (sigma[a], sigma[b]); ABSENT stays in slot 4."""
    out = 0
    for a in range(5):
        for b in range(5):
            if mask >> (5 * a + b) & 1:
                a2 = sigma[a] if a < 4 else NO_COPY
                b2 = sigma[b] if b < 4 else NO_COPY
                out |= 1 << (5 * a2 + b2)
    return out


def test_relabel_mask_moves_pairs():
    sigma = [1, 2, 3, 0]  # A->C, C->G, G->T, T->A
    assert relabel_mask(AG, sigma) == CT  # (A,G),(G,A) -> (C,T),(T,C)
    assert relabel_mask(A_ONE, sigma) == C_ONE  # (A,ABSENT) -> (C,ABSENT)
    assert relabel_mask(TWO_ANY, sigma) == TWO_ANY
    assert relabel_mask(ONE_ANY, sigma) == ONE_ANY


@pytest.mark.parametrize(
    ("seed", "sex"), [(21, Sex.MALE), (22, Sex.FEMALE), (23, None)]
)
def test_relabeling_bases_permutes_the_masks(seed, sex):
    rng = np.random.default_rng(seed)
    sigma = [int(v) for v in rng.permutation(4)]
    if sigma == [0, 1, 2, 3]:
        sigma = [1, 2, 3, 0]
    letters = {"ACGT"[i]: "ACGT"[sigma[i]] for i in range(4)}
    rows = random_rows(rng, 300)
    relabeled = [
        (
            r[0],
            r[1],
            r[2],
            None if r[3] is None else "".join(letters.get(ch, ch) for ch in r[3]),
            *r[4:],
        )
        for r in rows
    ]
    blocks, diagnostics = run_encode(rows, sex=sex)
    relabeled_blocks, relabeled_diagnostics = run_encode(relabeled, sex=sex)
    expected = [
        (chrom, {**s, "allowed": [relabel_mask(m, sigma) for m in s["allowed"]]})
        for chrom, s in summaries(blocks)
    ]
    assert summaries(relabeled_blocks) == expected
    assert diagnostic_rows(relabeled_diagnostics) == diagnostic_rows(diagnostics)


NATURAL_CHROMS = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]


def recast(frame: pd.DataFrame, kind: str) -> pd.DataFrame:
    """Recast ``chrom`` and ``genotype`` of a frame built by make_frame (object columns)."""
    out = frame.copy()
    chroms = out["chrom"].tolist()
    genotypes = out["genotype"].tolist()
    if kind == "object":
        pass
    elif kind == "str":
        out["chrom"] = out["chrom"].astype("str")
        out["genotype"] = out["genotype"].astype("str")
    else:
        others = sorted(set(chroms) - set(NATURAL_CHROMS))
        out["chrom"] = pd.Categorical(
            chroms, categories=NATURAL_CHROMS + others, ordered=True
        )
        present = {g for g in genotypes if isinstance(g, str)}
        out["genotype"] = pd.Categorical(
            genotypes, categories=sorted(present | {"CG", "TT"})
        )
    return out


@pytest.mark.parametrize("sex", [Sex.MALE, Sex.FEMALE, None])
def test_object_str_and_categorical_columns_agree(sex):
    rng = np.random.default_rng(31)
    rows = random_rows(rng, 300)
    frame, phased, sets, conflicted = build_inputs(rows)
    results = {}
    for kind in ("object", "str", "categorical"):
        recast_frame = recast(frame, kind)
        if kind == "object":
            assert recast_frame["chrom"].dtype == object
            assert recast_frame["genotype"].dtype == object
        elif kind == "str":
            assert recast_frame["chrom"].dtype == "str"
            assert recast_frame["genotype"].dtype == "str"
        else:
            assert isinstance(recast_frame["chrom"].dtype, pd.CategoricalDtype)
            assert isinstance(recast_frame["genotype"].dtype, pd.CategoricalDtype)
        assert recast_frame["genotype"].isna().tolist() == [r[3] is None for r in rows]
        blocks, diagnostics = encode(
            recast_frame,
            phased=phased,
            phase_set=sets,
            genotype_conflicted=conflicted,
            sex=sex,
        )
        results[kind] = (summaries(blocks), diagnostic_rows(diagnostics))
    assert results["str"] == results["object"]
    assert results["categorical"] == results["object"]


def test_categorical_columns_match_expected_masks():
    rows = [
        row("rs1", "X", 60_000, "A"),
        row("rs2", "X", 60_001, "AG"),  # PAR1
        row("rs3", "1", 100, "GA"),
        row("rs4", "1", 100, "AG"),
        row("rs5", "1", 200, None),
    ]
    frame, phased, sets, conflicted = build_inputs(rows)
    blocks, diagnostics = encode(
        recast(frame, "categorical"),
        phased=phased,
        phase_set=sets,
        genotype_conflicted=conflicted,
        sex=Sex.MALE,
    )
    assert list(blocks) == ["1", "X"]
    assert_block(
        blocks["1"],
        chrom="1",
        pos=[100, 200],
        rsid=["rs3;rs4", "rs5"],
        allowed=[AG, TWO_ANY],
        kinds=[OBSERVED, MISSING],
    )
    assert_block(
        blocks["X"],
        chrom="X",
        pos=[60_000],
        rsid=["rs1"],
        allowed=[A_ONE],
        kinds=[OBSERVED],
    )
    assert diagnostic_rows(diagnostics) == [
        region_diag("excluded_region", "X", 60_001, 60_001, 1)
    ]


# --------------------------------------------------------------------------------------
# 10. Ownership
# --------------------------------------------------------------------------------------


def test_encode_leaves_its_inputs_alone_and_ignores_later_changes():
    rows = [
        row("rs1", "1", 100, "AG", phased=True, phase_set=4),
        row("rs2", "1", 200, "CT", phased=True, phase_set=4),
        row("rs3", "1", 300, None, conflicted=True),
        row("rs4", "1", 300, "GG"),
        row("rs5", "X", X_POS, "A"),
    ]
    frame, phased, sets, conflicted = build_inputs(rows)
    copies = [
        frame.copy(deep=True),
        phased.copy(deep=True),
        sets.copy(deep=True),
        conflicted.copy(deep=True),
    ]
    blocks, diagnostics = encode(
        frame,
        phased=phased,
        phase_set=sets,
        genotype_conflicted=conflicted,
        sex=Sex.MALE,
    )
    pd.testing.assert_frame_equal(frame, copies[0])
    pd.testing.assert_series_equal(phased, copies[1])
    pd.testing.assert_series_equal(sets, copies[2])
    pd.testing.assert_series_equal(conflicted, copies[3])

    expected = [
        (
            "1",
            {
                "chrom": "1",
                "pos": [100, 200, 300],
                "rsid": ["rs1", "rs2", "rs3;rs4"],
                "allowed": [AG_ORDERED, CT_ORDERED, TWO_ANY],
                "phase_set": [0, 0, -1],
                "kinds": [list(ASSERTED), list(ASSERTED), list(MISSING)],
            },
        ),
        (
            "X",
            {
                "chrom": "X",
                "pos": [X_POS],
                "rsid": ["rs5"],
                "allowed": [A_ONE],
                "phase_set": [-1],
                "kinds": [list(OBSERVED)],
            },
        ),
    ]
    assert summaries(blocks) == expected
    expected_diagnostics = [locus_diag("position_conflict", "1", 300)]
    assert diagnostic_rows(diagnostics) == expected_diagnostics
    diagnostics_copy = diagnostics.copy(deep=True)

    frame.iloc[0, frame.columns.get_loc("genotype")] = "TT"
    frame.iloc[1, frame.columns.get_loc("pos")] = 250
    frame.iloc[4, frame.columns.get_loc("chrom")] = "2"
    frame.index = pd.Index(["rs9"] * len(frame), dtype="str", name="rsid")
    phased.iloc[:] = False
    sets.iloc[:] = 9
    conflicted.iloc[:] = False

    assert summaries(blocks) == expected
    pd.testing.assert_frame_equal(diagnostics, diagnostics_copy)


def test_calls_keeps_owned_copies_of_its_arrays():
    loci = Loci("1", [100, 200])
    allowed = np.array([AG_ORDERED, CT_ORDERED], dtype=np.uint32)
    phase_set = np.array([0, 0], dtype=np.int32)
    kinds = np.array([[1, 2, 2], [1, 2, 2]], dtype=np.uint8)
    calls = _Calls(loci, allowed, phase_set, kinds)

    for stored, given in [
        (calls.allowed, allowed),
        (calls.phase_set, phase_set),
        (calls.kinds, kinds),
    ]:
        assert not np.shares_memory(stored, given)
        assert given.flags.writeable  # the caller's array is left as it was

    allowed[:] = AG
    phase_set[:] = 5
    kinds[:] = 3
    assert calls.allowed.tolist() == [AG_ORDERED, CT_ORDERED]
    assert calls.phase_set.tolist() == [0, 0]
    assert calls.kinds.tolist() == [[1, 2, 2], [1, 2, 2]]


def _encoded_block() -> _Calls:
    blocks, _ = run_encode(
        [row("rs1", "1", 100, "AG", phased=True), row("rs2", "1", 200, None)],
        sex=Sex.FEMALE,
    )
    return blocks["1"]


def _aligned_block() -> _Calls:
    out = align(
        {
            "p": {"X": _male_x_block()},
            "q": {"X": block("X", [4_000_000], ["rs3"], [AG], [OBSERVED])},
        },
        sexes={"p": Sex.MALE, "q": Sex.FEMALE},
    )
    return out["p"]["X"]


def _constructed_block() -> _Calls:
    return make_calls(**autosome_args())


@pytest.mark.parametrize(
    "build",
    [
        pytest.param(_encoded_block, id="encode"),
        pytest.param(_aligned_block, id="align"),
        pytest.param(_constructed_block, id="constructor"),
    ],
)
def test_stored_arrays_are_read_only(build):
    calls = build()
    for array in (
        calls.allowed,
        calls.phase_set,
        calls.kinds,
        calls.loci.pos,
        calls.loci.rsid,
    ):
        assert_read_only(array)
    with pytest.raises(AttributeError):
        calls.allowed = calls.allowed  # type: ignore[misc]


def test_returned_mappings_cannot_be_assigned_to():
    blocks, _ = run_encode([row("rs1", "1", 100, "AG")], sex=Sex.FEMALE)
    assert_read_only_mapping(blocks, "1")
    empty, _ = encode(make_frame([]), sex=None)
    with pytest.raises(TypeError):
        empty["1"] = blocks["1"]  # type: ignore[index]

    aligned = align({"p": blocks}, sexes={"p": Sex.FEMALE, "q": Sex.MALE})
    assert_read_only_mapping(aligned, "p")
    assert_read_only_mapping(aligned["p"], "1")
    assert_read_only_mapping(aligned["q"], "1")


def test_align_ignores_later_changes_to_its_inputs():
    p_blocks = {"1": block("1", [100], ["rs1"], [AG], [OBSERVED])}
    people = {"p": p_blocks}
    sexes = {"p": Sex.FEMALE, "q": Sex.MALE}
    out = align(people, sexes=sexes)
    before = {name: summaries(blocks) for name, blocks in out.items()}

    p_blocks["2"] = block("2", [50], ["rs2"], [CC], [OBSERVED])
    p_blocks["1"] = block("1", [300], ["rs3"], [GG], [OBSERVED])
    people["r"] = {"1": block("1", [400], ["rs4"], [TT], [OBSERVED])}
    sexes["r"] = Sex.FEMALE
    sexes["q"] = None

    assert list(out) == ["p", "q"]
    assert {name: summaries(blocks) for name, blocks in out.items()} == before
    assert before["q"] == [
        (
            "1",
            {
                "chrom": "1",
                "pos": [100],
                "rsid": ["rs1"],
                "allowed": [TWO_ANY],
                "phase_set": [-1],
                "kinds": [list(MISSING)],
            },
        )
    ]


# --------------------------------------------------------------------------------------
# 11. Further input rules, type checks and copies
# --------------------------------------------------------------------------------------

# Invalid positions are checked first and count nowhere else; rs9 names none of them.
INVALID_FIRST_ROWS = [
    row("rs1", "1", 100, "AG"),
    row("rs9", "1", 249_250_622, "II"),  # non-base, but beyond chr1
    row("rs9", "1", 249_250_623, "A"),  # one character, but beyond chr1
    row("rs4", "X", 3_000_000, "AG"),
    row("rs9", "X", 155_270_561, "AG"),  # beyond X
]


@pytest.mark.parametrize(
    ("sex", "chroms", "x_diagnostics"),
    [
        pytest.param(
            None,
            ["1"],
            [
                region_diag("excluded_region", "X", 3_000_000, 3_000_000, 1),
                locus_diag("invalid_position", "X", 155_270_561),
            ],
            id="unknown-sex",
        ),
        pytest.param(
            Sex.FEMALE,
            ["1", "X"],
            [locus_diag("invalid_position", "X", 155_270_561)],
            id="female",
        ),
    ],
)
def test_invalid_positions_are_checked_first_and_count_nowhere_else(
    sex, chroms, x_diagnostics
):
    blocks, diagnostics = run_encode(INVALID_FIRST_ROWS, sex=sex)
    assert list(blocks) == chroms
    assert blocks["1"].loci.pos.tolist() == [100]
    assert blocks["1"].loci.rsid.tolist() == ["rs1"]
    if "X" in blocks:
        assert blocks["X"].loci.pos.tolist() == [3_000_000]
        assert blocks["X"].loci.rsid.tolist() == ["rs4"]
    assert diagnostic_rows(diagnostics) == [
        locus_diag("invalid_position", "1", 249_250_622),
        locus_diag("invalid_position", "1", 249_250_623),
        *x_diagnostics,
    ]


def _three_rows() -> pd.DataFrame:
    return make_frame(
        [
            row("rs1", "1", 100, "AG"),
            row("rs1", "1", 200, "CT"),
            row("rs2", "1", 300, "GG"),
        ]
    )


def _series(values, dtype, index) -> pd.Series:
    return pd.Series(pd.array(values, dtype=dtype), index=pd.Index(index, dtype="str"))


MISALIGNED_OR_MISTYPED = [
    pytest.param("phased", [False, False], bool, ["rs1", "rs1"], id="too-short"),
    pytest.param(
        "phased", [False] * 3, bool, ["rs1", "rs2", "rs3"], id="other-index-values"
    ),
    # The same identifiers in another order: indexes are compared row by row.
    pytest.param(
        "genotype_conflicted",
        [False] * 3,
        bool,
        ["rs1", "rs2", "rs1"],
        id="repeated-rsids-reordered",
    ),
    pytest.param(
        "phase_set", [1, 1, 1], "Int64", ["rs2", "rs1", "rs1"], id="phase-set-index"
    ),
    pytest.param(
        "phased", [True, None, False], "boolean", ["rs1", "rs1", "rs2"], id="boolean-na"
    ),
    pytest.param(
        "genotype_conflicted",
        [False, True, None],
        "boolean",
        ["rs1", "rs1", "rs2"],
        id="conflicted-boolean-na",
    ),
    pytest.param(
        "phased", [0.0, 1.0, 0.0], "float64", ["rs1", "rs1", "rs2"], id="float-phased"
    ),
    pytest.param(
        "phased", [True, False, True], object, ["rs1", "rs1", "rs2"], id="object-bools"
    ),
    pytest.param(
        "phase_set", [1.0, 2.0, 3.0], "float64", ["rs1", "rs1", "rs2"], id="float-set"
    ),
    pytest.param(
        "phase_set", [True, True, True], bool, ["rs1", "rs1", "rs2"], id="bool-set"
    ),
    # Negative on an unphased row still raises.
    pytest.param(
        "phase_set", [1, None, -1], "Int64", ["rs1", "rs1", "rs2"], id="negative-set"
    ),
    pytest.param(
        "phase_set", [1, 2, -3], "int64", ["rs1", "rs1", "rs2"], id="negative-np-set"
    ),
]


@pytest.mark.parametrize(("name", "values", "dtype", "index"), MISALIGNED_OR_MISTYPED)
def test_encode_rejects_misaligned_or_mistyped_series(name, values, dtype, index):
    with pytest.raises(ValueError):
        encode(_three_rows(), sex=Sex.FEMALE, **{name: _series(values, dtype, index)})


def test_series_are_aligned_row_by_row_not_by_index():
    # Both rows are named rs1; only the second is flagged.
    frame = make_frame([row("rs1", "1", 100, "AG"), row("rs1", "1", 200, None)])
    flags = bool_series(frame, [False, True])
    blocks, diagnostics = encode(frame, genotype_conflicted=flags, sex=Sex.FEMALE)
    assert_block(
        blocks["1"],
        chrom="1",
        pos=[100, 200],
        rsid=["rs1", "rs1"],
        allowed=[AG, TWO_ANY],
        kinds=[OBSERVED, MISSING],
    )
    assert diagnostic_rows(diagnostics) == [locus_diag("position_conflict", "1", 200)]


@pytest.mark.parametrize("set_dtype", ["Int64", "Int32", "UInt8", "int64", "uint16"])
@pytest.mark.parametrize("flag_dtype", [bool, "boolean"])
def test_encode_accepts_integer_sets_and_boolean_flags(set_dtype, flag_dtype):
    frame = _three_rows()
    index = frame.index.tolist()
    blocks, diagnostics = encode(
        frame,
        phased=_series([True, True, False], flag_dtype, index),
        phase_set=_series([7, 7, 2], set_dtype, index),
        genotype_conflicted=_series([False, False, False], flag_dtype, index),
        sex=Sex.FEMALE,
    )
    assert_block(
        blocks["1"],
        chrom="1",
        pos=[100, 200, 300],
        rsid=["rs1", "rs1", "rs2"],
        allowed=[AG_ORDERED, CT_ORDERED, GG],
        kinds=[ASSERTED, ASSERTED, OBSERVED],
        phase_set=[0, 0, -1],
    )
    assert diagnostic_rows(diagnostics) == []


@pytest.mark.parametrize(
    "other",
    [
        pytest.param(row("rs2", "1", 200, ["A", "G"]), id="list"),
        pytest.param(row("rs2", "1", 200, 5), id="int"),
        pytest.param(row("rs2", "MT", 5, 5), id="int-on-mt"),
    ],
)
def test_non_str_genotypes_raise(other):
    with pytest.raises(ValueError):
        run_encode([row("rs1", "1", 100, "AG"), other], sex=Sex.FEMALE)


def test_every_diagnostic_has_a_detail():
    rows = [
        row("rs1", "1", 100, "AG"),
        row("rs2", "1", 100, "CC"),  # position_conflict
        row("rs3", "1", 200, "AG", phased=True, phase_set=1),
        row("rs4", "1", 200, "GA", phased=True, phase_set=1),  # phase_conflict
        row("rs5", "1", 300, "II"),  # non_base_genotype
        row("rs6", "1", 400, "A"),  # ploidy_mismatch
        row("rs7", "1", 249_250_622, "AG"),  # invalid_position
        row("rs8", "Y", 10, "A"),  # excluded_region
    ]
    _, diagnostics = run_encode(rows, sex=Sex.FEMALE)
    assert diagnostics["code"].tolist() == [
        "position_conflict",
        "phase_conflict",
        "non_base_genotype",
        "ploidy_mismatch",
        "invalid_position",
        "excluded_region",
    ]
    for detail in diagnostics["detail"]:
        assert isinstance(detail, str)
        assert detail


@pytest.mark.parametrize(
    "chroms",
    [
        pytest.param([1, 1.0], id="int-then-float"),
        pytest.param([1.0, 1], id="float-then-int"),
        pytest.param([["1"], "1"], id="list"),
    ],
)
def test_non_str_chromosomes_raise(chroms):
    rows = [row(f"rs{i}", chrom, 100 * (i + 1), "AG") for i, chrom in enumerate(chroms)]
    with pytest.raises(ValueError):
        run_encode(rows, sex=Sex.FEMALE)


@pytest.mark.parametrize(
    ("values", "categories"),
    [
        pytest.param(["1", "Y"], ["1", "Y", 1], id="unused-int-category"),
        pytest.param([1, 24], [1, 24], id="int-categories"),
    ],
)
def test_non_str_chromosome_categories_raise(values, categories):
    frame = make_frame([row("rs1", "1", 100, "AG"), row("rs2", "Y", 10, "A")])
    frame["chrom"] = pd.Categorical(values, categories=categories)
    with pytest.raises(ValueError):
        encode(frame, sex=Sex.MALE)


@pytest.mark.parametrize(
    "values",
    [
        pytest.param([100.0, 200.0], id="float64"),
        pytest.param(pd.array([100, 200], dtype="Int64"), id="Int64"),
    ],
)
def test_encode_accepts_whole_positions_of_other_dtypes(values):
    frame = make_frame([row("rs1", "1", 1, "AG"), row("rs2", "1", 2, "CC")])
    frame["pos"] = values
    blocks, _ = encode(frame, sex=Sex.FEMALE)
    assert blocks["1"].loci.pos.tolist() == [100, 200]


@pytest.mark.skipif(
    np.finfo(np.longdouble).eps >= np.finfo(np.float64).eps,
    reason="longdouble is no more precise than float64 here",
)
def test_a_longdouble_position_just_above_a_whole_number_is_rejected():
    # 100 + 2**-50 is 100 in float64, but not in a wider longdouble.
    frame = _one_row_frame()
    frame["pos"] = np.array([np.longdouble(100) + np.longdouble(2) ** -50])
    assert frame["pos"].dtype == np.longdouble
    assert frame["pos"].iloc[0] != 100
    with pytest.raises(ValueError):
        encode(frame, sex=Sex.FEMALE)


def _two_rows() -> pd.DataFrame:
    return make_frame([row("rs1", "1", 100, "AG"), row("rs2", "1", 200, "CC")])


def _with_index(values) -> pd.DataFrame:
    frame = _two_rows()
    frame.index = pd.Index(values, dtype="str")
    return frame


@pytest.mark.parametrize(
    "build",
    [
        pytest.param(lambda: _with_index(["rs1", None]), id="missing-identifier"),
        pytest.param(lambda: _with_index(["rs1", ""]), id="empty-identifier"),
        pytest.param(
            lambda: pd.concat([_one_row_frame(), _one_row_frame()["chrom"]], axis=1),
            id="repeated-column",
        ),
        pytest.param(lambda: _one_row_frame().assign(pos=["100"]), id="text-positions"),
        pytest.param(lambda: _one_row_frame().assign(pos=[True]), id="bool-positions"),
        pytest.param(
            lambda: _two_rows().assign(genotype=pd.Categorical([5, None])),
            id="genotype-int-category",
        ),
        pytest.param(
            lambda: _two_rows().set_axis(pd.CategoricalIndex([1, 2])),
            id="index-int-categories",
        ),
    ],
)
def test_encode_rejects_bad_identifiers_columns_and_positions(build):
    with pytest.raises(ValueError):
        encode(build(), sex=Sex.FEMALE)


@pytest.mark.parametrize(
    "index",
    [
        pytest.param(
            pd.MultiIndex.from_tuples([("rs1", "a"), ("rs2", "b")]), id="two-levels"
        ),
        pytest.param(pd.MultiIndex.from_arrays([["rs1", "rs2"]]), id="one-level"),
    ],
)
def test_a_multi_index_raises_type_error(index):
    with pytest.raises(TypeError):
        encode(_two_rows().set_axis(index), sex=Sex.FEMALE)


@pytest.mark.parametrize(
    "build",
    [
        pytest.param(
            lambda: pd.DataFrame({"chrom": [], "pos": [], "genotype": []}),
            id="float64-columns",
        ),
        pytest.param(
            lambda: pd.DataFrame(columns=["chrom", "pos", "genotype"]),
            id="object-columns",
        ),
    ],
)
def test_an_untyped_empty_frame_is_valid(build):
    frame = build()
    frame.index = pd.Index([], dtype="str")
    blocks, diagnostics = encode(frame, sex=Sex.FEMALE)
    assert len(blocks) == 0
    assert len(diagnostics) == 0
    assert_diagnostics_schema(diagnostics)


@pytest.mark.parametrize(
    "genotypes",
    [
        pytest.param(np.full(2, np.nan), id="float-nan"),
        pytest.param(pd.Categorical([None, None]), id="categorical-no-categories"),
    ],
)
def test_an_all_null_genotype_column_is_valid(genotypes):
    frame = _two_rows()
    frame["genotype"] = genotypes
    blocks, diagnostics = encode(frame, sex=Sex.FEMALE)
    assert_block(
        blocks["1"],
        chrom="1",
        pos=[100, 200],
        rsid=["rs1", "rs2"],
        allowed=[TWO_ANY, TWO_ANY],
        kinds=[MISSING, MISSING],
    )
    assert diagnostic_rows(diagnostics) == []


def test_calls_rejects_a_bool_in_plain_input():
    args = autosome_args()
    with pytest.raises(TypeError):
        _Calls(**{**args, "allowed": [AG, True]})
    _Calls(**args)  # the same lists without the bool are accepted


def test_calls_copies_stay_read_only():
    calls = make_calls(
        Loci("1", [100, 200], ["rs1", None]),
        [AG_ORDERED, CT_ORDERED],
        [0, 0],
        [[1, 2, 2], [1, 2, 2]],
    )
    for copied in (
        copy.copy(calls),
        copy.deepcopy(calls),
        pickle.loads(pickle.dumps(calls)),
    ):
        assert isinstance(copied, _Calls)
        assert summary(copied) == summary(calls)
        for name in ("allowed", "phase_set", "kinds"):
            assert_read_only(getattr(copied, name))


@pytest.mark.parametrize(
    ("identifiers", "expected"),
    [
        pytest.param(["rs2;rs1", "rs1", "rs10"], "rs1;rs10;rs2", id="split-merged"),
        pytest.param(
            ["rs2;rs1", "rs1", "1:100;rs10"], "rs1;rs10;rs2", id="placeholder-piece"
        ),
        pytest.param(["1:100;rs5"], "rs5", id="own-placeholder-dropped"),
        pytest.param(["1:101;rs5"], "1:101;rs5", id="other-placeholder-kept"),
        pytest.param(["rs1;;rs2;"], "rs1;rs2", id="empty-pieces-dropped"),
        pytest.param([";"], None, id="only-separators"),
        pytest.param(["1:100;", "1:100"], None, id="only-placeholders"),
    ],
)
def test_identifiers_are_lists_of_names(identifiers, expected):
    rows = [row(name, "1", 100, "AG") for name in identifiers]
    rows.append(row("rs99", "1", 200, "CC"))
    blocks, diagnostics = run_encode(rows, sex=Sex.FEMALE)
    assert blocks["1"].loci.pos.tolist() == [100, 200]
    assert blocks["1"].loci.rsid.tolist() == [expected, "rs99"]
    assert diagnostic_rows(diagnostics) == []


def test_align_of_one_person_leaves_every_name_unchanged():
    rows = [
        row("rs2;rs1", "1", 100, "AG"),
        row("1:100;rs10", "1", 100, "AG"),
        row("1:200", "1", 200, "CC"),
        row("rs7;;rs3", "1", 300, "GT"),
        row("rs5;X:3000000", "X", X_POS, "CT"),
    ]
    blocks, _ = run_encode(rows, sex=Sex.FEMALE)
    assert [blocks[c].loci.rsid.tolist() for c in blocks] == [
        ["rs1;rs10;rs2", None, "rs3;rs7"],
        ["rs5"],
    ]
    aligned = align({"p": blocks}, sexes={"p": Sex.FEMALE})
    assert summaries(aligned["p"]) == summaries(blocks)


def test_align_splits_names_like_encode():
    # "rs9;1:100" holds the placeholder of its own position; ";" holds no name.
    p = {
        "1": block("1", [100, 200], ["rs9;1:100", "rs3;;rs3"], [AG, CC], [OBSERVED] * 2)
    }
    q = {"1": block("1", [100, 300], [";", "rs1"], [GG, TT], [OBSERVED] * 2)}
    out = align({"p": p, "q": q}, sexes={"p": Sex.FEMALE, "q": Sex.MALE})
    assert out["p"]["1"].loci.pos.tolist() == [100, 200, 300]
    assert out["p"]["1"].loci.rsid.tolist() == ["rs9", "rs3", "rs1"]


def test_rows_at_a_suppressed_position_still_count_in_row_diagnostics():
    rows = [
        row("rs1", "1", 100, None, conflicted=True),
        row("rs2", "1", 100, "II"),
        row("rs3", "1", 100, "A"),
    ]
    blocks, diagnostics = run_encode(rows, sex=Sex.FEMALE)
    assert_block(
        blocks["1"],
        chrom="1",
        pos=[100],
        rsid=["rs1;rs2;rs3"],
        allowed=[TWO_ANY],
        kinds=[MISSING],
    )
    # One position, so the order is by code.
    assert diagnostic_rows(diagnostics) == [
        region_diag("non_base_genotype", "1", 100, 100, 1),
        region_diag("ploidy_mismatch", "1", 100, 100, 1),
        locus_diag("position_conflict", "1", 100),
    ]


def test_a_non_str_identifier_on_an_excluded_chromosome_raises():
    frame = make_frame([row("rs1", "1", 100, "AG"), row("rs2", "Y", 10, "A")])
    frame.index = pd.Index(["rs1", 5], dtype=object)
    with pytest.raises(ValueError):
        encode(frame, sex=Sex.MALE)


def test_two_bases_beyond_a_male_x_are_only_an_invalid_position():
    rows = [row("rs1", "X", X_POS, "A"), row("rs2", "X", 155_270_561, "AG")]
    blocks, diagnostics = run_encode(rows, sex=Sex.MALE)
    assert_block(
        blocks["X"],
        chrom="X",
        pos=[X_POS],
        rsid=["rs1"],
        allowed=[A_ONE],
        kinds=[OBSERVED],
    )
    assert diagnostic_rows(diagnostics) == [
        locus_diag("invalid_position", "X", 155_270_561)
    ]


@pytest.mark.parametrize(
    ("sex", "allowed"), [("Male", A_ONE), ("Female", TWO_ANY)], ids=["male", "female"]
)
def test_a_sex_given_as_text_is_that_sex(sex, allowed):
    blocks, _ = run_encode([row("rs1", "X", X_POS, "A")], sex=sex)
    assert blocks["X"].allowed.tolist() == [allowed]
    out = align({"p": blocks}, sexes={"p": sex})
    assert out["p"]["X"].allowed.tolist() == [allowed]


@pytest.mark.parametrize(
    "sex",
    ["male", "M", "", 1, True],
    ids=["lower-case", "letter", "empty", "int", "bool"],
)
def test_an_invalid_sex_raises(sex):
    with pytest.raises(ValueError):
        run_encode([row("rs1", "1", 100, "AG")], sex=sex)
    with pytest.raises(ValueError):
        align({"p": {"1": _autosome_block()}}, sexes={"p": sex})
    with pytest.raises(ValueError):
        align({}, sexes={"q": sex})


@pytest.mark.parametrize("origin", [2, 3], ids=["forced", "estimated"])
def test_origin_may_be_known_on_a_one_copy_row(origin):
    calls = make_calls(
        Loci("X", [X_POS, 4_000_000]),
        [A_ONE, ONE_ANY],
        [-1, -1],
        [[1, 0, origin], [0, 0, origin]],
    )
    assert calls.kinds.tolist() == [[1, 0, origin], [0, 0, origin]]


def test_origin_may_be_known_throughout_a_phase_set():
    calls = make_calls(
        Loci("1", [100, 200, 300]),
        [AG_ORDERED, CT_ORDERED, AG],
        [0, 0, -1],
        [[1, 2, 3], [1, 2, 3], [1, 0, 0]],
    )
    assert calls.kinds[:, 2].tolist() == [3, 3, 0]
