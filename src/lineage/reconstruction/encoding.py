"""Genotype calls encoded for reconstruction: :class:`_Calls` blocks and their builders."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from types import MappingProxyType
from typing import Any, Final, NamedTuple, Self, TypeAlias

import numpy as np
import numpy.typing as npt
import pandas as pd

from lineage import _grch37
from lineage.types import ABSENT, UNKNOWN, Haplotypes, Loci, Sex

_Bools: TypeAlias = npt.NDArray[np.bool_]
_Indices: TypeAlias = npt.NDArray[np.intp]
_Ints: TypeAlias = npt.NDArray[np.integer[Any]]
_Masks: TypeAlias = npt.NDArray[np.uint32]
_Objects: TypeAlias = npt.NDArray[np.object_]
_Positions: TypeAlias = npt.NDArray[np.uint32]
_Fields: TypeAlias = tuple[
    Loci,
    npt.NDArray[np.uint32],
    npt.NDArray[np.int32],
    npt.NDArray[np.uint8],
]

_BASES: Final = "ACGT"
_NO_COPY: Final = 4
_SLOTS: Final = 5


def _pair(a: int, b: int) -> int:
    """Return the bit of a mask that stands for the pair of slots (`a`, `b`)."""
    return 1 << (_SLOTS * a + b)


_TWO_ANY: Final = sum(_pair(a, b) for a in range(_NO_COPY) for b in range(_NO_COPY))
_ONE_ANY: Final = sum(_pair(a, _NO_COPY) for a in range(_NO_COPY))

#: The mask that allows every pair, by number of copies.
_UNRESTRICTED: Final[Mapping[int, int]] = MappingProxyType({1: _ONE_ANY, 2: _TWO_ANY})

_ALL_BITS: Final = (1 << _SLOTS * _SLOTS) - 1

#: The pairs without a copy 0, which no mask allows.
_NO_COPY_0: Final = sum(_pair(_NO_COPY, b) for b in range(_SLOTS))

#: For each copy and slot, every pair that puts that slot on that copy.
_HOLDS: Final = tuple(
    tuple(
        sum(
            _pair(a, b)
            for a in range(_SLOTS)
            for b in range(_SLOTS)
            if (a, b)[copy] == slot
        )
        for slot in range(_SLOTS)
    )
    for copy in range(2)
)

_MISSING: Final = 0
_OBSERVED: Final = 1
_FORCED: Final = 2
_ESTIMATED: Final = 3
_ASSERTED: Final = 4

#: The columns of a block's kinds, in the order of _ATTRIBUTES.
_GENOTYPE: Final = 0
_PHASE: Final = 1
_ORIGIN: Final = 2

#: The attribute of each column of a block's kinds, with the codes it admits.
_ATTRIBUTES: Final = (
    ("genotype", (_MISSING, _OBSERVED, _FORCED, _ESTIMATED)),
    ("phase", (_MISSING, _FORCED, _ESTIMATED, _ASSERTED)),
    ("origin", (_MISSING, _FORCED, _ESTIMATED)),
)

#: The position of each chromosome name in natural order; other names follow, by text.
_NATURAL: Final[Mapping[str, int]] = MappingProxyType(
    {name: i for i, name in enumerate((*_grch37.LENGTHS, "Y", "MT"))}
)

#: The diagnostics a row's genotype can raise, each with its detail.
_ROW_CODES: Final[Mapping[str, str]] = MappingProxyType(
    {
        "non_base_genotype": (
            "a genotype other than one or two of the bases A, C, G and T is a no-call"
        ),
        "ploidy_mismatch": "a genotype whose bases do not fit {copies} is a no-call",
    }
)

_MAX_POS: Final = int(np.iinfo(np.uint32).max)
_MAX_PHASE_SET: Final = int(np.iinfo(np.int32).max)


@dataclass(frozen=True, init=False, eq=False, slots=True)
class _Calls:
    """One person's calls on one chromosome.

    Bit ``5*a + b`` of a mask stands for copy 0 holding `a` and copy 1 holding `b`,
    where 0 to 3 are A, C, G and T and 4 is no such copy. Where origin is known, copy 0
    is from the mother and copy 1 from the father.

    Each locus has a genotype (which alleles are present), a phase (which of them share
    a chromosome copy with the other loci of its phase set) and an origin (which parent
    each copy came from). Each has a kind: 0 missing, 1 observed (read from a test,
    genotype only), 2 forced (required by inheritance), 3 estimated (inferred) or 4
    asserted (given by the input, phase only).

    Parameters
    ----------
    loci : Loci
        The locus axis, shared rather than copied.
    allowed : array_like of int, shape (n,)
        The pairs possible at each locus, as a mask.
    phase_set : array_like of int, shape (n,)
        The phase set of each locus, or -1 for none.
    kinds : array_like of int, shape (n, 3)
        The kind code of the genotype, the phase and the origin at each locus.

    Attributes
    ----------
    allowed, phase_set, kinds : numpy.ndarray
        Read-only copies, as uint32, int32 and uint8.

    Raises
    ------
    TypeError
        If `loci` is not a Loci or an array does not hold integers.
    ValueError
        If `loci` is empty or an array is malformed.
    """

    loci: Loci
    allowed: npt.NDArray[np.uint32]
    phase_set: npt.NDArray[np.int32]
    kinds: npt.NDArray[np.uint8]

    def __init__(
        self,
        loci: Loci,
        allowed: npt.ArrayLike,
        phase_set: npt.ArrayLike,
        kinds: npt.ArrayLike,
    ) -> None:
        if not isinstance(loci, Loci):
            raise TypeError(f"loci must be a Loci, not {type(loci).__name__}")
        n = len(loci.pos)
        masks = _integers(allowed, "allowed", (n,))
        sets = _integers(phase_set, "phase_set", (n,))
        codes = _integers(kinds, "kinds", (n, len(_ATTRIBUTES)))
        if n == 0:
            raise ValueError("A block must have at least one locus; loci is empty")
        checked_allowed, copies = _check_masks(masks, loci.chrom)
        checked_sets = _check_phase_sets(sets, copies)
        checked_kinds = _check_kinds(
            codes, (checked_allowed == _UNRESTRICTED[copies], checked_sets == -1)
        )
        _check_origins(checked_sets, checked_kinds)
        for array in (checked_allowed, checked_sets, checked_kinds):
            array.setflags(write=False)
        object.__setattr__(self, "loci", loci)
        object.__setattr__(self, "allowed", checked_allowed)
        object.__setattr__(self, "phase_set", checked_sets)
        object.__setattr__(self, "kinds", checked_kinds)

    def __reduce__(self) -> tuple[type[Self], _Fields]:
        return (type(self), (self.loci, self.allowed, self.phase_set, self.kinds))


def _integers(value: object, label: str, shape: tuple[int, ...]) -> _Ints:
    """Return `value` as an integer array of `shape`, without copying it."""
    if _grch37.holds_bool(value):
        raise TypeError(f"{label} must hold integers, not bool")
    array = np.asarray(value)
    if array.dtype.kind not in "iu":
        raise TypeError(f"{label} must hold integers, not {array.dtype}")
    if array.shape != shape:
        raise ValueError(
            f"{label} must have shape {shape}, one entry per locus, not {array.shape}"
        )
    return array


def _reject(invalid: _Bools, values: npt.NDArray[Any], message: str) -> None:
    """Raise a ValueError naming the first locus of a block that is invalid, if any."""
    if invalid.ndim > 1:
        invalid = invalid.any(axis=tuple(range(1, invalid.ndim)))
    if invalid.any():
        i = int(np.argmax(invalid))
        raise ValueError(f"{message}; locus {i} has {values[i]}")


def _check_masks(masks: _Ints, chrom: str) -> tuple[_Masks, int]:
    """Return a block's masks as a new uint32 array, and the copies they are for."""
    _reject(
        (masks < 1) | (masks > _ALL_BITS),
        masks,
        "Each mask in allowed must be non-zero and use only bits 0 to 24",
    )
    checked = masks.astype(np.uint32)
    _reject(
        (checked & _NO_COPY_0) != 0,
        masks,
        "No mask may allow copy 0 to be absent (bits 20 to 24)",
    )
    copies = _copies_of(checked)
    _reject(
        (copies == 1) & ((checked & _TWO_ANY) != 0),
        masks,
        "No mask may allow pairs of both one copy and two copies",
    )
    _reject(
        copies != copies[0],
        masks,
        "Every mask in a block must be for as many copies as the first",
    )
    if copies[0] == 1 and chrom != "X":
        raise ValueError(f"Only X can have one copy, not chromosome {chrom!r}")
    return checked, int(copies[0])


def _copies_of(masks: _Masks) -> _Indices:
    """Return the copies each mask is for: one if it allows a pair without copy 1."""
    return np.where((masks & _ONE_ANY) != 0, 1, 2)


def _check_phase_sets(sets: _Ints, copies: int) -> npt.NDArray[np.int32]:
    """Return a block's phase sets as a new int32 array."""
    _reject(
        (sets < -1) | (sets > _MAX_PHASE_SET),
        sets,
        f"Each phase_set must be -1 or a number from 0 to {_MAX_PHASE_SET:,}",
    )
    checked = sets.astype(np.int32)
    if copies == 1:
        _reject(checked != -1, sets, "One copy has no phase; phase_set must be -1")
    return checked


def _check_kinds(codes: _Ints, missing: tuple[_Bools, _Bools]) -> npt.NDArray[np.uint8]:
    """Return a block's kinds as uint8, given where genotype and phase are missing."""
    _reject(
        (codes < _MISSING) | (codes > _ASSERTED),
        codes,
        f"Each kind must be a code from {_MISSING} to {_ASSERTED}",
    )
    checked = codes.astype(np.uint8)
    for column, (attribute, admitted) in enumerate(_ATTRIBUTES):
        _reject(
            ~np.isin(checked[:, column], admitted),
            codes,
            f"The {attribute} kind must be one of {', '.join(map(str, admitted))}",
        )
    for column, absent in enumerate(missing):
        attribute = _ATTRIBUTES[column][0]
        _reject(
            (checked[:, column] == _MISSING) != absent,
            codes,
            f"The {attribute} kind must be 0 exactly when the {attribute} is missing",
        )
    return checked


def _check_origins(sets: npt.NDArray[np.int32], kinds: npt.NDArray[np.uint8]) -> None:
    """Check that origin is known at every locus of a phase set or at none."""
    linked = np.flatnonzero(sets >= 0)
    by_set = linked[np.argsort(sets[linked], kind="stable")]
    known = kinds[by_set, _ORIGIN] != _MISSING
    differs = (sets[by_set[1:]] == sets[by_set[:-1]]) & (known[1:] != known[:-1])
    if differs.any():
        k = int(np.argmax(differs))
        i, j = int(by_set[k]), int(by_set[k + 1])
        raise ValueError(
            f"Origin must be known at every locus of a phase set or at none; loci {i} "
            f"and {j} of phase set {sets[i]} differ"
        )


def haplotypes(calls: _Calls) -> Haplotypes:
    """Return the allele of each copy at each locus of a block.

    Parameters
    ----------
    calls : _Calls

    Returns
    -------
    Haplotypes
        A new int8 array of shape (2, n): copy 0 and copy 1 at each locus.
    """
    alleles = np.full((2, len(calls.allowed)), UNKNOWN, dtype=np.int8)
    for copy, holds in enumerate(_HOLDS):
        for slot, mask in enumerate(holds):
            only = (calls.allowed & mask) == calls.allowed
            alleles[copy, only] = slot if slot != _NO_COPY else ABSENT
    return alleles


class _Rows(NamedTuple):
    """Arrays with an entry per row of a frame, as :func:`encode` reads them."""

    pos: _Positions
    genotype: _Indices  # indexes a _Domain; -1, a null genotype, is its last entry
    phased: _Bools
    flagged: _Bools
    phase_set: _Indices  # identifies a set; -1 is the chromosome's implicit set
    name: _Indices  # indexes the identifiers given to _split

    def take(self, which: slice | _Indices | _Bools) -> _Rows:
        """Return the rows `which` selects."""
        return _Rows(
            self.pos[which],
            self.genotype[which],
            self.phased[which],
            self.flagged[which],
            self.phase_set[which],
            self.name[which],
        )


class _Domain(NamedTuple):
    """What each distinct genotype of a frame allows, for one number of copies."""

    unordered: _Masks
    ordered: _Masks
    code: npt.NDArray[np.int8]  # indexes _ROW_CODES; -1 for none


class _Pieces(NamedTuple):
    """The names in each of a list of identifiers, from :func:`_split`."""

    start: _Indices  # identifier i holds text[name[start[i]:start[i] + length[i]]]
    length: _Indices
    name: _Indices
    text: _Objects
    colon: _Bools


class _Finding(NamedTuple):
    """Diagnostics of one code and chromosome."""

    code: str
    chrom: str
    start: npt.NDArray[np.int64]
    end: npt.NDArray[np.int64]
    per_locus: bool
    detail: str


class _Report:
    """Normalization diagnostics, gathered a code and chromosome at a time."""

    def __init__(self) -> None:
        self._findings: list[_Finding] = []

    def loci(self, code: str, chrom: str, positions: _Positions, detail: str) -> None:
        """Report a finding at each of `positions`."""
        if len(positions):
            values = positions.astype(np.int64)
            self._findings.append(_Finding(code, chrom, values, values, True, detail))

    def rows(self, code: str, chrom: str, positions: _Positions, detail: str) -> None:
        """Report the input rows at `positions` as one finding, with their count."""
        n = len(positions)
        if n:
            count = "1 row" if n == 1 else f"{n} rows"
            self._findings.append(
                _Finding(
                    code,
                    chrom,
                    np.array([np.min(positions)], dtype=np.int64),
                    np.array([np.max(positions)], dtype=np.int64),
                    False,
                    f"{count}: {detail}",
                )
            )

    def frame(self) -> pd.DataFrame:
        """Return the diagnostics by chromosome, then position or start, then code."""
        findings = self._findings
        sizes = np.array([len(f.start) for f in findings], dtype=np.intp)

        def repeated(values: list[Any], dtype: npt.DTypeLike = object) -> Any:
            return np.repeat(np.array(values, dtype=dtype), sizes)

        chroms = sorted({f.chrom for f in findings}, key=_natural_key)
        chrom_rank = repeated([chroms.index(f.chrom) for f in findings], np.intp)
        code = repeated([f.code for f in findings])
        _, code_rank = np.unique(code, return_inverse=True)
        start = np.concatenate([np.empty(0, np.int64), *(f.start for f in findings)])
        end = np.concatenate([np.empty(0, np.int64), *(f.end for f in findings)])
        per_locus = repeated([f.per_locus for f in findings], np.bool_)
        order = np.lexsort((code_rank, start, chrom_rank))
        per_locus = per_locus[order]
        return pd.DataFrame(
            {
                "code": pd.array(code[order], dtype="str"),
                "chrom": pd.array(repeated([f.chrom for f in findings])[order], "str"),
                "pos": pd.arrays.IntegerArray(start[order], ~per_locus),
                "start": pd.arrays.IntegerArray(start[order], per_locus),
                "end": pd.arrays.IntegerArray(end[order], per_locus),
                "detail": pd.array(
                    repeated([f.detail for f in findings])[order], "str"
                ),
            }
        )


def encode(
    frame: pd.DataFrame,
    *,
    phased: pd.Series[bool] | None = None,
    phase_set: pd.Series[int] | None = None,
    genotype_conflicted: pd.Series[bool] | None = None,
    sex: Sex | None,
) -> tuple[Mapping[str, _Calls], pd.DataFrame]:
    """Return one person's blocks and the diagnostics of normalizing their calls.

    Each Series is aligned with `frame` row by row and has the frame's index.

    Parameters
    ----------
    frame : pandas.DataFrame
        A normalized snps frame on GRCh37: rsIDs as the index, and columns "chrom",
        "pos" and "genotype".
    phased : pandas.Series of bool, optional
        Whether each row is phased: its first allele is on the same chromosome copy
        as the first allele of the other phased rows in its phase set. By default,
        no row is phased.
    phase_set : pandas.Series of int, optional
        The phase set of each row; NA is the chromosome's implicit set. By default,
        each chromosome's phased rows form one set.
    genotype_conflicted : pandas.Series of bool, optional
        Which rows snps nulled because their calls conflicted. By default, none.
    sex : Sex or None
        The person's sex; without one, X is excluded.

    Returns
    -------
    blocks : mapping of str to _Calls
        A read-only block per chromosome with a kept row, 1 to 22 then X.
    diagnostics : pandas.DataFrame
        A row per finding, with a RangeIndex and columns "code", "chrom", "pos",
        "start", "end" and "detail".

    Raises
    ------
    TypeError
        If the index of `frame` is a MultiIndex.
    ValueError
        If the frame, a Series or the sex is malformed.
    """
    sex = _sex(sex)
    for column in ("chrom", "pos", "genotype"):
        count = frame.columns.tolist().count(column)
        if count != 1:
            raise ValueError(f"frame must have one {column!r} column, not {count}")
    if isinstance(frame.index, pd.MultiIndex):
        raise TypeError("The index of frame must hold rsIDs, not a MultiIndex")
    _check_text(frame.index, "The index of frame", null=False)
    _check_text(frame["chrom"], "Column 'chrom'", null=False)
    _check_text(frame["genotype"], "Column 'genotype'", null=True)
    chrom_code, chroms = _row_chromosomes(frame["chrom"])
    name_code, identifiers = _identifiers(frame.index)
    pieces = _split(identifiers)
    genotype_code, domains = _genotype_domains(frame["genotype"])
    rows = _Rows(
        _row_positions(frame["pos"]),
        genotype_code,
        _flags(phased, frame, "phased"),
        _flags(genotype_conflicted, frame, "genotype_conflicted"),
        _set_identities(phase_set, frame),
        name_code,
    )

    natural = sorted(range(len(chroms)), key=lambda i: _natural_key(chroms[i]))
    rank = np.empty(len(chroms), dtype=np.intp)
    rank[natural] = np.arange(len(chroms))
    row_rank = rank[chrom_code]
    rows = rows.take(np.lexsort((rows.pos, row_rank)))
    ends = np.cumsum(np.bincount(row_rank, minlength=len(chroms)))

    report = _Report()
    blocks: dict[str, _Calls] = {}
    for k, i in enumerate(natural):
        chrom = chroms[i]
        start = int(ends[k - 1]) if k else 0
        kept = _kept_rows(chrom, rows.take(slice(start, int(ends[k]))), sex, report)
        if kept is not None and len(kept[0].pos):
            blocks[chrom] = _block(chrom, kept[1], kept[0], domains, pieces, report)
    return MappingProxyType(blocks), report.frame()


def _sex(value: Sex | str | None) -> Sex | None:
    """Return a given sex as a Sex, or None."""
    return None if value is None else Sex(value)


def _natural_key(chrom: str) -> tuple[int, str]:
    """Return the sort key of a chromosome name in natural order."""
    return (_NATURAL.get(chrom, len(_NATURAL)), chrom)


def _copies(chrom: str, sex: Sex | None) -> int | None:
    """Return a person's copies of a modeled chromosome, or None for X of no sex."""
    if chrom != "X":
        return 2
    if sex is None:
        return None
    return 1 if sex == Sex.MALE else 2


def _check_text(values: pd.Series | pd.Index, label: str, *, null: bool) -> None:
    """Check that every value is a str, or null where `null` allows, before hashing."""
    if len(values) == 0 or (null and bool(values.isna().all())):
        return
    if not null:
        _check_present(values, label)
    dtype = values.dtype
    if isinstance(dtype, pd.CategoricalDtype):
        _check_text(dtype.categories, label, null=False)
    elif not isinstance(dtype, pd.StringDtype):
        found = pd.api.types.infer_dtype(values, skipna=True)
        if dtype != np.dtype(object) or found not in ("string", "empty"):
            nulls = " or nulls" if null else ""
            raise ValueError(f"{label} must hold str values{nulls}, not {found} values")


def _check_present(values: pd.Series | pd.Index, label: str) -> None:
    """Check that no value is missing."""
    if values.hasnans:
        raise ValueError(f"{label} has missing values")


def _row_chromosomes(column: pd.Series) -> tuple[_Indices, list[str]]:
    """Return the chromosome of each row of an :func:`encode` frame, and the names."""
    codes, uniques = pd.factorize(column)
    return codes, list(uniques)


def _identifiers(index: pd.Index) -> tuple[_Indices, _Objects]:
    """Return the identifier of each row of an :func:`encode` frame, and the text."""
    codes, uniques = pd.factorize(index)
    text = np.asarray(uniques, dtype=object)
    if (text == "").any():
        raise ValueError("Each identifier in the index of frame must be non-empty")
    return codes, text


def _genotype_domains(column: pd.Series) -> tuple[_Indices, Mapping[int, _Domain]]:
    """Return each row's genotype code and what each code allows, by copies."""
    codes, uniques = pd.factorize(column)
    # A null genotype has code -1, which indexes the last entry.
    values: list[str | None] = [*uniques, None]
    domains = {}
    for copies in _UNRESTRICTED:
        rules = [_genotype_rule(value, copies) for value in values]
        domains[copies] = _Domain(
            np.array([rule[0] for rule in rules], dtype=np.uint32),
            np.array([rule[1] for rule in rules], dtype=np.uint32),
            np.array(
                [
                    -1 if rule[2] is None else list(_ROW_CODES).index(rule[2])
                    for rule in rules
                ],
                dtype=np.int8,
            ),
        )
    return codes, MappingProxyType(domains)


def _genotype_rule(value: str | None, copies: int) -> tuple[int, int, str | None]:
    """Return what a genotype allows, unordered and ordered, and its diagnostic code.

    None stands for a null genotype. A no-call allows every pair of its copies.
    """
    any_pair = _UNRESTRICTED[copies]
    if value is None:
        return any_pair, any_pair, None
    if not 1 <= len(value) <= 2 or any(character not in _BASES for character in value):
        return any_pair, any_pair, "non_base_genotype"
    bases = [_BASES.index(character) for character in value]
    if copies == 1:
        if len(set(bases)) != 1:
            return any_pair, any_pair, "ploidy_mismatch"
        one = _pair(bases[0], _NO_COPY)
        return one, one, None
    if len(bases) != 2:
        return any_pair, any_pair, "ploidy_mismatch"
    a, b = bases
    return _pair(a, b) | _pair(b, a), _pair(a, b), None


def _row_positions(column: pd.Series) -> _Positions:
    """Return the positions of an :func:`encode` frame as a new uint32 array."""
    if len(column) == 0:
        return np.empty(0, dtype=np.uint32)
    label = "Column 'pos'"
    dtype = column.dtype
    if pd.api.types.is_bool_dtype(dtype) or not (
        pd.api.types.is_integer_dtype(dtype) or pd.api.types.is_float_dtype(dtype)
    ):
        raise ValueError(f"{label} must hold whole numbers, not {dtype}")
    _check_present(column, label)
    values = _grch37.check_whole_numbers(column.to_numpy(), label, entry="row")
    outside = (values < 1) | (values > _MAX_POS)
    if outside.any():
        i = int(np.argmax(outside))
        raise ValueError(
            f"{label} must be from 1 to {_MAX_POS:,}; {values[i]} at row {i} is not"
        )
    return values.astype(np.uint32)


def _check_aligned(series: pd.Series, frame: pd.DataFrame, label: str) -> None:
    """Check that `series` is aligned with `frame` row by row."""
    if len(series) != len(frame) or not series.index.equals(frame.index):
        raise ValueError(
            f"{label} must be aligned with frame row by row: it needs the frame's "
            "length and an equal index"
        )


def _flags(series: pd.Series[bool] | None, frame: pd.DataFrame, label: str) -> _Bools:
    """Return a flag per row of an :func:`encode` frame; none is set by default."""
    if series is None:
        return np.zeros(len(frame), dtype=np.bool_)
    _check_aligned(series, frame, label)
    dtype = series.dtype
    if not (
        dtype == np.dtype(np.bool_)
        or (isinstance(dtype, pd.BooleanDtype) and not series.hasnans)
    ):
        raise ValueError(f"{label} must have dtype bool, or boolean with no NA")
    return series.to_numpy(dtype=np.bool_)


def _set_identities(series: pd.Series[int] | None, frame: pd.DataFrame) -> _Indices:
    """Return which phase set each row of an :func:`encode` frame names.

    The numbers only identify sets, so a code stands for each; NA is -1.
    """
    if series is None:
        return np.full(len(frame), -1, dtype=np.intp)
    _check_aligned(series, frame, "phase_set")
    if not pd.api.types.is_integer_dtype(series.dtype):
        raise ValueError(f"phase_set must have an integer dtype, not {series.dtype}")
    if bool((series < 0).any()):
        raise ValueError("phase_set must not hold negative numbers")
    codes, _ = pd.factorize(series)
    return codes


def _kept_rows(
    chrom: str, rows: _Rows, sex: Sex | None, report: _Report
) -> tuple[_Rows, int] | None:
    """Return the rows of a chromosome that are kept, and the copies they are on.

    The rows are ordered by position; those left out are reported.
    """
    length = _grch37.LENGTHS.get(chrom)
    if length is None:
        report.rows(
            "excluded_region", chrom, rows.pos, f"chromosome {chrom} is not modeled"
        )
        return None
    valid = int(np.searchsorted(rows.pos, length, side="right"))
    report.loci(
        "invalid_position",
        chrom,
        np.unique(rows.pos[valid:]),
        f"Beyond the end of chromosome {chrom} ({length:,} bases); the rows here are "
        "left out",
    )
    rows = rows.take(slice(0, valid))
    copies = _copies(chrom, sex)
    if copies is None:
        report.rows("excluded_region", chrom, rows.pos, "X needs a known sex")
        return None
    if chrom == "X":
        outside = np.ones(len(rows.pos), dtype=np.bool_)
        for region, (first, last) in enumerate((_grch37.PAR1, _grch37.PAR2), start=1):
            inside = (rows.pos >= first) & (rows.pos <= last)
            report.rows(
                "excluded_region",
                chrom,
                rows.pos[inside],
                f"pseudoautosomal region {region} ({first:,} to {last:,}) is not "
                "modeled",
            )
            outside &= ~inside
        rows = rows.take(outside)
    return rows, copies


def _block(
    chrom: str,
    copies: int,
    rows: _Rows,
    domains: Mapping[int, _Domain],
    pieces: _Pieces,
    report: _Report,
) -> _Calls:
    """Return the block of one chromosome's kept rows, ordered by position."""
    new = np.r_[True, rows.pos[1:] != rows.pos[:-1]]
    starts = np.flatnonzero(new)
    locus = np.cumsum(new) - 1
    pos = rows.pos[starts]
    n = len(starts)
    domain = domains[copies]
    any_pair = _UNRESTRICTED[copies]

    # A no-call allows every pair, so it never narrows the intersection.
    unordered = domain.unordered[rows.genotype]
    common = np.bitwise_and.reduceat(unordered, starts)
    conflict = np.logical_or.reduceat(rows.flagged, starts) | (common == 0)

    kept, dropped, order, phase_set = _reconcile_asserted_phase(
        rows, starts, copies, domain, unordered, common, conflict
    )
    allowed = np.where(conflict, any_pair, np.where(kept, order, common))
    kinds = np.zeros((n, len(_ATTRIBUTES)), dtype=np.uint8)
    kinds[allowed != any_pair, _GENOTYPE] = _OBSERVED
    kinds[kept, _PHASE] = _ASSERTED

    report.loci(
        "position_conflict",
        chrom,
        pos[conflict],
        "A row here is flagged as conflicted or the calls share no pair; the locus is "
        "a no-call",
    )
    report.loci(
        "phase_conflict",
        chrom,
        pos[dropped],
        "The phased rows here disagree on order or phase set; the genotype is kept "
        "unordered",
    )
    code = domain.code[rows.genotype]
    copies_text = "one copy" if copies == 1 else "two copies"
    for i, (row_code, detail) in enumerate(_ROW_CODES.items()):
        report.rows(
            row_code, chrom, rows.pos[code == i], detail.format(copies=copies_text)
        )

    rsid = _locus_names(chrom, rows.pos, locus, rows.name, pieces, n)
    return _Calls(
        Loci(chrom, pos, rsid),
        allowed,
        phase_set,
        kinds,
    )


def _reconcile_asserted_phase(
    rows: _Rows,
    starts: _Indices,
    copies: int,
    domain: _Domain,
    unordered: _Masks,
    common: _Masks,
    conflict: _Bools,
) -> tuple[_Bools, _Bools, _Masks, _Indices]:
    """Return which loci keep or drop an asserted order, and their masks and sets.

    A locus keeps its order only if all its phased rows name one set and agree on order.
    """
    any_pair = _UNRESTRICTED[copies]
    if copies == 2:
        phased = rows.phased & (unordered != any_pair)
    else:
        phased = np.zeros(len(rows.pos), dtype=np.bool_)
    tried = np.logical_or.reduceat(phased, starts) & ~conflict
    order = np.bitwise_and.reduceat(
        np.where(phased, domain.ordered[rows.genotype], any_pair), starts
    )
    order &= common
    first_set = np.minimum.reduceat(
        np.where(phased, rows.phase_set, np.iinfo(np.intp).max), starts
    )
    last_set = np.maximum.reduceat(np.where(phased, rows.phase_set, -2), starts)
    # Equal least and greatest set identities: the phased rows all name one set.
    kept = tried & (first_set == last_set) & (order != 0)
    phase_set = np.full(len(starts), -1, dtype=np.intp)
    sets, first_at, which = np.unique(
        first_set[kept], return_index=True, return_inverse=True
    )
    number = np.empty(len(sets), dtype=np.intp)
    number[np.argsort(first_at)] = np.arange(len(sets))
    phase_set[kept] = number[which]
    return kept, tried & ~kept, order, phase_set


def _split(identifiers: _Objects) -> _Pieces:
    """Return the non-empty ";"-separated names in each of the distinct identifiers."""
    several = pd.Series(identifiers, dtype=object).str.contains(";", regex=False)
    single = np.flatnonzero(~several.to_numpy(dtype=np.bool_))
    split = np.flatnonzero(several.to_numpy(dtype=np.bool_))
    parts = (
        pd.Series(identifiers[split], index=split, dtype=object)
        .str.split(";")
        .explode()
    )
    parts = parts[parts != ""]
    owner = np.concatenate([single, parts.index.to_numpy(dtype=np.intp)])
    found = np.concatenate([identifiers[single], parts.to_numpy(dtype=object)])
    name, text = pd.factorize(found[np.argsort(owner, kind="stable")])
    length = np.bincount(owner, minlength=len(identifiers))
    distinct = np.asarray(text, dtype=object)
    colon = pd.Series(distinct, dtype=object).str.contains(":", regex=False)
    return _Pieces(
        np.cumsum(length) - length, length, name, distinct, colon.to_numpy(np.bool_)
    )


def _locus_names(
    chrom: str,
    pos: _Positions,
    locus: _Indices,
    ident: _Indices,
    pieces: _Pieces,
    n: int,
) -> _Objects:
    """Return the names of `n` loci on `chrom`, without their own placeholders.

    Entry i gives locus ``locus[i]``, at position ``pos[i]``, identifier ``ident[i]``.
    """
    length = pieces.length[ident]
    entry = np.repeat(np.arange(len(ident)), length)
    within = np.arange(len(entry)) - np.repeat(np.cumsum(length) - length, length)
    name = pieces.name[pieces.start[ident][entry] + within]
    own = np.zeros(len(name), dtype=np.bool_)
    colon = pieces.colon[name]
    if colon.any():
        placeholder = f"{chrom}:" + pos[entry[colon]].astype(str).astype(object)
        own[colon] = pieces.text[name[colon]] == placeholder
    return _joined_names(locus[entry[~own]], name[~own], pieces.text, n)


def _joined_names(locus: _Indices, name: _Indices, text: _Objects, n: int) -> _Objects:
    """Return each of `n` loci's distinct names, sorted and joined by ";", or None.

    Locus ``locus[i]`` has the name ``text[name[i]]``.
    """
    joined = np.full(n, None, dtype=object)
    if len(locus) == 0:
        return joined
    m = len(text)
    pairs = np.unique(locus.astype(np.int64) * m + name)
    at, name = np.divmod(pairs, m)
    starts = np.flatnonzero(np.r_[True, at[1:] != at[:-1]])
    ends = np.r_[starts[1:], len(at)]
    one = ends - starts == 1
    joined[at[starts[one]]] = text[name[starts[one]]]
    for s, e in zip(starts[~one].tolist(), ends[~one].tolist(), strict=True):
        joined[at[s]] = ";".join(sorted(text[name[s:e]].tolist()))
    return joined


def align(
    people: Mapping[str, Mapping[str, _Calls]],
    *,
    sexes: Mapping[str, Sex | None],
) -> Mapping[str, Mapping[str, _Calls]]:
    """Return a group's blocks on one shared locus axis per chromosome.

    Parameters
    ----------
    people : mapping of str to mapping of str to _Calls
        Each tested person's blocks, by chromosome.
    sexes : mapping of str to Sex or None
        Every participant's sex, in the order of the result.

    Returns
    -------
    mapping of str to mapping of str to _Calls
        Read-only, by name, then by chromosome in natural order; the blocks on a
        chromosome share one Loci.

    Raises
    ------
    ValueError
        If a sex is not valid, or a person or a block does not fit `sexes` or its key.
    """
    sexes = {name: _sex(sex) for name, sex in sexes.items()}
    for name, blocks in people.items():
        if name not in sexes:
            raise ValueError(
                f"{name!r} has blocks but no sex; sexes must name everyone in people"
            )
        for chrom, calls in blocks.items():
            if calls.loci.chrom != chrom:
                raise ValueError(
                    f"The block of {name!r} keyed {chrom!r} is on chromosome "
                    f"{calls.loci.chrom!r}"
                )
            if chrom == "X":
                copies = _copies(chrom, sexes[name])
                if copies is None:
                    raise ValueError(f"{name!r} has an X block but no known sex")
                given_copies = int(_copies_of(calls.allowed[:1])[0])
                if given_copies != copies:
                    raise ValueError(
                        f"The X block of {name!r} has {given_copies} copies, but a "
                        f"{sexes[name]} has {copies}"
                    )
    chroms = sorted(
        {chrom for blocks in people.values() for chrom in blocks}, key=_natural_key
    )
    aligned: dict[str, dict[str, _Calls]] = {name: {} for name in sexes}
    for chrom in chroms:
        given = {
            name: blocks[chrom] for name, blocks in people.items() if chrom in blocks
        }
        loci = _union(chrom, [calls.loci for calls in given.values()])
        for name, sex in sexes.items():
            copies = _copies(chrom, sex)
            if copies is not None:
                aligned[name][chrom] = _placed(loci, copies, given.get(name))
    return MappingProxyType(
        {name: MappingProxyType(blocks) for name, blocks in aligned.items()}
    )


def _union(chrom: str, axes: list[Loci]) -> Loci:
    """Return the sorted union of locus axes, with the names at each position merged."""
    pos = np.concatenate([axis.pos for axis in axes])
    rsid = np.concatenate([axis.rsid for axis in axes])
    union, at = np.unique(pos, return_inverse=True)
    named = np.flatnonzero(pd.notna(rsid))
    ident, identifiers = pd.factorize(rsid[named])
    names = _locus_names(
        chrom,
        pos[named],
        at[named],
        ident,
        _split(np.asarray(identifiers, dtype=object)),
        len(union),
    )
    return Loci(chrom, union, names)


def _placed(loci: Loci, copies: int, calls: _Calls | None) -> _Calls:
    """Return a block on `loci` holding `calls`, if given, and nothing elsewhere."""
    n = len(loci.pos)
    allowed = np.full(n, _UNRESTRICTED[copies], dtype=np.uint32)
    phase_set = np.full(n, -1, dtype=np.int32)
    kinds = np.zeros((n, len(_ATTRIBUTES)), dtype=np.uint8)
    if calls is not None:
        at = np.searchsorted(loci.pos, calls.loci.pos)
        allowed[at] = calls.allowed
        phase_set[at] = calls.phase_set
        kinds[at] = calls.kinds
    return _Calls(loci, allowed, phase_set, kinds)
