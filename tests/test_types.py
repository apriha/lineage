"""Tests of the loci of a chromosome, on invented positions and names."""

from __future__ import annotations

import copy
import dataclasses
import pickle
import warnings

import numpy as np
import pytest

from lineage.types import Loci


def test_loci():
    loci = Loci("1", [10, 20, 30], ["rs1", None, "rs2;rs3"])
    assert loci.chrom == "1"
    assert loci.pos.dtype == np.uint32
    assert loci.pos.tolist() == [10, 20, 30]
    assert loci.rsid.dtype == object
    assert loci.rsid.tolist() == ["rs1", None, "rs2;rs3"]


def test_rsid_left_out():
    loci = Loci("2", [10, 20, 30])
    assert loci.rsid.dtype == object
    assert loci.rsid.tolist() == [None, None, None]


@pytest.mark.parametrize("rsid", [None, []])
def test_empty(rsid):
    loci = Loci("3", [], rsid)
    assert loci.pos.shape == (0,)
    assert loci.pos.dtype == np.uint32
    assert loci.rsid.shape == (0,)
    assert loci.rsid.dtype == object


@pytest.mark.parametrize(
    "pos",
    [
        [10, 20, 30],
        (10, 20, 30),
        np.array([10, 20, 30], dtype=np.int16),
        np.array([10, 20, 30], dtype=np.int64),
        np.array([10, 20, 30], dtype=np.uint16),
        np.array([10, 20, 30], dtype=np.uint32),
        np.array([10, 20, 30], dtype=np.uint64),
        np.array([10.0, 20.0, 30.0], dtype=np.float32),
        [10.0, 20.0, 30.0],
    ],
    ids=[
        "list",
        "tuple",
        "int16",
        "int64",
        "uint16",
        "uint32",
        "uint64",
        "float32",
        "floats",
    ],
)
def test_position_types(pos):
    loci = Loci("4", pos)
    assert loci.pos.dtype == np.uint32
    assert loci.pos.tolist() == [10, 20, 30]


def test_rsid_from_a_str_array():
    loci = Loci("5", [10, 20], np.array(["rs1", "rs2"]))
    assert loci.rsid.dtype == object
    assert loci.rsid.tolist() == ["rs1", "rs2"]


@pytest.mark.parametrize(
    ("pos", "accepted"),
    [
        (1, True),
        (60_000, True),
        (60_001, False),
        (2_699_520, False),
        (2_699_521, True),
        (154_931_043, True),
        (154_931_044, False),
        (155_260_560, False),
        (155_260_561, True),
        (155_270_560, True),
        (155_270_561, False),
    ],
)
def test_x_boundaries(pos, accepted):
    if accepted:
        assert Loci("X", [pos]).pos.tolist() == [pos]
    else:
        with pytest.raises(ValueError, match="outside the pseudoautosomal regions"):
            Loci("X", [pos])


@pytest.mark.parametrize(
    ("chrom", "length"),
    [("1", 249_250_621), ("21", 48_129_895), ("22", 51_304_566)],
)
def test_last_base(chrom, length):
    assert Loci(chrom, [1, length]).pos.tolist() == [1, length]
    with pytest.raises(ValueError, match=f"must be 1 to {length:,}"):
        Loci(chrom, [1, length + 1])


@pytest.mark.parametrize(
    ("chrom", "pos", "dtype", "accepted"),
    [
        ("1", 249_250_624, np.float32, False),
        ("22", 51_304_568, np.float32, False),
        ("X", 154_931_040, np.float32, True),
        ("X", 60_000, np.float16, True),
        ("1", 60_000, np.float16, True),
    ],
)
def test_float32_and_float16_boundaries(chrom, pos, dtype, accepted):
    values = np.array([pos], dtype=dtype)
    assert values[0] == pos
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        if accepted:
            assert Loci(chrom, values).pos.tolist() == [pos]
        else:
            with pytest.raises(ValueError, match=f"; {pos}.0 at index 0 is not"):
                Loci(chrom, values)


@pytest.mark.parametrize(
    ("pos", "message"),
    [
        ([0], "must be 1 to 249,250,621; 0 at index 0 is not"),
        ([-5], "must be 1 to 249,250,621; -5 at index 0 is not"),
        ([10.0, 20.5], "must be finite whole numbers; 20.5 at index 1 is not"),
        ([10.0, np.nan], "must be finite whole numbers; nan at index 1 is not"),
        ([np.inf], "must be finite whole numbers; inf at index 0 is not"),
        ([True, False], "must be integers, not bool"),
        (np.array([True, False]), "must be integers, not bool"),
        ([10, True], "must be integers, not bool"),
        ((10, np.True_), "must be integers, not bool"),
        ([[10, True]], "must be integers, not bool"),
        (
            np.array([2**32 + 100], dtype=np.int64),
            "must be 1 to 249,250,621; 4294967396 at index 0 is not",
        ),
        (["10", "20"], "must be integers, not <U2"),
        (np.array([10, 20], dtype=object), "must be integers, not object"),
        ([10 + 0j], "must be integers, not complex128"),
        ([10, 20, 20], "must be strictly increasing; 20 at index 2 does not follow 20"),
        ([10, 30, 20], "must be strictly increasing; 20 at index 2 does not follow 30"),
        (
            np.array([30, 20], dtype=np.uint64),
            "must be strictly increasing; 20 at index 1",
        ),
        ([[10, 20], [30, 40]], "must be a 1-D array, not 2-D"),
        (10, "must be a 1-D array, not 0-D"),
    ],
)
def test_invalid_positions(pos, message):
    with pytest.raises(ValueError, match=f"Positions on chromosome '1' {message}"):
        Loci("1", pos)


@pytest.mark.parametrize("chrom", ["Y", "MT", "PAR", "23", "chr1", "x", ""])
def test_unknown_chromosome(chrom):
    with pytest.raises(ValueError, match='must be "1" to "22" or "X"'):
        Loci(chrom, [10])


@pytest.mark.parametrize("chrom", [1, None, b"1"])
def test_chromosome_not_a_str(chrom):
    with pytest.raises(TypeError, match="A chromosome must be a str"):
        Loci(chrom, [10])


@pytest.mark.parametrize(
    ("rsid", "message"),
    [
        (["rs1"], "1-D array of 2 names, one per position, not shape \\(1,\\)"),
        (["rs1", "rs2", "rs3"], "not shape \\(3,\\)"),
        ([["rs1", "rs2"]], "not shape \\(1, 2\\)"),
        ("rs1", "not shape \\(\\)"),
        (["rs1", 2], "non-empty str or None; 2 at index 1 is not"),
        (["", "rs2"], "non-empty str or None; '' at index 0 is not"),
        (["rs1", np.nan], "non-empty str or None; nan at index 1 is not"),
    ],
)
def test_invalid_rsid(rsid, message):
    with pytest.raises(ValueError, match=message):
        Loci("1", [10, 20], rsid)


def test_caller_arrays_are_copied():
    pos = np.array([10, 20, 30], dtype=np.uint32)
    rsid = np.array(["rs1", "rs2", "rs3"], dtype=object)
    loci = Loci("1", pos, rsid)
    pos[0] = 15
    rsid[0] = "rs9"
    assert loci.pos.tolist() == [10, 20, 30]
    assert loci.rsid.tolist() == ["rs1", "rs2", "rs3"]
    assert not np.shares_memory(loci.pos, pos)
    assert not np.shares_memory(loci.rsid, rsid)


@pytest.mark.parametrize("rsid", [None, ["rs1", "rs2", "rs3"]])
def test_arrays_are_read_only(rsid):
    loci = Loci("1", [10, 20, 30], rsid)
    for array, value in ((loci.pos, 15), (loci.rsid, "rs9")):
        with pytest.raises(ValueError, match="read-only"):
            array[0] = value
        if array.base is not None:
            with pytest.raises(ValueError, match="read-only"):
                array.base[0] = value
    assert loci.pos.tolist() == [10, 20, 30]


@pytest.mark.parametrize(
    "copied",
    [copy.copy, copy.deepcopy, lambda loci: pickle.loads(pickle.dumps(loci))],
    ids=["copy", "deepcopy", "pickle"],
)
def test_copies_stay_read_only(copied):
    loci = copied(Loci("X", [10, 20], ["rs1", None]))
    assert loci.chrom == "X"
    assert loci.pos.dtype == np.uint32
    assert loci.pos.tolist() == [10, 20]
    assert loci.rsid.tolist() == ["rs1", None]
    assert not loci.pos.flags.writeable
    assert not loci.rsid.flags.writeable


def test_frozen():
    loci = Loci("1", [10])
    with pytest.raises(dataclasses.FrozenInstanceError):
        loci.chrom = "2"
    with pytest.raises(dataclasses.FrozenInstanceError):
        loci.pos = np.array([20], dtype=np.uint32)
