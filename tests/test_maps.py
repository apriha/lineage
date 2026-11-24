"""Tests of genetic maps, on invented tables and uniform maps."""

from __future__ import annotations

import copy
import dataclasses
import pickle

import numpy as np
import pandas as pd
import pytest

from lineage.maps import GeneticMap

REGIONS = (
    *(str(i) for i in range(1, 23)),
    "X_par1",
    "X",
    "X_par2",
)

SHA256 = "0123456789abcdef" * 4


def table(pos, rate=None, cm=None):
    """Return a table; rates default to 1 and genetic positions to 0, 1, 2, ..."""
    n = len(pos)
    return pd.DataFrame(
        {
            "pos": pos,
            "rate": [1.0] * n if rate is None else rate,
            "map": [float(i) for i in range(n)] if cm is None else cm,
        }
    )


def one_region(frame, region="1"):
    return GeneticMap.from_frames("test", {region: frame})


@pytest.fixture
def small_map():
    return GeneticMap.from_frames(
        "test",
        {"1": table([100, 200, 400], rate=[1.0, 1.0, 1.0], cm=[0.0, 1.0, 3.0])},
        source="maps/test.txt",
        sha256=SHA256,
    )


def test_uniform_regions():
    genetic_map = GeneticMap.uniform()
    assert genetic_map.regions == REGIONS
    assert tuple(genetic_map.tables) == REGIONS
    assert len(genetic_map.tables) == 25
    assert genetic_map.name == "uniform"
    assert genetic_map.build == 37
    assert genetic_map.source is None
    assert genetic_map.sha256 is None


@pytest.mark.parametrize(
    ("region", "first", "last", "cm"),
    [
        ("1", 1, 249_250_621, 249.25062),
        ("22", 1, 51_304_566, 51.304565),
        ("X_par1", 60_001, 2_699_520, 2.639519),
        ("X", 2_699_521, 154_931_043, 152.231522),
        ("X_par2", 154_931_044, 155_260_560, 0.329516),
    ],
)
def test_uniform_tables(region, first, last, cm):
    frame = GeneticMap.uniform().tables[region]
    assert frame["pos"].tolist() == [first, last]
    assert frame["rate"].tolist() == [1.0, 1.0]
    assert frame["map"].iloc[0] == 0.0
    assert frame["map"].iloc[-1] == pytest.approx(cm, abs=1e-12)


@pytest.mark.parametrize(
    ("chrom", "length"),
    [
        ("1", 249_250_621),
        ("2", 243_199_373),
        ("3", 198_022_430),
        ("4", 191_154_276),
        ("5", 180_915_260),
        ("6", 171_115_067),
        ("7", 159_138_663),
        ("8", 146_364_022),
        ("9", 141_213_431),
        ("10", 135_534_747),
        ("11", 135_006_516),
        ("12", 133_851_895),
        ("13", 115_169_878),
        ("14", 107_349_540),
        ("15", 102_531_392),
        ("16", 90_354_753),
        ("17", 81_195_210),
        ("18", 78_077_248),
        ("19", 59_128_983),
        ("20", 63_025_520),
        ("21", 48_129_895),
        ("22", 51_304_566),
    ],
)
def test_autosome_lengths(chrom, length):
    assert GeneticMap.uniform().tables[chrom]["pos"].tolist() == [1, length]
    one_region(table([length]), chrom)
    with pytest.raises(ValueError, match=f"must be 1 to {length:,};"):
        one_region(table([length + 1]), chrom)


def test_uniform_rate():
    genetic_map = GeneticMap.uniform(rate=2)
    frame = genetic_map.tables["1"]
    assert frame["rate"].tolist() == [2.0, 2.0]
    assert frame["map"].iloc[-1] == pytest.approx(498.50124, abs=1e-12)
    assert genetic_map.cm("1", 125_000_001).item() == pytest.approx(250.0)


def test_uniform_rate_zero():
    genetic_map = GeneticMap.uniform(rate=0)
    for region in genetic_map.regions:
        assert genetic_map.tables[region]["map"].tolist() == [0.0, 0.0]


@pytest.mark.parametrize(
    ("chromosomes", "regions"),
    [
        (["X", "2"], ("2", "X")),
        (("X_par2", "X_par1", "22"), ("22", "X_par1", "X_par2")),
        (["1"], ("1",)),
    ],
)
def test_uniform_chromosomes(chromosomes, regions):
    assert GeneticMap.uniform(chromosomes=chromosomes).regions == regions


@pytest.mark.parametrize(
    ("kwargs", "error", "message"),
    [
        ({"rate": True}, TypeError, "rate must be an int or float, not bool"),
        ({"rate": "1"}, TypeError, "rate must be an int or float, not str"),
        ({"rate": None}, TypeError, "rate must be an int or float, not NoneType"),
        ({"rate": -0.5}, ValueError, "rate must be finite and at least 0, not -0.5"),
        ({"rate": np.nan}, ValueError, "rate must be finite and at least 0, not nan"),
        ({"rate": np.inf}, ValueError, "rate must be finite and at least 0, not inf"),
        (
            {"rate": 1e307},
            ValueError,
            "Column 'map' of the table of region '1' has values that are not finite",
        ),
        (
            {"chromosomes": "1"},
            TypeError,
            "must be a sequence of region names, not str",
        ),
        (
            {"chromosomes": {"1"}},
            TypeError,
            "must be a sequence of region names, not set",
        ),
        ({"chromosomes": []}, ValueError, "chromosomes is empty"),
        ({"chromosomes": ["1", "2", "1"]}, ValueError, "chromosomes repeats a region"),
        ({"chromosomes": ["1", "Y"]}, ValueError, "'Y' is not a region"),
        ({"chromosomes": ["chrX"]}, ValueError, "'chrX' is not a region"),
        ({"chromosomes": [1]}, ValueError, "1 is not a region"),
    ],
)
def test_invalid_uniform(kwargs, error, message):
    with pytest.raises(error, match=message):
        GeneticMap.uniform(**kwargs)


def test_cm(small_map):
    result = small_map.cm("1", [50, 100, 150, 300, 400, 500])
    assert result.dtype == np.float64
    assert result.tolist() == [0.0, 0.0, 0.5, 2.0, 3.0, 3.0]


def test_cm_keeps_shape(small_map):
    scalar = small_map.cm("1", 150)
    assert scalar.shape == ()
    assert scalar.dtype == np.float64
    assert scalar.item() == 0.5
    grid = small_map.cm("1", np.array([[50, 150], [300, 500]], dtype=np.uint32))
    assert grid.shape == (2, 2)
    assert grid.tolist() == [[0.0, 0.5], [2.0, 3.0]]
    assert small_map.cm("1", []).shape == (0,)
    assert small_map.cm("1", [150.5, 250.0]).tolist() == pytest.approx([0.505, 1.5])


def test_cm_of_one_row():
    genetic_map = one_region(table([1000], cm=[5.0]))
    assert genetic_map.cm("1", [1, 1000, 10**9]).tolist() == [5.0, 5.0, 5.0]


@pytest.mark.parametrize("region", ["2", "X", "Y", "chr1", 1, None])
def test_cm_of_missing_region(small_map, region):
    with pytest.raises(KeyError) as info:
        small_map.cm(region, [100])
    assert info.value.args[0] == (
        f"The map has no region {region!r}; use one of its regions ('1')"
    )


def test_cm_lists_regions():
    genetic_map = GeneticMap.uniform(chromosomes=["X", "2", "1"])
    with pytest.raises(KeyError, match=r"use one of its regions \('1', '2', 'X'\)"):
        genetic_map.cm("X_par1", [100])


@pytest.mark.parametrize(
    ("positions", "message"),
    [
        ([100, np.nan], "must be finite numbers"),
        ([np.inf], "must be finite numbers"),
        ([True], "must be numbers, not bool"),
        (True, "must be numbers, not bool"),
        ([True, 150], "must be numbers, not bool"),
        ((150, np.True_), "must be numbers, not bool"),
        ([[150, 300], [400, False]], "must be numbers, not bool"),
        (np.array([True, False]), "must be numbers, not bool"),
        (["100"], "must be numbers, not <U3"),
        ([100 + 0j], "must be numbers, not complex128"),
        (np.array([100], dtype=object), "must be numbers, not object"),
    ],
)
def test_cm_of_invalid_positions(small_map, positions, message):
    with pytest.raises(ValueError, match=f"Positions {message}"):
        small_map.cm("1", positions)


def test_from_frames_orders_regions():
    keys = ["X", "2", "1", "X_par2", "X_par1"]
    genetic_map = GeneticMap.from_frames(
        "test", {key: GeneticMap.uniform().tables[key] for key in keys}
    )
    assert genetic_map.regions == ("1", "2", "X_par1", "X", "X_par2")
    assert list(genetic_map.tables) == ["1", "2", "X_par1", "X", "X_par2"]


def test_from_frames_fields():
    genetic_map = GeneticMap.from_frames(
        "HapMap2", {"1": table([100])}, source="maps/test.txt", sha256=SHA256.upper()
    )
    assert genetic_map.name == "HapMap2"
    assert genetic_map.build == 37
    assert genetic_map.source == "maps/test.txt"
    assert genetic_map.sha256 == SHA256


def test_direct_construction():
    genetic_map = GeneticMap("test", 37, {"1": table([100, 200])}, sha256=SHA256)
    assert genetic_map.regions == ("1",)
    assert genetic_map.sha256 == SHA256
    with pytest.raises(ValueError, match="'Y' is not a region"):
        GeneticMap("test", 37, {"Y": table([100])})


def test_replace_checks_again():
    genetic_map = GeneticMap.uniform(chromosomes=["1"])
    renamed = dataclasses.replace(genetic_map, name="renamed")
    assert renamed.name == "renamed"
    pd.testing.assert_frame_equal(renamed.tables["1"], genetic_map.tables["1"])
    with pytest.raises(ValueError, match="name is empty"):
        dataclasses.replace(genetic_map, name="")


def test_stored_tables():
    frame = pd.DataFrame(
        {
            "map": [0.0, 1.5],
            "extra": ["a", "b"],
            "pos": [100.0, 300.0],
            "rate": [1, 2],
        },
        index=[7, 3],
    )
    stored = one_region(frame).tables["1"]
    expected = pd.DataFrame(
        {
            "pos": np.array([100, 300], dtype=np.uint32),
            "rate": [1.0, 2.0],
            "map": [0.0, 1.5],
        }
    )
    pd.testing.assert_frame_equal(stored, expected)
    assert isinstance(stored.index, pd.RangeIndex)
    assert stored.columns.tolist() == ["pos", "rate", "map"]


def test_nullable_columns():
    frame = pd.DataFrame(
        {
            "pos": pd.array([100, 200], dtype="Int64"),
            "rate": pd.array([1.0, 1.0], dtype="Float64"),
            "map": pd.array([0.0, 1.0], dtype="Float64"),
        }
    )
    stored = one_region(frame).tables["1"]
    assert stored.dtypes.tolist() == [np.uint32, np.float64, np.float64]
    assert stored["pos"].tolist() == [100, 200]


@pytest.mark.parametrize(
    ("region", "pos", "accepted"),
    [
        ("X_par1", 60_000, False),
        ("X_par1", 60_001, True),
        ("X_par1", 2_699_520, True),
        ("X_par1", 2_699_521, False),
        ("X", 1, True),
        ("X", 60_000, True),
        ("X", 60_001, False),
        ("X", 2_699_520, False),
        ("X", 2_699_521, True),
        ("X", 154_931_043, True),
        ("X", 154_931_044, False),
        ("X", 155_260_560, False),
        ("X", 155_260_561, True),
        ("X", 155_270_560, True),
        ("X", 155_270_561, False),
        ("X_par2", 154_931_043, False),
        ("X_par2", 154_931_044, True),
        ("X_par2", 155_260_560, True),
        ("X_par2", 155_260_561, False),
        ("1", 1, True),
        ("1", 249_250_621, True),
        ("1", 249_250_622, False),
        ("1", 0, False),
    ],
)
def test_position_boundaries(region, pos, accepted):
    if accepted:
        assert one_region(table([pos]), region).tables[region]["pos"].tolist() == [pos]
    else:
        with pytest.raises(ValueError, match=f"region '{region}' must be [0-9,]+ to"):
            one_region(table([pos]), region)


def test_equal_genetic_positions():
    genetic_map = one_region(
        table([100, 200, 300], rate=[0.0, 0.0, 0.5], cm=[1.0, 1.0, 1.0])
    )
    assert genetic_map.tables["1"]["map"].tolist() == [1.0, 1.0, 1.0]
    assert genetic_map.tables["1"]["rate"].tolist() == [0.0, 0.0, 0.5]


@pytest.mark.parametrize(
    ("frame", "message"),
    [
        (table([0, 100]), "Positions in the table of region '1' must be 1 to"),
        (table([0, 100]), "must be 1 to 249,250,621; 0 at index 0 is not"),
        (table([-100]), "must be 1 to 249,250,621; -100 at index 0 is not"),
        (
            table([100.0, 100.5]),
            "must be finite whole numbers; 100.5 at index 1 is not",
        ),
        (
            table(np.array([100, 2**32 + 100], dtype=np.int64)),
            "must be 1 to 249,250,621; 4294967396 at index 1 is not",
        ),
        (
            table([100, 200, 200]),
            "must be strictly increasing; 200 at index 2 does not follow 200",
        ),
        (
            table([100, 300, 200]),
            "must be strictly increasing; 200 at index 2 does not follow 300",
        ),
        (
            table([100, 200], cm=[1.0, 0.5]),
            "genetic positions that decrease; row 1 has 0.5",
        ),
        (table([100, 200], rate=[1.0, -0.1]), "negative rates; row 1 has -0.1"),
        (
            table([100, np.nan]),
            "Positions .* finite whole numbers; nan at index 1 is not",
        ),
        (table([100, 200], cm=[0.0, np.nan]), "Column 'map' .* has missing values"),
        (table([100, 200], rate=[np.nan, 1.0]), "Column 'rate' .* has missing values"),
        (
            table([100, np.inf]),
            "Positions .* finite whole numbers; inf at index 1 is not",
        ),
        (
            table([100, 200], cm=[0.0, np.inf]),
            "Column 'map' .* not finite; row 1 has inf",
        ),
        (
            table([100, 200], rate=[-np.inf, 1.0]),
            "Column 'rate' .* not finite; row 0 has -inf",
        ),
        (
            table(pd.array([100, None], dtype="Int64")),
            "Positions .* finite whole numbers; nan at index 1 is not",
        ),
        (
            table([True, False]),
            "Positions in the table of region '1' must be integers, not bool",
        ),
        (
            table([100, 200], cm=[False, True]),
            "Column 'map' .* must be numeric, not bool",
        ),
        (table([100, 200], rate=[True, True]), "Column 'rate' .* must be numeric"),
        (table(["100", "200"]), "Positions .* must be integers, not object"),
        (table([100, 200], cm=["0", "1"]), "Column 'map' .* must be numeric"),
        (table([100], cm=[1 + 0j]), "Column 'map' .* must be numeric, not complex128"),
        (
            table([100]).drop(columns="rate"),
            "has no 'rate' column; it needs columns 'pos', 'rate', 'map'",
        ),
        (pd.DataFrame({"position": [100]}), "has no 'pos' or 'rate' or 'map' column"),
        (
            pd.concat([table([100]), table([100])[["pos"]]], axis=1),
            "has more than one 'pos' column",
        ),
        (table([]), "region '1' is empty; give at least one row"),
    ],
)
def test_invalid_table(frame, message):
    with pytest.raises(ValueError, match=message):
        one_region(frame)


@pytest.mark.parametrize("region", ["Y", "chrX", "23", "x", "X_PAR1", "", 1, None])
def test_unknown_region(region):
    with pytest.raises(ValueError, match="is not a region of a genetic map") as info:
        GeneticMap.from_frames("test", {"1": table([100]), region: table([100])})
    assert str(info.value).endswith('use "1" to "22", "X_par1", "X" or "X_par2"')


@pytest.mark.parametrize(
    ("tables", "error", "message"),
    [
        ({}, ValueError, "tables is empty"),
        ([("1", table([100]))], TypeError, "tables must be a mapping"),
        (table([100]), TypeError, "tables must be a mapping"),
        ({"1": {"pos": [100], "rate": [1.0], "map": [0.0]}}, TypeError, "DataFrame"),
        ({"1": table([100])["pos"]}, TypeError, "must be a DataFrame, not Series"),
    ],
)
def test_invalid_tables(tables, error, message):
    with pytest.raises(error, match=message):
        GeneticMap.from_frames("test", tables)


@pytest.mark.parametrize(
    ("kwargs", "error", "message"),
    [
        ({"name": 5}, TypeError, "name must be a str, not int"),
        ({"name": ""}, ValueError, "name is empty"),
        ({"build": 36}, ValueError, "must be 37, not 36"),
        ({"build": "37"}, ValueError, "must be 37, not '37'"),
        ({"build": 37.0}, ValueError, "must be 37, not 37.0"),
        ({"build": True}, ValueError, "must be 37, not True"),
        ({"build": np.int64(37)}, ValueError, "must be 37"),
        ({"source": 5}, TypeError, "source must be a str or None, not int"),
        ({"sha256": SHA256[:-1]}, ValueError, "64 hexadecimal characters"),
        ({"sha256": SHA256 + "0"}, ValueError, "64 hexadecimal characters"),
        ({"sha256": "g" * 64}, ValueError, "64 hexadecimal characters"),
        ({"sha256": " " + SHA256[1:]}, ValueError, "64 hexadecimal characters"),
        ({"sha256": SHA256 + "\n"}, ValueError, "64 hexadecimal characters"),
        ({"sha256": SHA256.encode()}, TypeError, "sha256 must be a str or None"),
    ],
)
def test_invalid_fields(kwargs, error, message):
    fields = {"name": "test", "build": 37, "tables": {"1": table([100])}} | kwargs
    with pytest.raises(error, match=message):
        GeneticMap(**fields)


def test_caller_frame_is_copied():
    frame = table([100, 200], cm=[0.0, 1.0])
    genetic_map = one_region(frame)
    frame.loc[0, "pos"] = 150
    frame["map"] = [5.0, 6.0]
    del frame["rate"]
    stored = genetic_map.tables["1"]
    assert stored["pos"].tolist() == [100, 200]
    assert stored["map"].tolist() == [0.0, 1.0]
    assert genetic_map.cm("1", 150).item() == 0.5


def test_caller_arrays_are_copied():
    pos = np.array([100, 200], dtype=np.uint32)
    cm = np.array([0.0, 1.0])
    frame = pd.DataFrame({"pos": pos, "rate": [1.0, 1.0], "map": cm}, copy=False)
    genetic_map = one_region(frame)
    pos[0] = 150
    cm[1] = 9.0
    assert genetic_map.tables["1"]["pos"].tolist() == [100, 200]
    assert genetic_map.tables["1"]["map"].tolist() == [0.0, 1.0]


def assert_unchanged(genetic_map):
    expected = pd.DataFrame(
        {
            "pos": np.array([100, 200, 400], dtype=np.uint32),
            "rate": [1.0, 1.0, 1.0],
            "map": [0.0, 1.0, 3.0],
        }
    )
    pd.testing.assert_frame_equal(genetic_map.tables["1"], expected)
    assert genetic_map.cm("1", [150, 300]).tolist() == [0.5, 2.0]


def test_column_and_loc_writes_to_a_returned_table(small_map):
    frame = small_map.tables["1"]
    frame["map"] = 9.0
    frame.loc[0, "pos"] = 50
    frame.iloc[1, 1] = 7.0
    del frame["rate"]
    assert frame["map"].tolist() == [9.0, 9.0, 9.0]
    assert frame["pos"].tolist() == [50, 200, 400]
    assert_unchanged(small_map)


def test_array_write_to_a_returned_table(small_map):
    frame = small_map.tables["1"]
    frame["map"].array[2] = 0.0
    assert frame["map"].tolist() == [0.0, 1.0, 0.0]
    assert_unchanged(small_map)


@pytest.mark.parametrize(
    "select",
    [
        lambda frame: frame[["rate", "map"]].to_numpy(),
        lambda frame: frame["map"].to_numpy(),
        lambda frame: frame["pos"].to_numpy(),
        lambda frame: frame.columns.to_numpy(),
        lambda frame: frame.index.to_numpy(),
    ],
    ids=["rate and map", "map", "pos", "columns", "index"],
)
def test_buffer_writes_to_a_returned_table(small_map, select):
    array = select(small_map.tables["1"])
    for buffer in (array, array.base):
        if buffer is not None and buffer.flags.writeable:
            buffer[...] = 0
    assert_unchanged(small_map)


def test_lookups_are_independent(small_map):
    first = small_map.tables["1"]
    second = small_map.tables["1"]
    assert first is not second
    first["map"].array[0] = 5.0
    assert second["map"].tolist() == [0.0, 1.0, 3.0]
    assert_unchanged(small_map)


@pytest.mark.parametrize(
    "copied",
    [copy.copy, copy.deepcopy, lambda value: pickle.loads(pickle.dumps(value))],
    ids=["copy", "deepcopy", "pickle"],
)
def test_copies(small_map, copied):
    result = copied(small_map)
    assert dataclasses.replace(result).regions == ("1",)
    assert (result.name, result.build, result.source, result.sha256) == (
        "test",
        37,
        "maps/test.txt",
        SHA256,
    )
    assert_unchanged(result)
    frame = result.tables["1"]
    frame["map"].array[0] = 5.0
    frame.loc[1, "pos"] = 150
    assert result.tables["1"] is not result.tables["1"]
    assert_unchanged(result)
    assert_unchanged(small_map)


def test_tables_compare_by_identity(small_map):
    tables = small_map.tables
    assert (tables == small_map.tables) is True
    assert (tables != small_map.tables) is False
    assert (tables == copy.deepcopy(small_map).tables) is False
    assert (tables == {"1": tables["1"]}) is False
    assert hash(tables) == hash(small_map.tables)


def test_tables_are_read_only(small_map):
    with pytest.raises(TypeError):
        small_map.tables["1"] = table([100])
    with pytest.raises(TypeError):
        small_map.tables["2"] = table([100])
    with pytest.raises(TypeError):
        del small_map.tables["1"]
    with pytest.raises(dataclasses.FrozenInstanceError):
        small_map.tables = {}
    with pytest.raises(dataclasses.FrozenInstanceError):
        small_map.name = "other"
    assert small_map.regions == ("1",)


def test_tables_mapping(small_map):
    tables = small_map.tables
    assert "1" in tables
    assert "2" not in tables
    assert ["1"] not in tables
    assert tables.get("2") is None
    assert list(tables.keys()) == ["1"]
    assert [len(frame) for frame in tables.values()] == [3]
    with pytest.raises(KeyError):
        tables["2"]
