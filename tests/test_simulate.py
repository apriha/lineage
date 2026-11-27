"""Tests for synthetic inheritance, unphased observations and encoding."""

from dataclasses import FrozenInstanceError

import numpy as np
import pandas as pd
import pytest

from lineage.maps import GeneticMap
from lineage.reconstruction import simulate
from lineage.reconstruction.encoding import align, encode
from lineage.reconstruction.pedigree import Pedigree, Person
from lineage.reconstruction.simulate import (
    Observation,
    Panel,
    simulate_pedigree,
    synthetic_panel,
)
from lineage.types import ABSENT, UNKNOWN, Loci, Sex
from tests.support.pedigrees import collateral, connected, maternal, nuclear


@pytest.fixture
def map_and_panel():
    genetic_map = GeneticMap.uniform(chromosomes=("1", "X"))
    panel = synthetic_panel(genetic_map, count=1100, seed=13)
    return genetic_map, panel


def _check_edges(result):
    for person in result.pedigree:
        child = result.haplotypes[person.name]
        for row, parent in enumerate((person.mother, person.father)):
            if parent is None:
                continue
            for chrom, copies in result.transmissions[parent, person.name].items():
                if chrom == "X" and row == 1:
                    expected = 0 if person.sex == Sex.FEMALE else ABSENT
                    np.testing.assert_array_equal(copies, expected)
                present = copies != ABSENT
                markers = np.flatnonzero(present)
                truth = result.haplotypes[parent][chrom]
                np.testing.assert_array_equal(
                    child[chrom][row, present], truth[copies[present], markers]
                )
                np.testing.assert_array_equal(child[chrom][row, ~present], ABSENT)


@pytest.mark.parametrize(
    "fixture", [lambda: nuclear(2), collateral, connected, maternal]
)
def test_invented_pedigrees(fixture, map_and_panel):
    pedigree, tested = fixture()
    genetic_map, panel = map_and_panel
    result = simulate_pedigree(pedigree, panel, genetic_map, tested=tested, seed=17)
    assert tuple(result.haplotypes) == pedigree.names
    assert tuple(result.observations) == tested
    assert set(result.transmissions) == {
        (parent, person.name)
        for person in pedigree
        for parent in (person.mother, person.father)
        if parent is not None
    }
    _check_edges(result)
    for name in tested:
        for chrom, observed in result.observations[name].items():
            assert observed.assayed.all()
            truth = result.haplotypes[name][chrom]
            wanted = (
                truth
                if chrom == "X" and pedigree[name].sex == Sex.MALE
                else np.sort(truth, axis=0)
            )
            np.testing.assert_array_equal(observed.calls, wanted)


def test_sibship(map_and_panel):
    pedigree, _ = nuclear()
    genetic_map, panel = map_and_panel
    children = pedigree.children("M1")
    result = simulate_pedigree(pedigree, panel, genetic_map, tested=children, seed=1)
    assert tuple(result.observations) == children
    _check_edges(result)


@pytest.mark.parametrize("model", ["independent", "correlated"])
def test_multigeneration_templates_and_inheritance(monkeypatch, map_and_panel, model):
    genetic_map, panel = map_and_panel
    pedigree, _ = maternal()
    partial = Person("Child5", Sex.FEMALE, mother="Mother")
    pedigree = Pedigree([partial, *reversed(pedigree.people)])
    captured = []
    draw_templates = simulate._templates

    def capture(rng, frequency):
        templates = draw_templates(rng, frequency)
        captured.append(templates)
        return templates

    monkeypatch.setattr(simulate, "_templates", capture)
    result = simulate_pedigree(pedigree, panel, genetic_map, founders=model, seed=29)
    _check_edges(result)
    if model == "independent":
        assert not captured
        return

    for (chrom, block), templates in zip(panel.items(), captured, strict=True):
        assert templates.shape == (12, 550)
        for person in pedigree:
            for copy, parent in enumerate((person.mother, person.father)):
                if parent is not None or (
                    chrom == "X" and copy == 1 and person.sex == Sex.MALE
                ):
                    continue
                allele = result.haplotypes[person.name][chrom][copy]
                choices = allele == block.alleles[1]
                for offset in (0, 256, 512):
                    stretch = slice(offset, offset + 256)
                    assert np.any(
                        np.all(templates[:, stretch] == choices[stretch], axis=1)
                    )


def test_encode_and_align_only_observed_assays(map_and_panel):
    pedigree, _ = nuclear()
    genetic_map, panel = map_and_panel
    settings = {
        "C1": Observation(chip_keep=0.25, missing_rate=1),
        "C2": Observation(chip_keep=0.25),
    }
    result = simulate_pedigree(
        pedigree, panel, genetic_map, tested=["C2", "C1"], observation=settings, seed=9
    )
    encoded = {}
    sexes = {name: pedigree[name].sex for name in result.observations}
    for name, blocks in result.observations.items():
        rows = []
        for chrom, observed in blocks.items():
            for pos, call in zip(
                panel[chrom].loci.pos[observed.assayed],
                observed.calls[:, observed.assayed].T,
                strict=True,
            ):
                genotype = (
                    None
                    if UNKNOWN in call
                    else "".join("ACGT"[code] for code in call if code >= 0)
                )
                rows.append((f"{chrom}:{pos}", chrom, pos, genotype))
        frame = pd.DataFrame(rows, columns=["rsid", "chrom", "pos", "genotype"])
        encoded[name], _ = encode(frame.set_index("rsid"), sex=sexes[name])

    aligned = align(encoded, sexes=sexes)
    assert tuple(aligned) == ("C1", "C2")
    for chrom, block in panel.items():
        c1 = result.observations["C1"][chrom].assayed
        c2 = result.observations["C2"][chrom].assayed
        union = c1 | c2
        assert np.any(c1 & ~c2) and np.any(~union)
        for name in aligned:
            calls = aligned[name][chrom]
            np.testing.assert_array_equal(calls.loci.pos, block.loci.pos[union])
            np.testing.assert_array_equal(calls.phase_set, -1)
            np.testing.assert_array_equal(calls.kinds[:, 1:], 0)
        np.testing.assert_array_equal(aligned["C1"][chrom].kinds[:, 0], 0)
        np.testing.assert_array_equal(aligned["C2"][chrom].kinds[:, 0], c2[union])


@pytest.mark.parametrize("bad", [-0.01, 1.01, float("nan"), float("inf"), 2**100])
def test_probability_range(bad):
    for name in ("error_rate", "missing_rate", "chip_keep"):
        with pytest.raises(ValueError, match=name):
            Observation(**{name: bad})
    with pytest.raises(ValueError, match="frequency"):
        Panel(Loci("1", [1]), [[0], [1]], np.array([bad], dtype=np.float64))


def test_frequency_uses_validated_array_values():
    class Frequency:
        def __array__(self, dtype=None, copy=None):
            value = (
                np.nan if dtype is not None and dtype == np.dtype(np.float64) else 0.5
            )
            return np.array([value], dtype=dtype)

    panel = Panel(Loci("1", [1]), [[0], [1]], Frequency())
    np.testing.assert_array_equal(panel.frequency, [0.5])


@pytest.mark.parametrize("bad", [True, np.bool_(False), "0.5", [0.5], None])
def test_probability_type(bad):
    with pytest.raises(TypeError, match="error_rate"):
        Observation(error_rate=bad)


@pytest.mark.parametrize(
    "bad",
    [
        [[0]],
        [[0], [0]],
        [[0], [257]],
        [[0], [2**100]],
        [[0], [1.5]],
        [[0], [np.nan]],
        [[True], [0]],
        [["1"], ["2"]],
    ],
)
def test_allele_validation(bad):
    with pytest.raises(ValueError, match="alleles"):
        Panel(Loci("1", [1]), bad, [0.5])


@pytest.mark.parametrize("bad", [0.5, np.nan])
def test_invalid_allele_in_last_column(bad):
    with pytest.raises(ValueError, match="alleles"):
        Panel(Loci("1", [1, 2]), [[0, 0], [1, bad]], [0.5, 0.5])


@pytest.mark.parametrize(
    "value, error",
    [(0, ValueError), (-1, ValueError), (True, TypeError), (2.0, TypeError)],
)
def test_count(value, error, map_and_panel):
    genetic_map, _ = map_and_panel
    with pytest.raises(error, match="count"):
        synthetic_panel(genetic_map, count=value)


@pytest.mark.parametrize(
    "seed, error",
    [(-1, ValueError), (False, TypeError), (3.0, TypeError), ("3", TypeError)],
)
def test_seed(seed, error, map_and_panel):
    genetic_map, panel = map_and_panel
    with pytest.raises(error, match="seed"):
        synthetic_panel(genetic_map, count=2, seed=seed)
    with pytest.raises(error, match="seed"):
        simulate_pedigree(Pedigree([]), panel, genetic_map, seed=seed)


def test_region_and_people_validation(map_and_panel):
    genetic_map, panel = map_and_panel
    pedigree, _ = nuclear()
    for regions in ([], ["1", "1"], ["2"], ["Y"], ["X_par1"]):
        with pytest.raises(ValueError):
            synthetic_panel(genetic_map, count=3, chromosomes=regions)
    with pytest.raises(TypeError):
        synthetic_panel(genetic_map, count=3, chromosomes="1")
    for blocks in ({"X": panel["1"]}, {"2": panel["1"]}, {"X_par1": panel["X"]}):
        with pytest.raises(ValueError):
            simulate_pedigree(pedigree, blocks, genetic_map)
    with pytest.raises(TypeError):
        simulate_pedigree(pedigree, {"1": None}, genetic_map)
    for tested in (["Nobody"], ["C1", "C1"]):
        with pytest.raises(ValueError, match=r"C1.*C4"):
            simulate_pedigree(pedigree, panel, genetic_map, tested=tested)
    with pytest.raises(ValueError, match="Nobody"):
        simulate_pedigree(
            pedigree, panel, genetic_map, observation={"Nobody": Observation()}
        )
    with pytest.raises(TypeError, match="C1"):
        simulate_pedigree(pedigree, panel, genetic_map, observation={"C1": None})
    with pytest.raises(ValueError, match="founders"):
        simulate_pedigree(pedigree, panel, genetic_map, founders="unsupported")


def test_frozen_results(map_and_panel):
    genetic_map, panel = map_and_panel
    pedigree, _ = nuclear()
    supplied = dict(panel)
    result = simulate_pedigree(pedigree, supplied, genetic_map, seed=7)
    assert result.pedigree is pedigree
    assert result.genetic_map is genetic_map
    supplied.clear()
    assert tuple(result.panel) == ("1", "X")
    alleles = np.array([[0], [1]], dtype=np.int8)
    frequency = np.array([0.5])
    block = Panel(Loci("1", [1]), alleles, frequency)
    alleles[:] = 3
    frequency[0] = 0
    np.testing.assert_array_equal(block.alleles, [[0], [1]])
    np.testing.assert_array_equal(block.frequency, [0.5])
    with pytest.raises(FrozenInstanceError):
        result.observations["C1"]["1"].calls = alleles
    for buffer in (
        block.alleles,
        block.frequency,
        result.haplotypes["C1"]["1"],
        result.transmissions["M1", "C1"]["1"],
        result.observations["C1"]["1"].calls,
        result.observations["C1"]["1"].assayed,
    ):
        while isinstance(buffer, np.ndarray):
            assert not buffer.flags.writeable
            with pytest.raises(ValueError):
                buffer.setflags(write=True)
            buffer = buffer.base
        assert isinstance(buffer, bytes)
    for mapping in (
        panel,
        result.panel,
        result.haplotypes,
        result.haplotypes["C1"],
        result.transmissions,
        result.transmissions["M1", "C1"],
        result.observations,
        result.observations["C1"],
        result.parameters,
        result.parameters["observation"],
    ):
        with pytest.raises(TypeError):
            mapping["new"] = None


@pytest.mark.parametrize("model", ["independent", "correlated"])
def test_replay_with_reordered_people_and_panel(map_and_panel, model):
    genetic_map, panel = map_and_panel
    pedigree, tested = nuclear()
    settings = {"C1": Observation(error_rate=0.4, missing_rate=0.2, chip_keep=0.5)}
    original = simulate_pedigree(
        pedigree,
        panel,
        genetic_map,
        tested=tested,
        observation=settings,
        founders=model,
    )
    seed = original.parameters["seed"]
    assert isinstance(seed, int)
    reordered = simulate_pedigree(
        Pedigree(reversed(pedigree.people)),
        dict(reversed(list(panel.items()))),
        genetic_map,
        tested=tested,
        observation=settings,
        founders=model,
        seed=seed,
    )
    np.testing.assert_array_equal(
        original.haplotypes["C1"]["1"], reordered.haplotypes["C1"]["1"]
    )
    np.testing.assert_array_equal(
        original.transmissions["F1", "C1"]["1"],
        reordered.transmissions["F1", "C1"]["1"],
    )
    before, after = original.observations["C1"]["1"], reordered.observations["C1"]["1"]
    np.testing.assert_array_equal(before.calls, after.calls)
    np.testing.assert_array_equal(before.assayed, after.assayed)
    assert tuple(reordered.panel) == ("1", "X")
    assert tuple(reordered.haplotypes) == tuple(reversed(pedigree.names))
    assert tuple(reordered.observations) == tuple(reversed(tested))
    assert original.parameters["tested"] == tested
    assert original.parameters["observation"]["F1"] == Observation()
    assert original.parameters["model_version"] == "1"
    assert original.parameters["map_name"] == genetic_map.name
    assert original.parameters["map_sha256"] == genetic_map.sha256


@pytest.mark.parametrize("model", ["independent", "correlated"])
def test_unlisted_parent_has_a_private_founding_copy(map_and_panel, model):
    genetic_map, panel = map_and_panel
    pedigree, _ = nuclear()
    original = simulate_pedigree(pedigree, panel, genetic_map, founders=model, seed=59)
    changed_people = [
        Person("C1", Sex.FEMALE, "M1") if p.name == "C1" else p for p in pedigree
    ]
    changed = simulate_pedigree(
        Pedigree(changed_people),
        panel,
        genetic_map,
        tested=[],
        founders=model,
        seed=59,
    )
    assert ("F1", "C1") not in changed.transmissions
    np.testing.assert_array_equal(
        changed.haplotypes["C1"]["1"][0], original.haplotypes["C1"]["1"][0]
    )
    assert np.any(changed.haplotypes["C1"]["1"][1] != original.haplotypes["C1"]["1"][1])


def test_empty_axis_and_one_row_map():
    genetic_map = GeneticMap.from_frames(
        "one site", {"1": pd.DataFrame({"pos": [10], "map": [10.0], "rate": [0.0]})}
    )
    one = synthetic_panel(genetic_map, count=1, seed=1)
    np.testing.assert_array_equal(one["1"].loci.pos, [10])
    with pytest.raises(ValueError, match="room for 1"):
        synthetic_panel(genetic_map, count=2, seed=1)
    pedigree, _ = nuclear()
    assert not simulate_pedigree(Pedigree([]), one, genetic_map).haplotypes
    empty = simulate_pedigree(pedigree, {}, genetic_map, tested=[])
    assert all(not blocks for blocks in empty.haplotypes.values())
    assert all(not blocks for blocks in empty.transmissions.values())
    assert not empty.observations

    clamped = Panel(Loci("1", [1, 20]), [[0, 1], [2, 3]], [0.5, 0.5])
    result = simulate_pedigree(pedigree, {"1": clamped}, genetic_map, seed=2)
    for edge in result.transmissions.values():
        copies = edge["1"]
        assert copies[0] == copies[1]
    _check_edges(result)


def test_empty_axis_with_large_map_span():
    genetic_map = GeneticMap.from_frames(
        "wide span",
        {"1": pd.DataFrame({"pos": [1, 2], "map": [0.0, 1e30], "rate": [0.0, 0.0]})},
    )
    pedigree, _ = nuclear()
    block = Panel(Loci("1", []), np.empty((2, 0), dtype=np.int8), [])
    result = simulate_pedigree(
        pedigree, {"1": block}, genetic_map, founders="correlated"
    )
    assert result.observations["C1"]["1"].calls.shape == (2, 0)
    copies = result.transmissions["M1", "C1"]["1"]
    assert copies.shape == (0,)
    assert copies.dtype == np.int8
