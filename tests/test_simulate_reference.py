import numpy as np
import pandas as pd
import pytest

from lineage.maps import GeneticMap
from lineage.reconstruction import Pedigree, Person
from lineage.reconstruction.simulate import (
    Observation,
    Panel,
    _meiosis,
    simulate_pedigree,
    synthetic_panel,
)
from lineage.types import ABSENT, UNKNOWN, Loci, Sex


def _map_table(start, end, start_cm, end_cm):
    return pd.DataFrame(
        {
            "pos": [start, end],
            "rate": (end_cm - start_cm) * 1_000_000 / (end - start),
            "map": [start_cm, end_cm],
        }
    )


def test_meiosis_uses_full_map_span_and_switches_at_crossovers():
    class FixedDraws:
        def __init__(self):
            self.poisson_means = []
            self.crossovers = iter((67.0, 100.0, 50.0, 62.0))

        def poisson(self, lam=None, size=None, *, mean=None):
            self.poisson_means.append(lam if mean is None else mean)
            return 4 if size is None else np.full(size, 4)

        def uniform(self, low=None, high=None, size=None, *, start=None, end=None):
            assert (low if start is None else start) == 40.0
            assert (high if end is None else end) == 140.0
            if size is None:
                return next(self.crossovers)
            count = int(np.prod(size))
            return np.fromiter(self.crossovers, dtype=np.float64, count=count).reshape(
                size
            )

        def integers(self, low, high=None, size=None, dtype=np.int64, endpoint=False):
            bounds = (0, low) if high is None else (low, high)
            assert (bounds[0], bounds[1] + endpoint) == (0, 2)
            return 1 if size is None else np.ones(size, dtype=dtype)

    rng = FixedDraws()
    markers = np.array([43.0, 50.0, 59.0, 64.0, 68.0, 99.0, 100.0, 135.0])
    copies = _meiosis(rng, markers, start_cm=40.0, end_cm=140.0)

    # The 100 cM map span gives a Poisson mean of one even though the markers span 92 cM.
    assert rng.poisson_means == [1.0]
    # From copy 1, events at 50 and 100 count at the marker. The two events
    # between 59 and 68 cancel one another.
    np.testing.assert_array_equal(copies, [1, 0, 0, 1, 0, 0, 1, 1])
    assert copies.dtype == np.int8


def test_x_inheritance_and_haploid_observations():
    mother = Person("M1", Sex.FEMALE)
    father = Person("F1", Sex.MALE)
    daughter = Person("C1", Sex.FEMALE, mother=mother, father=father)
    son = Person("C2", Sex.MALE, mother=mother, father=father)
    pedigree = Pedigree([mother, father, daughter, son])
    map_x = GeneticMap.from_frames(
        "non-PAR X", {"X": _map_table(3_000_000, 3_000_100, 5.0, 105.0)}
    )
    pos = np.arange(3_000_000, 3_000_100, 10)
    alleles = np.array(
        [[2, 3, 1, 3, 3, 0, 0, 1, 1, 0], [0, 1, 2, 0, 1, 3, 2, 2, 3, 1]],
        dtype=np.int8,
    )
    panel_x = Panel(
        loci=Loci("X", pos),
        alleles=alleles,
        frequency=np.array([0, 1, 0, 1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]),
    )
    panel = {"X": panel_x}
    result = simulate_pedigree(
        pedigree,
        panel,
        map_x,
        tested={"C1", "C2"},
        observation={"C1": Observation(), "C2": Observation(error_rate=1.0)},
        seed=928,
    )

    mother_x = result.haplotypes["M1"]["X"]
    father_x = result.haplotypes["F1"]["X"]
    daughter_x = result.haplotypes["C1"]["X"]
    son_x = result.haplotypes["C2"]["X"]
    # At the first four sites the frequency fixes both founders' existing copies.
    np.testing.assert_array_equal(mother_x[:, :4], [[2, 1, 1, 0], [2, 1, 1, 0]])
    np.testing.assert_array_equal(father_x[0, :4], [2, 1, 1, 0])
    np.testing.assert_array_equal(father_x[1], ABSENT)
    for child, child_x in [("C1", daughter_x), ("C2", son_x)]:
        maternal = result.transmissions[("M1", child)]["X"]
        assert np.all(np.isin(maternal, [0, 1]))
        np.testing.assert_array_equal(
            child_x[0], mother_x[maternal, np.arange(pos.size)]
        )

    np.testing.assert_array_equal(result.transmissions[("F1", "C1")]["X"], 0)
    np.testing.assert_array_equal(daughter_x[1], father_x[0])
    np.testing.assert_array_equal(result.transmissions[("F1", "C2")]["X"], ABSENT)
    np.testing.assert_array_equal(son_x[1], ABSENT)
    np.testing.assert_array_equal(
        result.observations["C1"]["X"].calls, np.sort(daughter_x, axis=0)
    )

    # An error on the son's single called copy has exactly one possible replacement.
    replacement = np.where(son_x[0] == alleles[0], alleles[1], alleles[0])
    flipped = result.observations["C2"]["X"]
    assert flipped.assayed.dtype == np.bool_
    np.testing.assert_array_equal(flipped.assayed, np.ones(pos.size, dtype=bool))
    np.testing.assert_array_equal(flipped.calls[0], replacement)
    np.testing.assert_array_equal(flipped.calls[1], ABSENT)

    for observation, assayed in [
        (Observation(missing_rate=1.0), True),
        (Observation(chip_keep=0.0), False),
    ]:
        rerun = simulate_pedigree(
            pedigree, panel, map_x, tested={"C2"}, observation=observation, seed=928
        )
        calls = rerun.observations["C2"]["X"]
        # Missingness keeps the assay flag, while the chip mask clears it.
        np.testing.assert_array_equal(calls.assayed, np.full(pos.size, assayed))
        np.testing.assert_array_equal(calls.calls[0], UNKNOWN)
        np.testing.assert_array_equal(calls.calls[1], ABSENT)
        np.testing.assert_array_equal(rerun.haplotypes["C2"]["X"], son_x)


@pytest.mark.parametrize("founders", ["independent", "correlated"])
def test_retained_streams_survive_pedigree_and_test_selection_changes(founders):
    mother = Person("M1", Sex.FEMALE)
    father = Person("F1", Sex.MALE)
    daughter = Person("C2", Sex.FEMALE, mother="M1", father="F1")
    son = Person("C1", Sex.MALE, mother="M1", father="F1")
    spare = Person("C3", Sex.FEMALE)
    original = Pedigree([spare, mother, father, daughter])
    updated = Pedigree([mother, father, son, daughter])

    region_starts = {"2": 1_000, "X": 3_000_000}
    genetic_map = GeneticMap.from_frames(
        "stream fixture",
        {
            chrom: _map_table(start, start + 400, 10.0, 160.0)
            for chrom, start in region_starts.items()
        },
    )
    # Each region spans two template stretches in correlated mode.
    n = 260
    panels = {
        chrom: Panel(
            loci=Loci(chrom, np.arange(start, start + n)),
            alleles=np.repeat(np.array([[0], [3]], dtype=np.int8), n, axis=1),
            frequency=np.linspace(0.1, 0.6, n),
        )
        for chrom, start in region_starts.items()
    }
    observation = Observation(error_rate=0.35, missing_rate=0.2, chip_keep=0.7)
    original_tested = {"C2", "F1", "M1"}
    updated_tested = {"C1", "C2", "F1"}
    before = simulate_pedigree(
        original,
        panels,
        genetic_map,
        tested=original_tested,
        observation=observation,
        founders=founders,
        seed=531,
    )
    after = simulate_pedigree(
        updated,
        panels,
        genetic_map,
        tested=updated_tested,
        observation=observation,
        founders=founders,
        seed=531,
    )

    assert set(before.observations) == original_tested
    assert set(after.observations) == updated_tested
    for chrom in panels:
        for person in ("M1", "F1", "C2"):
            np.testing.assert_array_equal(
                before.haplotypes[person][chrom], after.haplotypes[person][chrom]
            )
        for edge in (("M1", "C2"), ("F1", "C2")):
            np.testing.assert_array_equal(
                before.transmissions[edge][chrom], after.transmissions[edge][chrom]
            )
        for person in original_tested & updated_tested:
            old_calls = before.observations[person][chrom]
            new_calls = after.observations[person][chrom]
            np.testing.assert_array_equal(old_calls.calls, new_calls.calls)
            np.testing.assert_array_equal(old_calls.assayed, new_calls.assayed)


@pytest.mark.parametrize(
    "count, expected", [(4, {"2": 2, "10": 1, "X": 1}), (2, {"2": 1, "10": 1})]
)
def test_synthetic_panel_allocates_total_in_natural_region_order(count, expected):
    genetic_map = GeneticMap.from_frames(
        "unordered",
        {
            "X": _map_table(3_000_010, 3_000_013, 0.0, 0.0),
            "10": _map_table(30, 33, 0.0, 0.0),
            "2": _map_table(20, 21, 0.0, 0.0),
        },
    )
    panels = synthetic_panel(genetic_map, count=count, seed=43)

    # Naturally 2 precedes 10 and X. Four gives 2 the extra site, and two skips X.
    assert tuple(panels) == tuple(expected)
    assert {chrom: panel.loci.pos.size for chrom, panel in panels.items()} == expected
    if count == 4:
        np.testing.assert_array_equal(panels["2"].loci.pos, [20, 21])
    for panel in panels.values():
        assert np.all(np.diff(panel.loci.pos) > 0)
        assert np.all(panel.alleles[0] != panel.alleles[1])
        assert np.all((panel.alleles >= 0) & (panel.alleles <= 3))
        assert np.all((panel.frequency >= 0.05) & (panel.frequency < 0.5))


@pytest.mark.parametrize(
    "start, end",
    [
        (60_000, 2_699_521),  # Immediately outside the two inclusive PAR1 boundaries.
        (
            154_931_043,
            155_260_561,
        ),  # Immediately outside the two inclusive PAR2 boundaries.
    ],
)
def test_synthetic_panel_excludes_par_across_sparse_x_map(start, end):
    genetic_map = GeneticMap.from_frames(
        "sparse X", {"X": _map_table(start, end, 3.0, 8.0)}
    )
    panels = synthetic_panel(genetic_map, count=2, chromosomes=["X"], seed=67)

    # Every bp in the interior is in PAR, so the endpoints are the only possible draw.
    assert set(panels) == {"X"}
    np.testing.assert_array_equal(panels["X"].loci.pos, [start, end])
    if start == 60_000:
        with pytest.raises(ValueError):
            synthetic_panel(genetic_map, count=3, chromosomes=["X"], seed=67)


def test_observation_error_and_missingness_preserve_homozygous_truth():
    pedigree = Pedigree([Person("C1", Sex.FEMALE)])
    genetic_map = GeneticMap.from_frames(
        "constant autosome", {"1": _map_table(100, 103, 0.0, 0.0)}
    )
    alleles = np.array([[3, 2, 0, 1], [0, 1, 2, 3]], dtype=np.int8)
    panel = {
        "1": Panel(
            loci=Loci("1", [100, 101, 102, 103]),
            alleles=alleles,
            frequency=np.array([0.0, 1.0, 0.0, 1.0]),
        )
    }
    errors = simulate_pedigree(
        pedigree, panel, genetic_map, observation=Observation(error_rate=1.0), seed=19
    )
    missing = simulate_pedigree(
        pedigree, panel, genetic_map, observation=Observation(missing_rate=1.0), seed=19
    )

    # Both copies must be the first, second, first and second panel allele respectively.
    expected_truth = np.array([[3, 1, 0, 3], [3, 1, 0, 3]], dtype=np.int8)
    np.testing.assert_array_equal(errors.haplotypes["C1"]["1"], expected_truth)
    np.testing.assert_array_equal(missing.haplotypes["C1"]["1"], expected_truth)
    corrupted = errors.observations["C1"]["1"]
    np.testing.assert_array_equal(
        corrupted.assayed, np.ones(alleles.shape[1], dtype=bool)
    )
    # Either other genotype contains the alternate allele at least once.
    assert corrupted.calls.shape == expected_truth.shape
    assert np.all(np.any(corrupted.calls != expected_truth, axis=0))
    assert np.all(corrupted.calls[0] <= corrupted.calls[1])
    assert np.all((corrupted.calls == alleles[0]) | (corrupted.calls == alleles[1]))

    no_calls = missing.observations["C1"]["1"]
    np.testing.assert_array_equal(no_calls.assayed, corrupted.assayed)
    np.testing.assert_array_equal(no_calls.calls, np.full_like(expected_truth, UNKNOWN))
