"""Integration tests for metric inputs and pooled results."""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pandas as pd
import pytest

from lineage.maps import GeneticMap
from lineage.reconstruction.encoding import _Calls, _pair
from lineage.reconstruction.metrics import (
    _CoupleCalls,
    _EvaluationInput,
    evaluate_reconstruction,
)
from lineage.reconstruction.pedigree import Pedigree, Person
from lineage.reconstruction.simulate import Panel, simulate_pedigree
from lineage.types import Loci, Sex


def _panel(chrom: str, positions: list[int]) -> Panel:
    n = len(positions)
    return Panel(Loci(chrom, positions), [[0] * n, [2] * n], [0.0] * n)


def _calls(chrom: str, positions: list[int], labels: list[int], allele: int) -> _Calls:
    return _Calls(
        Loci(chrom, positions, ["other"] * len(positions)),
        [_pair(allele, allele)] * len(positions),
        labels,
        [[2, 2 if label >= 0 else 0, 0] for label in labels],
    )


def test_pooled_n50_alignment_and_copy_on_write() -> None:
    pedigree = Pedigree([Person("C1", Sex.FEMALE)])
    short = [1, 100_001, 200_001, 900_001]
    long = [1, 1_000_001]
    panel = {"10": _panel("10", long), "2": _panel("2", short)}
    truth = simulate_pedigree(
        pedigree, panel, GeneticMap.uniform(chromosomes=tuple(panel)), seed=1
    )
    named = {
        "2": _calls("2", [*short[:3], 400_001], [0, 0, 0, 0], 0),
        "10": _calls("10", long, [0, 0], 0),
    }
    source = {"C1": named}
    inputs = _EvaluationInput(pedigree, {}, source)
    named.clear()
    source.clear()

    result = evaluate_reconstruction(inputs, truth)
    table = result.table
    assert list(table["chrom"]) == ["2", "10"]
    assert list(table["phase_n50_cm"]) == pytest.approx([0.2, 1.0])
    assert list(table["reported_copies"]) == [6, 4]
    summary = result.summary()
    assert summary.at[0, "phase_sets"] == 2
    assert summary.at[0, "phase_n50_cm"] == pytest.approx(0.2)
    assert summary.at[0, "coverage"] == pytest.approx(10 / 12)

    table.loc[0, "known_copies"] = 0
    summary.loc[0, "coverage"] = 0
    assert result.table.loc[0, "known_copies"] == 6
    assert result.summary().loc[0, "coverage"] == pytest.approx(10 / 12)
    assert isinstance(table.index, pd.RangeIndex)


def test_couple_metadata_and_off_panel_identity_anchor() -> None:
    mother = Person("M1", Sex.FEMALE)
    father = Person("F1", Sex.MALE)
    pedigree = Pedigree([father, mother])
    panel = {"1": _panel("1", [1, 2, 3])}
    truth = simulate_pedigree(
        pedigree, panel, GeneticMap.uniform(chromosomes=("1",)), seed=1
    )
    labels = np.array([7, 7, 7], dtype=np.int32)
    roles = np.array([-1, -1, 0], dtype=np.int8)
    a = _calls("1", [1, 2, 4], [-1, -1, -1], 0)
    b = _calls("1", [1, 2, 4], [-1, -1, -1], 2)
    couple = _CoupleCalls("M1", "F1", a, b, labels, roles)
    labels[:] = -1
    roles[:] = 1
    assert list(couple.label_set) == [7, 7, 7]
    assert list(couple.a_is) == [-1, -1, 0]
    with pytest.raises(ValueError, match="WRITEABLE"):
        couple.a_is.setflags(write=True)
    with pytest.raises(ValueError, match="share an autosomal chromosome and positions"):
        _CoupleCalls("M1", "F1", a, _calls("1", [1, 2, 5], [-1] * 3, 2), labels, roles)
    with pytest.raises(ValueError, match="both parents"):
        _CoupleCalls(
            "M1", "F1", a, b, couple.label_set, np.array([0, 1, -1], dtype=np.int8)
        )

    evaluated = evaluate_reconstruction(
        _EvaluationInput(pedigree, {}, {}, (couple,)), truth
    )
    table = evaluated.table.set_index("person")
    assert table.at["M1+F1", "eligible_copies"] == 12
    assert table.at["M1+F1", "reported_copies"] == 8
    assert table.at["M1", "reported_copies"] == table.at["F1", "reported_copies"] == 0


def test_anonymous_tied_recovery_is_invariant_to_labels_and_bases() -> None:
    pedigree = Pedigree([Person("M1", Sex.FEMALE), Person("F1", Sex.MALE)])
    positions = [1, 2, 3, 4, 5]
    genetic_map = GeneticMap.uniform(chromosomes=("1",))
    seed_truth = simulate_pedigree(
        pedigree, {"1": _panel("1", positions)}, genetic_map, seed=1
    )
    label_set = np.array([0, 0, 1, 1, -1], dtype=np.int32)
    a_is = np.full(5, -1, dtype=np.int8)
    loci = Loci("1", positions)

    def table(perm: tuple[int, int, int, int], swap_first: bool) -> pd.DataFrame:
        first, other, last, _ = perm
        truth = replace(
            seed_truth,
            panel={
                "1": Panel(
                    loci,
                    [
                        [first, first, first, first, first],
                        [last, last, last, last, last],
                    ],
                    [0.5] * 5,
                )
            },
            haplotypes={
                "M1": {"1": np.array([[first] * 5, [last] * 5], dtype=np.int8)},
                "F1": {"1": np.full((2, 5), first, dtype=np.int8)},
            },
        )
        forced, estimated = (_pair(last, last), 2), (_pair(other, last), 3)
        a = (
            [estimated, estimated, forced, forced, forced]
            if swap_first
            else [forced] * 5
        )
        b = (
            [forced, forced, estimated, estimated, estimated]
            if swap_first
            else [estimated] * 5
        )

        def block(calls: list[tuple[int, int]]) -> _Calls:
            return _Calls(
                loci,
                [domain for domain, _ in calls],
                [-1] * 5,
                [[kind, 0, 0] for _, kind in calls],
            )

        couple = _CoupleCalls("M1", "F1", block(a), block(b), label_set, a_is)
        return evaluate_reconstruction(
            _EvaluationInput(pedigree, {}, {}, (couple,)), truth
        ).table

    original = table((0, 1, 2, 3), False)
    row = original.set_index("person").loc["M1+F1"]
    assert row["correct_copies"] == row["forced_correct"] == 5
    assert row["complete_correct"] == row["estimated_correct"] == 0
    pd.testing.assert_frame_equal(original, table((0, 1, 2, 3), True))
    pd.testing.assert_frame_equal(original, table((2, 3, 0, 1), False))


@pytest.mark.parametrize("labels", [[0, 1], [-1, -1]])
def test_anonymous_phase_set_requires_linked_identity(labels: list[int]) -> None:
    phased = _calls("1", [1, 2], [0, 0], 0)
    unphased = _calls("1", [1, 2], [-1, -1], 2)
    with pytest.raises(ValueError, match="multi-locus anonymous phase set"):
        _CoupleCalls(
            "M1",
            "F1",
            phased,
            unphased,
            np.asarray(labels, dtype=np.int32),
            np.full(2, -1, dtype=np.int8),
        )


def test_duplicate_and_colliding_couple_rows() -> None:
    mother = Person("M1", Sex.FEMALE)
    father = Person("F1", Sex.MALE)
    a = _calls("1", [1], [-1], 0)
    b = _calls("1", [1], [-1], 2)
    couple = _CoupleCalls(
        "M1", "F1", a, b, np.array([-1], dtype=np.int32), np.array([-1], dtype=np.int8)
    )
    pedigree = Pedigree([mother, father])
    panel = {"1": _panel("1", [1])}
    genetic_map = GeneticMap.uniform(chromosomes=("1",))
    truth = simulate_pedigree(pedigree, panel, genetic_map, seed=1)
    with pytest.raises(ValueError, match="Multiple blocks for couple"):
        evaluate_reconstruction(
            _EvaluationInput(pedigree, {}, {}, (couple, couple)), truth
        )

    pedigree = Pedigree([mother, father, Person("M1+F1", Sex.FEMALE)])
    truth = simulate_pedigree(pedigree, panel, genetic_map, seed=1)
    with pytest.raises(ValueError, match="Anonymous couple name.*is ambiguous"):
        evaluate_reconstruction(_EvaluationInput(pedigree, {}, {}, (couple,)), truth)


def test_anonymous_linked_phase_breaks_equal_genotype_scores() -> None:
    pedigree = Pedigree([Person("M1", Sex.FEMALE), Person("F1", Sex.MALE)])
    panel = {"1": _panel("1", [1, 2, 3])}
    seed_truth = simulate_pedigree(
        pedigree, panel, GeneticMap.uniform(chromosomes=("1",)), seed=1
    )
    truth = replace(
        seed_truth,
        haplotypes={
            "M1": {"1": np.asarray([[0, 0, 0], [2, 2, 2]], dtype=np.int8)},
            "F1": {"1": np.asarray([[2, 2, 0], [0, 0, 2]], dtype=np.int8)},
        },
    )
    loci = panel["1"].loci
    phased = _Calls(loci, [_pair(0, 2)] * 3, [0, 0, 0], [[3, 3, 0]] * 3)
    unordered = _Calls(
        loci, [_pair(0, 2) | _pair(2, 0)] * 3, [-1, -1, -1], [[2, 0, 0]] * 3
    )
    label_set = np.zeros(3, dtype=np.int32)
    a_is = np.full(3, -1, dtype=np.int8)

    def evaluate(a: _Calls, b: _Calls) -> pd.DataFrame:
        couple = _CoupleCalls("M1", "F1", a, b, label_set, a_is)
        return evaluate_reconstruction(
            _EvaluationInput(pedigree, {}, {}, (couple,)), truth
        ).table

    actual = evaluate(phased, unordered)
    row = actual.set_index("person").loc["M1+F1"]
    assert row["complete_correct"] == 6
    assert row["phase_correct"] == 6
    assert row["switch_opportunities"] == 2
    assert row["switch_errors"] == 0
    pd.testing.assert_frame_equal(actual, evaluate(unordered, phased))


def test_unordered_heterozygotes_do_not_contribute_to_phase_coverage() -> None:
    pedigree = Pedigree([Person("C1", Sex.FEMALE)])
    panel = {"1": _panel("1", [1, 2])}
    seed_truth = simulate_pedigree(
        pedigree, panel, GeneticMap.uniform(chromosomes=("1",)), seed=1
    )
    truth = replace(
        seed_truth,
        haplotypes={"C1": {"1": np.asarray([[0, 0], [2, 2]], dtype=np.int8)}},
    )
    calls = _Calls(
        panel["1"].loci, [_pair(0, 2) | _pair(2, 0)] * 2, [0, 0], [[2, 2, 0]] * 2
    )
    evaluated = evaluate_reconstruction(
        _EvaluationInput(pedigree, {}, {"C1": {"1": calls}}), truth
    ).table.iloc[0]
    assert evaluated["phase_eligible_heterozygotes"] == 2
    assert evaluated["complete_correct"] == 2
    assert evaluated["phased_copies"] == evaluated["phased_heterozygotes"] == 0
    assert evaluated["phase_coverage"] == 0
