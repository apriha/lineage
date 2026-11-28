from __future__ import annotations

import numpy as np
import pytest

from lineage.maps import GeneticMap
from lineage.reconstruction.encoding import _Calls
from lineage.reconstruction.metrics import (
    _CoupleCalls,
    _EvaluationInput,
    evaluate_reconstruction,
)
from lineage.reconstruction.pedigree import Pedigree, Person
from lineage.reconstruction.simulate import Panel, _SimulatedArrays
from lineage.types import Loci, Sex

A, G = 0, 2
ABSENT = -2
ABSENT_SLOT = 4
MISSING, OBSERVED, FORCED, ESTIMATED = 0, 1, 2, 3


def _domain(*pairs):
    return sum(1 << (5 * first + second) for first, second in set(pairs))


DIPLOID_MISSING = _domain(*((a, b) for a in range(4) for b in range(4)))
HAPLOID_MISSING = _domain(*((a, ABSENT_SLOT) for a in range(4)))
AA = _domain((A, A))
AG = _domain((A, G))
GA = _domain((G, A))
GG = _domain((G, G))
UNORDERED_AG = AG | GA
UNORDERED_A = _domain(*[(A, b) for b in range(4)], *[(a, A) for a in range(4)])


def _loci(chrom, positions, rsid=None):
    if rsid is not None:
        rsid = np.asarray(rsid)
    return Loci(chrom, np.asarray(positions, dtype=np.int64), rsid=rsid)


def _calls(loci, domains, genotype=FORCED, phase_set=None, origin=MISSING):
    allowed = np.asarray(domains, dtype=np.uint32)
    if phase_set is None:
        phase_set = np.full(len(allowed), -1, dtype=np.int32)
    else:
        phase_set = np.asarray(phase_set, dtype=np.int32)
    kinds = np.zeros((len(allowed), 3), dtype=np.uint8)
    kinds[:, 0] = genotype
    kinds[np.isin(allowed, [DIPLOID_MISSING, HAPLOID_MISSING]), 0] = MISSING
    kinds[phase_set >= 0, 1] = ESTIMATED
    kinds[:, 2] = origin
    return _Calls(loci, allowed, phase_set, kinds)


def _truth(pedigree, axes, haplotypes):
    panel = {}
    for chrom, loci in axes.items():
        n = next(len(chromosomes[chrom][0]) for chromosomes in haplotypes.values())
        panel[chrom] = Panel(
            loci,
            np.asarray([[A] * n, [G] * n], dtype=np.int8),
            np.full(n, 0.5),
        )
    arrays = {
        person: {
            chrom: np.asarray(copies, dtype=np.int8)
            for chrom, copies in chromosomes.items()
        }
        for person, chromosomes in haplotypes.items()
    }
    genetic_map = GeneticMap.uniform(rate=1.0, chromosomes=("1", "X"))
    return _SimulatedArrays(pedigree, panel, genetic_map, arrays, {}, {}, {})


def _row(evaluation, person, chrom):
    table = evaluation.table
    matching = table[(table["person"] == person) & (table["chrom"] == chrom)]
    assert len(matching) == 1
    return matching.iloc[0]


def _check(row, **expected):
    for column, value in expected.items():
        assert row[column] == pytest.approx(value, nan_ok=True), column


@pytest.mark.parametrize(
    ("domains", "origin_correct", "phase_correct", "switch_errors"),
    [
        pytest.param([AG, AG, AG], 3, 6, 0, id="perfect"),
        pytest.param([GA, GA, GA], 0, 6, 0, id="whole-set-flip"),
        pytest.param([AG, GA, GA], 1, 4, 1, id="one-switch"),
    ],
)
def test_phase_orientation_and_parental_origin(
    domains, origin_correct, phase_correct, switch_errors
):
    loci = _loci("1", [1_000_000, 2_000_000, 4_000_000])
    pedigree = Pedigree([Person("C1", Sex.FEMALE)])
    truth = _truth(pedigree, {"1": loci}, {"C1": {"1": [[A] * 3, [G] * 3]}})
    calls = _calls(loci, domains, phase_set=[7, 7, 7], origin=FORCED)
    result = _EvaluationInput(pedigree, {}, {"C1": {"1": calls}})

    # All three genotypes match. Phase takes one orientation across the set,
    # while origin uses the reported order. The span is 3 Mb at 1 cM per Mb.
    _check(
        _row(evaluate_reconstruction(result, truth), "C1", "1"),
        eligible_copies=6,
        observed_copies=0,
        reported_copies=6,
        known_copies=6,
        correct_copies=6,
        coverage=1,
        consistency=1,
        forced_known=6,
        forced_correct=6,
        estimated_known=0,
        estimated_correct=0,
        complete_known=3,
        complete_correct=3,
        complete_accuracy=1,
        phased_copies=6,
        phase_correct=phase_correct,
        phase_sets=1,
        phase_n50_cm=3,
        phase_eligible_heterozygotes=3,
        phased_heterozygotes=3,
        phase_coverage=1,
        switch_opportunities=2,
        switch_errors=switch_errors,
        origin_known=3,
        origin_correct=origin_correct,
    )


def test_wrong_genotype_and_no_call_keep_the_truth_denominator():
    loci = _loci("1", [1_000_000, 2_000_000, 3_000_000])
    pedigree = Pedigree([Person("C1", Sex.FEMALE)])
    truth = _truth(pedigree, {"1": loci}, {"C1": {"1": [[A] * 3, [G] * 3]}})
    calls = _calls(loci, [AA, DIPLOID_MISSING, UNORDERED_AG])
    result = _EvaluationInput(pedigree, {}, {"C1": {"1": calls}})

    # AA contributes two known copies and one correct copy against AG.
    # The no-call contributes two eligible copies and no reported copies.
    _check(
        _row(evaluate_reconstruction(result, truth), "C1", "1"),
        eligible_copies=6,
        observed_copies=0,
        reported_copies=4,
        known_copies=4,
        correct_copies=3,
        coverage=4 / 6,
        consistency=3 / 4,
        complete_known=2,
        complete_correct=1,
        complete_accuracy=1 / 2,
        phase_eligible_heterozygotes=3,
        phase_coverage=0,
    )


def test_unordered_partial_domains_contribute_minimum_dosage():
    loci = _loci("1", [1_000_000, 2_000_000])
    pedigree = Pedigree([Person("C1", Sex.FEMALE)])
    truth = _truth(pedigree, {"1": loci}, {"C1": {"1": [[A, A], [G, G]]}})
    observed = _calls(loci, [DIPLOID_MISSING, UNORDERED_A], genotype=OBSERVED)
    calls = _calls(loci, [UNORDERED_A, UNORDERED_AG], genotype=[FORCED, ESTIMATED])
    result = _EvaluationInput(pedigree, {"C1": {"1": observed}}, {"C1": {"1": calls}})

    # The first domain guarantees one A across all orders. The second site
    # already had one A, so completing AG adds just the estimated G.
    _check(
        _row(evaluate_reconstruction(result, truth), "C1", "1"),
        eligible_copies=4,
        observed_copies=1,
        reported_copies=3,
        known_copies=2,
        correct_copies=2,
        coverage=2 / 3,
        consistency=1,
        forced_known=1,
        forced_correct=1,
        estimated_known=1,
        estimated_correct=1,
        complete_known=1,
        complete_correct=1,
        complete_accuracy=1,
    )


def test_parental_origin_counts_partial_diploid_calls_by_locus():
    loci = _loci("1", [1_000_000, 2_000_000])
    pedigree = Pedigree([Person("C1", Sex.FEMALE)])
    truth = _truth(pedigree, {"1": loci}, {"C1": {"1": [[A, A], [G, G]]}})
    maternal_a = _domain(*((A, b) for b in range(4)))
    paternal_a = _domain(*((a, A) for a in range(4)))
    calls = _calls(loci, [maternal_a, paternal_a], phase_set=[4, 4], origin=FORCED)
    result = _EvaluationInput(pedigree, {}, {"C1": {"1": calls}})

    # Both reported As match the truth dosage. Only the first is assigned
    # to the maternal copy, so one of the two known origin loci is correct.
    _check(
        _row(evaluate_reconstruction(result, truth), "C1", "1"),
        known_copies=2,
        correct_copies=2,
        origin_known=2,
        origin_correct=1,
    )


@pytest.mark.parametrize(
    ("observed_base", "origin_correct"),
    [
        pytest.param(G, 2, id="correct-observation"),
        pytest.param(A, 1, id="wrong-observation"),
    ],
)
def test_male_x_no_call_and_physically_absent_copy(observed_base, origin_correct):
    loci = _loci("X", [10_000_000, 11_000_000, 12_000_000])
    pedigree = Pedigree([Person("C1", Sex.MALE)])
    truth = _truth(pedigree, {"X": loci}, {"C1": {"X": [[A, G, A], [ABSENT] * 3]}})
    a = _domain((A, ABSENT_SLOT))
    observed_allele = _domain((observed_base, ABSENT_SLOT))
    observed = _calls(
        loci, [HAPLOID_MISSING, observed_allele, HAPLOID_MISSING], genotype=OBSERVED
    )
    calls = _calls(
        loci,
        [a, observed_allele, HAPLOID_MISSING],
        genotype=[FORCED, OBSERVED, MISSING],
        origin=[FORCED, MISSING, MISSING],
    )
    result = _EvaluationInput(pedigree, {"C1": {"X": observed}}, {"C1": {"X": calls}})

    # Three non-PAR sites have one physical copy each. One is observed, one
    # newly known, and one unknown. The absent paternal copies never enter.
    # Both reported alleles have maternal origin, including the raw observation.
    # An incorrect observation reduces origin accuracy without changing fill.
    _check(
        _row(evaluate_reconstruction(result, truth), "C1", "X"),
        eligible_copies=3,
        observed_copies=1,
        reported_copies=2,
        known_copies=1,
        correct_copies=1,
        coverage=1 / 2,
        consistency=1,
        origin_known=2,
        origin_correct=origin_correct,
        phase_eligible_heterozygotes=0,
        phase_coverage=np.nan,
        phase_n50_cm=np.nan,
    )


def test_new_phase_on_fully_observed_genotypes_adds_no_genotype_copies():
    loci = _loci("1", [1_000_000, 2_000_000])
    pedigree = Pedigree([Person("C1", Sex.FEMALE)])
    truth = _truth(pedigree, {"1": loci}, {"C1": {"1": [[A, A], [G, G]]}})
    observed = _calls(loci, [UNORDERED_AG] * 2, genotype=OBSERVED)
    calls = _calls(loci, [AG] * 2, genotype=OBSERVED, phase_set=[3, 3])
    result = _EvaluationInput(pedigree, {"C1": {"1": observed}}, {"C1": {"1": calls}})

    # Both genotypes were already complete. The new phase spans two truth
    # heterozygotes, even though all three genotype rates have zero denominators.
    _check(
        _row(evaluate_reconstruction(result, truth), "C1", "1"),
        eligible_copies=4,
        observed_copies=4,
        reported_copies=4,
        known_copies=0,
        correct_copies=0,
        coverage=np.nan,
        consistency=np.nan,
        complete_known=0,
        complete_correct=0,
        complete_accuracy=np.nan,
        phase_eligible_heterozygotes=2,
        phased_heterozygotes=2,
        phase_coverage=1,
        phased_copies=4,
        phase_correct=4,
        phase_sets=1,
        phase_n50_cm=1,
    )


def test_position_alignment_and_summary_pool_counts_and_phase_sets():
    truth_one = _loci("1", [1_000_000, 3_000_000, 5_000_000], ["rs1", "rs2", "rs3"])
    call_one = _loci("1", [1_000_000, 3_000_000], ["rs4", "rs5"])
    loci_x = _loci("X", [10_000_000, 10_500_000, 11_000_000])
    pedigree = Pedigree([Person("C1", Sex.FEMALE)])
    truth = _truth(
        pedigree,
        {"1": truth_one, "X": loci_x},
        {"C1": {"1": [[A] * 3, [G] * 3], "X": [[A] * 3, [G] * 3]}},
    )
    observed_x = _calls(
        loci_x, [UNORDERED_AG, DIPLOID_MISSING, DIPLOID_MISSING], genotype=OBSERVED
    )
    calls_one = _calls(call_one, [AG, AG], phase_set=[0, 0])
    calls_x = _calls(
        loci_x, [GA, GA, GA], genotype=[OBSERVED, FORCED, FORCED], phase_set=[0, 0, 0]
    )
    result = _EvaluationInput(
        pedigree, {"C1": {"X": observed_x}}, {"C1": {"1": calls_one, "X": calls_x}}
    )
    evaluation = evaluate_reconstruction(result, truth)

    # Chromosome 1 calls align despite different rsIDs, leaving the third
    # truth site unknown. The X block independently takes the flipped order.
    _check(
        _row(evaluation, "C1", "1"),
        eligible_copies=6,
        known_copies=4,
        coverage=4 / 6,
        phase_n50_cm=2,
    )
    _check(
        _row(evaluation, "C1", "X"),
        eligible_copies=6,
        observed_copies=2,
        known_copies=4,
        coverage=1,
        phase_n50_cm=1,
    )

    # Pooling gives 8 new copies over 10 unobserved copies. The longer set
    # contains 2 of 5 phased loci, so the shorter span is the pooled N50.
    summary = evaluation.summary()
    assert "chrom" not in summary.columns
    assert list(summary["person"]) == ["C1"]
    _check(
        summary.iloc[0],
        eligible_copies=12,
        observed_copies=2,
        reported_copies=10,
        known_copies=8,
        correct_copies=8,
        coverage=8 / 10,
        consistency=1,
        complete_known=4,
        complete_correct=4,
        phase_sets=2,
        phase_n50_cm=1,
        phase_eligible_heterozygotes=6,
        phased_heterozygotes=5,
        phase_coverage=5 / 6,
        phased_copies=10,
        phase_correct=10,
        switch_opportunities=3,
        switch_errors=0,
    )


@pytest.mark.parametrize(
    ("label_sets", "a_is", "correct_copies"),
    [
        pytest.param([0, 0, 0], [-1, -1, -1], 8, id="one-shared-assignment"),
        pytest.param([0, 0, 1], [-1, -1, -1], 12, id="independent-label-sets"),
        pytest.param([-1, -1, -1], [-1, -1, -1], 12, id="independent-unlinked-sites"),
        pytest.param([0, 0, 0], [1, 1, 1], 4, id="fixed-a-is-father"),
    ],
)
def test_anonymous_couple_assigns_a_and_b_for_the_whole_label_set(
    label_sets, a_is, correct_copies
):
    loci = _loci("1", [1_000_000, 2_000_000, 3_000_000])
    pedigree = Pedigree(
        [
            Person("M1", Sex.FEMALE),
            Person("F1", Sex.MALE),
            Person("C1", Sex.MALE, mother="M1", father="F1"),
        ]
    )
    truth = _truth(
        pedigree,
        {"1": loci},
        {
            "M1": {"1": [[A] * 3, [A] * 3]},
            "F1": {"1": [[G] * 3, [G] * 3]},
            "C1": {"1": [[A] * 3, [G] * 3]},
        },
    )
    a = _calls(loci, [AA, AA, GG])
    b = _calls(loci, [GG, GG, AA])
    couple = _CoupleCalls(
        "M1",
        "F1",
        a,
        b,
        np.asarray(label_sets, dtype=np.int32),
        np.asarray(a_is, dtype=np.int32),
    )
    result = _EvaluationInput(pedigree, {}, {}, couples=(couple,))

    # Each site contributes four correct copies if its labels match and zero
    # otherwise. A as mother wins the first two sites and loses the third.
    _check(
        _row(evaluate_reconstruction(result, truth), "M1+F1", "1"),
        eligible_copies=12,
        reported_copies=12,
        known_copies=12,
        correct_copies=correct_copies,
        consistency=correct_copies / 12,
        coverage=1,
        complete_known=6,
        complete_correct=correct_copies // 2,
        complete_accuracy=correct_copies / 12,
    )
