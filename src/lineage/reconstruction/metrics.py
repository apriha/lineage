"""Score reconstruction calls against the retained simulation panel."""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterator, Mapping
from dataclasses import dataclass
from types import MappingProxyType
from typing import NamedTuple, TypeAlias

import numpy as np
import numpy.typing as npt
import pandas as pd

from lineage.reconstruction.encoding import (
    _BASES,
    _ESTIMATED,
    _FORCED,
    _GENOTYPE,
    _HOLDS,
    _ORIGIN,
    _Calls,
    _copies,
    _copies_of,
    _natural_key,
    _pair,
    haplotypes,
)
from lineage.reconstruction.pedigree import Pedigree
from lineage.reconstruction.simulate import _readonly, _SimulatedArrays
from lineage.types import UNKNOWN, Haplotypes

_Blocks: TypeAlias = Mapping[str, Mapping[str, _Calls]]
_Indices: TypeAlias = npt.NDArray[np.intp]
_Labels: TypeAlias = npt.NDArray[np.int32]
_Bytes: TypeAlias = npt.NDArray[np.int8]
_Counts: TypeAlias = npt.NDArray[np.int64]
_Bools: TypeAlias = npt.NDArray[np.bool_]
_Positions: TypeAlias = npt.NDArray[np.uint32]
_Centimorgans: TypeAlias = npt.NDArray[np.float64]
_PhaseLength: TypeAlias = tuple[float, int]


@dataclass(frozen=True, slots=True, eq=False)
class _CoupleCalls:
    """Genotype calls labelled A and B for two parents of the same child.

    Each label refers to one parent throughout its label set.
    ``a_is`` identifies A as the mother (0), father (1), or unknown (-1).
    """

    mother: str
    father: str
    a: _Calls
    b: _Calls
    label_set: _Labels
    a_is: _Bytes

    def __post_init__(self) -> None:
        chrom = self.a.loci.chrom
        if (
            chrom == "X"
            or chrom != self.b.loci.chrom
            or not np.array_equal(self.a.loci.pos, self.b.loci.pos)
        ):
            raise ValueError(
                "Couple calls must share an autosomal chromosome and positions"
            )
        n = len(self.a.loci.pos)
        labels = np.asarray(self.label_set)
        identity = np.asarray(self.a_is)
        if (
            labels.dtype.kind not in "iu"
            or labels.shape != (n,)
            or np.any((labels < -1) | (labels > np.iinfo(np.int32).max))
        ):
            raise ValueError(
                "label_set must give an int32 label or -1 for each couple locus"
            )
        if (
            identity.dtype.kind not in "iu"
            or identity.shape != (n,)
            or np.any((identity < -1) | (identity > 1))
        ):
            raise ValueError("a_is must give 0, 1 or -1 for each couple locus")
        mothers = labels[(labels >= 0) & (identity == 0)]
        fathers = labels[(labels >= 0) & (identity == 1)]
        if np.intersect1d(mothers, fathers).size:
            raise ValueError("a_is assigns both parents to A within a label set")
        for block in (self.a, self.b):
            linked = np.flatnonzero(block.phase_set >= 0)
            ordered = linked[np.argsort(block.phase_set[linked], kind="stable")]
            earlier, later = ordered[:-1], ordered[1:]
            same_phase = block.phase_set[earlier] == block.phase_set[later]
            crosses = (labels[earlier] < 0) | (labels[earlier] != labels[later])
            if np.any(same_phase & crosses):
                raise ValueError(
                    "Each multi-locus anonymous phase set must stay within "
                    "one linked label_set"
                )
        object.__setattr__(
            self, "label_set", _readonly(labels.astype(np.int32, copy=False))
        )
        object.__setattr__(
            self, "a_is", _readonly(identity.astype(np.int8, copy=False))
        )


@dataclass(frozen=True, slots=True, eq=False)
class _EvaluationInput:
    """Pedigree, raw observations and final named and anonymous calls.

    The producer supplies distinct named and anonymous facts. Mapping layers
    and couple metadata are owned here. Run-owned call blocks are shared as
    immutable inputs.
    """

    pedigree: Pedigree
    observations: _Blocks
    calls: _Blocks
    couples: tuple[_CoupleCalls, ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "observations", _freeze_blocks(self.pedigree, self.observations)
        )
        object.__setattr__(self, "calls", _freeze_blocks(self.pedigree, self.calls))
        couples = tuple(self.couples)
        names = set(self.pedigree.names)
        for couple in couples:
            if (
                couple.mother == couple.father
                or couple.mother not in names
                or couple.father not in names
            ):
                raise ValueError(
                    "A couple must name two distinct people in the pedigree"
                )
        object.__setattr__(self, "couples", couples)


def _freeze_blocks(pedigree: Pedigree, blocks: _Blocks) -> _Blocks:
    names = set(pedigree.names)
    frozen: dict[str, Mapping[str, _Calls]] = {}
    for person, chromosomes in blocks.items():
        if person not in names:
            raise ValueError(f"Calls name {person!r}, who is missing from the pedigree")
        stored = dict(chromosomes)
        for chrom, block in stored.items():
            if chrom != block.loci.chrom:
                raise ValueError(
                    f"Calls for {person!r} keyed by {chrom!r} contain chromosome "
                    f"{block.loci.chrom!r}"
                )
        frozen[person] = MappingProxyType(stored)
    return MappingProxyType(frozen)


_COLUMNS = (
    "eligible_copies",
    "observed_copies",
    "reported_copies",
    "known_copies",
    "correct_copies",
    "coverage",
    "consistency",
    "forced_known",
    "forced_correct",
    "estimated_known",
    "estimated_correct",
    "complete_known",
    "complete_correct",
    "complete_accuracy",
    "phased_copies",
    "phase_correct",
    "phase_sets",
    "phase_n50_cm",
    "phase_eligible_heterozygotes",
    "phased_heterozygotes",
    "phase_coverage",
    "switch_opportunities",
    "switch_errors",
    "origin_known",
    "origin_correct",
)
_RATES = {
    "consistency": ("correct_copies", "known_copies"),
    "complete_accuracy": ("complete_correct", "complete_known"),
    "phase_coverage": ("phased_heterozygotes", "phase_eligible_heterozygotes"),
}


class _Score(NamedTuple):
    counts: Counter[str]
    lengths: list[_PhaseLength]


def _ratio(num: int, den: int) -> float:
    return num / den if den else np.nan


def _n50(lengths: list[_PhaseLength]) -> float:
    """Return the cM span at which the longest phase sets cover half the phased sites."""
    remaining = sum(weight for _, weight in lengths)
    for span, weight in sorted(lengths, reverse=True):
        remaining -= 2 * weight
        if remaining <= 0:
            return span
    return np.nan


def _frame(
    scores: Mapping[tuple[str, ...], _Score], keys: tuple[str, ...]
) -> pd.DataFrame:
    columns: dict[str, pd.Series[str] | pd.Series[int] | pd.Series[float]] = {
        key: pd.Series([row[i] for row in scores], dtype="string")
        for i, key in enumerate(keys)
    }
    for column in _COLUMNS:
        if column == "coverage":
            values = [
                _ratio(c["known_copies"], c["eligible_copies"] - c["observed_copies"])
                for c, _ in scores.values()
            ]
        elif column == "phase_n50_cm":
            values = [_n50(lengths) for _, lengths in scores.values()]
        elif column in _RATES:
            num, den = _RATES[column]
            values = [_ratio(c[num], c[den]) for c, _ in scores.values()]
        else:
            columns[column] = pd.Series(
                [c[column] for c, _ in scores.values()], dtype="int64"
            )
            continue
        columns[column] = pd.Series(values, dtype="float64")
    return pd.DataFrame(columns, index=pd.RangeIndex(len(scores)))


@dataclass(frozen=True, slots=True, init=False, eq=False)
class Evaluation:
    """Genotype, phase and origin scores by person and chromosome."""

    _table: pd.DataFrame
    _summary: pd.DataFrame

    def __init__(self, scores: Mapping[tuple[str, ...], _Score]) -> None:
        pooled: dict[tuple[str, ...], _Score] = {}
        for row, (counts, lengths) in scores.items():
            total = pooled.setdefault((row[0],), _Score(Counter(), []))
            total.counts.update(counts)
            total.lengths.extend(lengths)
        object.__setattr__(self, "_table", _frame(scores, ("person", "chrom")))
        object.__setattr__(self, "_summary", _frame(pooled, ("person",)))

    @property
    def table(self) -> pd.DataFrame:
        """Return the per-person, per-chromosome scores."""
        return self._table.copy(deep=False)

    def summary(self) -> pd.DataFrame:
        """Return per-person scores pooled across chromosomes."""
        return self._summary.copy(deep=False)


def _group_rows(labels: _Labels, rows: _Indices) -> tuple[_Indices, _Indices]:
    ordered = rows[np.argsort(labels[rows], kind="stable")]
    bounds = np.r_[0, np.flatnonzero(np.diff(labels[ordered])) + 1, len(ordered)]
    return ordered, bounds


def _groups(labels: _Labels, selected: _Bools) -> Iterator[_Indices]:
    ordered, bounds = _group_rows(labels, np.flatnonzero(selected & (labels >= 0)))
    for group in map(int, np.flatnonzero(np.diff(bounds) >= 2)):
        yield ordered[bounds[group] : bounds[group + 1]]


def _match_positions(
    panel: _Positions, block: _Calls | None
) -> tuple[_Indices, _Indices]:
    if block is None:
        empty = np.empty(0, dtype=np.intp)
        return empty, empty
    _, panel_rows, call_rows = np.intersect1d(
        panel, block.loci.pos, assume_unique=True, return_indices=True
    )
    return panel_rows, call_rows


def _dosage(allowed: npt.NDArray[np.uint32]) -> _Bytes:
    """Find the minimum count of each base across allowed genotypes at each site."""
    dosage = np.empty((len(_BASES), len(allowed)), dtype=np.int8)
    for base in range(len(_BASES)):
        carries = _HOLDS[0][base] | _HOLDS[1][base]
        dose = (allowed & carries) == allowed
        dosage[base] = dose.astype(np.int8) + (allowed == _pair(base, base))
    return dosage


def _truth_dosage(target: Haplotypes) -> _Bytes:
    return np.stack(
        [(target == base).sum(axis=0, dtype=np.int8) for base in range(len(_BASES))]
    )


class _Reported(NamedTuple):
    dosage: _Bytes
    ordered: Haplotypes
    phase_set: _Labels
    kinds: npt.NDArray[np.uint8]

    def take(self, rows: _Indices) -> _Reported:
        return _Reported(
            self.dosage[:, rows],
            self.ordered[:, rows],
            self.phase_set[rows],
            self.kinds[rows],
        )


def _aligned_dosage(
    block: _Calls | None, n: int, dest: _Indices, source: _Indices
) -> _Bytes:
    dosage = np.zeros((len(_BASES), n), dtype=np.int8)
    if block is not None:
        dosage[:, dest] = _dosage(block.allowed[source])
    return dosage


def _reported(block: _Calls | None, panel: _Positions) -> _Reported:
    n = len(panel)
    dest, source = _match_positions(panel, block)
    dosage = _aligned_dosage(block, n, dest, source)
    ordered = np.full((2, n), UNKNOWN, dtype=np.int8)
    sets = np.full(n, -1, dtype=np.int32)
    kinds = np.zeros((n, 3), dtype=np.uint8)
    if block is not None:
        ordered[:, dest] = haplotypes(block)[:, source]
        sets[dest] = block.phase_set[source]
        kinds[dest] = block.kinds[source]
    return _Reported(dosage, ordered, sets, kinds)


def _observation(block: _Calls | None, panel: _Positions) -> _Bytes:
    dest, source = _match_positions(panel, block)
    return _aligned_dosage(block, len(panel), dest, source)


class _Genotype(NamedTuple):
    new: _Counts
    correct: _Counts
    complete: _Bools
    matches: _Bools


def _genotype(
    observed: _Bytes, reported: _Bytes, target: _Bytes, copies: int
) -> _Genotype:
    new = np.maximum(reported - observed, 0)
    correct = np.minimum(new, np.maximum(target - observed, 0))
    complete = (reported.sum(axis=0) == copies) & (observed.sum(axis=0) < copies)
    return _Genotype(
        new.sum(axis=0), correct.sum(axis=0), complete, (reported == target).all(axis=0)
    )


def _orientation(reported: Haplotypes, target: Haplotypes) -> tuple[_Counts, _Counts]:
    present = reported >= 0
    given = (present & (reported == target)).sum(axis=0)
    reversed_ = (present & (reported == target[::-1])).sum(axis=0)
    return given, reversed_


def _phase(
    score: _Score,
    reported: _Reported,
    target: Haplotypes,
    matches: _Bools,
    cm: _Centimorgans,
    heterozygous: _Bools,
) -> None:
    counts, lengths = score
    positional = (reported.ordered >= 0).sum(axis=0)
    direct, flipped = _orientation(reported.ordered, target)
    for group in _groups(reported.phase_set, positional > 0):
        counts["phased_copies"] += int(positional[group].sum())
        counts["phase_correct"] += max(
            int(direct[group].sum()), int(flipped[group].sum())
        )
        lengths.append((float(cm[group[-1]] - cm[group[0]]), len(group)))
    counts["phase_sets"] = len(lengths)

    complete_het = (reported.dosage > 0).sum(axis=0) == 2
    for group in _groups(
        reported.phase_set, heterozygous & complete_het & (positional == 2)
    ):
        counts["phased_heterozygotes"] += len(group)
        informative = group[matches[group]]
        if len(informative) >= 2:
            agrees_with_truth_order = direct[informative] == 2
            counts["switch_opportunities"] += len(informative) - 1
            counts["switch_errors"] += int(
                np.count_nonzero(np.diff(agrees_with_truth_order))
            )


def _score(
    reported: _Reported,
    observed: _Bytes,
    target: Haplotypes,
    copies: int,
    cm: _Centimorgans,
    *,
    origin: bool,
) -> _Score:
    genotype = _genotype(observed, reported.dosage, _truth_dosage(target), copies)
    new, correct, complete, matches = genotype
    counts: Counter[str] = Counter(
        eligible_copies=copies * target.shape[1],
        observed_copies=int(observed.sum()),
        reported_copies=int(reported.dosage.sum()),
        known_copies=int(new.sum()),
        correct_copies=int(correct.sum()),
        complete_known=int(complete.sum()),
        complete_correct=int((complete & matches).sum()),
    )
    for kind, label in ((_FORCED, "forced"), (_ESTIMATED, "estimated")):
        selected = reported.kinds[:, _GENOTYPE] == kind
        counts[f"{label}_known"] = int(new[selected].sum())
        counts[f"{label}_correct"] = int(correct[selected].sum())
    score = _Score(counts, [])
    heterozygous = (copies == 2) & (target[0] != target[1])
    counts["phase_eligible_heterozygotes"] = int(heterozygous.sum())
    if copies == 2:
        _phase(score, reported, target, matches, cm, heterozygous)
    if origin:
        present = reported.ordered >= 0
        known_origin = (copies == 1) | (
            heterozygous & (reported.kinds[:, _ORIGIN] != 0)
        )
        known = present.any(axis=0) & known_origin
        correct_origin = ((~present) | (reported.ordered == target)).all(axis=0)
        counts["origin_known"] = int(known.sum())
        counts["origin_correct"] = int((known & correct_origin).sum())
    return score


def _phase_counts(
    reported: _Reported, target: Haplotypes, matches: _Bools, cm: _Centimorgans
) -> Counter[str]:
    score: _Score = _Score(Counter(), [])
    _phase(score, reported, target, matches, cm, target[0] != target[1])
    return score.counts


def _phase_key(counts: Counter[str]) -> tuple[int, int, int, int]:
    return (
        counts["phase_correct"],
        counts["phased_heterozygotes"],
        counts["switch_opportunities"],
        -counts["switch_errors"],
    )


def _choose_mother(*differences: _Counts) -> _Bools:
    """Choose A as mother when the first nonzero difference is positive, or on a tie."""
    choice = np.ones(len(differences[0]), dtype=np.bool_)
    for differences_by_key in reversed(differences):
        choice = np.where(differences_by_key == 0, choice, differences_by_key > 0)
    return choice


def _assign_couple(
    couple: _CoupleCalls,
    a: _Reported,
    b: _Reported,
    mother: Haplotypes,
    father: Haplotypes,
    panel: _Positions,
    empty: _Bytes,
    cm: _Centimorgans,
) -> _Bools:
    dest, source = _match_positions(panel, couple.a)
    label = couple.label_set[source]
    identity = couple.a_is[source]

    # An identity anchor outside the panel still fixes its entire label set.
    anchors = (couple.label_set >= 0) & (couple.a_is >= 0)
    mother_labels = couple.label_set[anchors & (couple.a_is == 0)]
    father_labels = couple.label_set[anchors & (couple.a_is == 1)]
    identity[np.isin(label, mother_labels)] = 0
    identity[np.isin(label, father_labels)] = 1

    truth_m = _truth_dosage(mother)
    truth_f = _truth_dosage(father)
    a_mother = _genotype(empty, a.dosage, truth_m, 2)
    a_father = _genotype(empty, a.dosage, truth_f, 2)
    b_mother = _genotype(empty, b.dosage, truth_m, 2)
    b_father = _genotype(empty, b.dosage, truth_f, 2)
    a_correct = a_mother.correct[dest] - a_father.correct[dest]
    b_correct = b_father.correct[dest] - b_mother.correct[dest]
    correct = a_correct + b_correct
    complete = (
        (a_mother.complete & a_mother.matches).astype(np.int64)
        + (b_father.complete & b_father.matches)
        - (a_father.complete & a_father.matches)
        - (b_mother.complete & b_mother.matches)
    )[dest]
    a_kinds = a.kinds[dest, _GENOTYPE]
    b_kinds = b.kinds[dest, _GENOTYPE]
    forced = a_correct * (a_kinds == _FORCED) + b_correct * (b_kinds == _FORCED)
    estimated = a_correct * (a_kinds == _ESTIMATED) + b_correct * (
        b_kinds == _ESTIMATED
    )
    # Unlinked loci contribute genotype scores only. With every key equal, either
    # assignment contributes the same counts.
    to_mother = _choose_mother(correct, complete, forced, estimated)

    unresolved = np.flatnonzero((label >= 0) & (identity == -1))
    if len(unresolved):
        ordered, bounds = _group_rows(label, unresolved)
        starts = bounds[:-1]
        widths = np.diff(bounds)
        correct_totals = np.add.reduceat(correct[ordered], starts)
        complete_totals = np.add.reduceat(complete[ordered], starts)
        choices = _choose_mother(
            correct_totals,
            complete_totals,
            np.add.reduceat(forced[ordered], starts),
            np.add.reduceat(estimated[ordered], starts),
        )
        # Phase sets belong to one label set, so a linked assignment can be
        # compared using the same phase scoring as the final evaluation.
        tied = np.flatnonzero(
            (correct_totals == 0) & (complete_totals == 0) & (widths >= 2)
        )
        for group in map(int, tied):
            rows = dest[ordered[bounds[group] : bounds[group + 1]]]
            if not np.any((a.phase_set[rows] >= 0) | (b.phase_set[rows] >= 0)):
                continue
            a_group, b_group = a.take(rows), b.take(rows)
            m_group, f_group = mother[:, rows], father[:, rows]
            distances = cm[rows]
            direct = _phase_counts(a_group, m_group, a_mother.matches[rows], distances)
            direct.update(
                _phase_counts(b_group, f_group, b_father.matches[rows], distances)
            )
            swapped = _phase_counts(a_group, f_group, a_father.matches[rows], distances)
            swapped.update(
                _phase_counts(b_group, m_group, b_mother.matches[rows], distances)
            )
            direct_key, swapped_key = _phase_key(direct), _phase_key(swapped)
            if direct_key != swapped_key:
                choices[group] = direct_key > swapped_key
        to_mother[ordered] = np.repeat(choices, widths)
    to_mother[identity >= 0] = identity[identity >= 0] == 0

    assignment = np.ones(len(panel), dtype=np.bool_)
    assignment[dest] = to_mother
    return assignment


def _couple_score(
    couple: _CoupleCalls,
    mother: Haplotypes,
    father: Haplotypes,
    panel: _Positions,
    cm: _Centimorgans,
) -> _Score:
    a = _reported(couple.a, panel)
    b = _reported(couple.b, panel)
    empty = _observation(None, panel)
    to_mother = _assign_couple(couple, a, b, mother, father, panel, empty, cm)
    a_truth = np.where(to_mother, mother, father)
    b_truth = np.where(to_mother, father, mother)
    a_score = _score(a, empty, a_truth, 2, cm, origin=False)
    b_score = _score(b, empty, b_truth, 2, cm, origin=False)
    return _Score(a_score.counts + b_score.counts, a_score.lengths + b_score.lengths)


def _check_ploidy(
    person: str,
    chrom: str,
    copies: int,
    result: _EvaluationInput,
    block: _Calls | None,
    raw: _Calls | None,
) -> None:
    expected = _copies(chrom, result.pedigree[person].sex)
    for value in (
        expected,
        *(int(_copies_of(c.allowed)[0]) for c in (block, raw) if c is not None),
    ):
        if value is not None and value != copies:
            raise ValueError(
                f"Ploidy mismatch for {person!r} on chromosome {chrom!r}: "
                f"reconstruction has {value} copies, truth has {copies}"
            )


def evaluate_reconstruction(
    result: _EvaluationInput, truth: _SimulatedArrays
) -> Evaluation:
    """Score reconstruction genotypes, phases and origins against a simulation.

    Parameters
    ----------
    result : _EvaluationInput
        The pedigree, raw test observations and final calls.
    truth : _SimulatedArrays
        The simulation's retained panel, map and inherited copies.

    Returns
    -------
    Evaluation
        Per-person and chromosome scores with a pooled summary.

    Raises
    ------
    ValueError
        If reconstruction and truth differ in their people or ploidy.
    """
    names = set(result.pedigree.names)
    truth_names = set(truth.pedigree.names)
    if names != truth_names:
        raise ValueError(
            f"Reconstruction and truth must name the same people: "
            f"only in reconstruction {sorted(names - truth_names)}, "
            f"only in truth {sorted(truth_names - names)}"
        )

    positions = {chrom: panel.loci.pos for chrom, panel in truth.panel.items()}
    cm = {chrom: truth.genetic_map.cm(chrom, pos) for chrom, pos in positions.items()}
    scores: dict[tuple[str, ...], _Score] = {}
    for person in names:
        for chrom, pos in positions.items():
            block = result.calls.get(person, {}).get(chrom)
            raw = result.observations.get(person, {}).get(chrom)
            target = truth.haplotypes[person][chrom]
            copies = _copies(chrom, truth.pedigree[person].sex) or 2
            _check_ploidy(person, chrom, copies, result, block, raw)
            scores[person, chrom] = _score(
                _reported(block, pos),
                _observation(raw, pos),
                target,
                copies,
                cm[chrom],
                origin=True,
            )
    couple_names: dict[str, tuple[str, str]] = {}
    for couple in result.couples:
        chrom = couple.a.loci.chrom
        if chrom not in positions:
            continue
        parents = (couple.mother, couple.father)
        name = "+".join(parents)
        if name in names or couple_names.setdefault(name, parents) != parents:
            raise ValueError(
                f"Anonymous couple name {name!r} is ambiguous; use distinct pedigree "
                "names that yield unique couple names"
            )
        key = (name, chrom)
        if key in scores:
            raise ValueError(
                f"Multiple blocks for couple {name!r} on chromosome {chrom!r}; "
                "provide one block per couple and chromosome"
            )
        scores[key] = _couple_score(
            couple,
            truth.haplotypes[couple.mother][chrom],
            truth.haplotypes[couple.father][chrom],
            positions[chrom],
            cm[chrom],
        )
    ordered = sorted(scores, key=lambda key: (key[0], _natural_key(key[1])))
    return Evaluation({key: scores[key] for key in ordered})
