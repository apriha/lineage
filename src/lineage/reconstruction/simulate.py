"""Synthetic pedigree truth and unphased observations."""

from __future__ import annotations

import hashlib
from collections.abc import Collection, Mapping, Sequence
from dataclasses import dataclass
from enum import IntEnum
from types import MappingProxyType
from typing import Literal, TypeAlias, TypeVar

import numpy as np
import numpy.typing as npt

from lineage import _grch37
from lineage.maps import GeneticMap
from lineage.reconstruction.pedigree import Pedigree, _known_names
from lineage.types import ABSENT, UNKNOWN, Haplotypes, Loci, Sex

__all__ = ["Observation", "Panel", "simulate_pedigree", "synthetic_panel"]

_Bytes: TypeAlias = npt.NDArray[np.int8]
_Floats: TypeAlias = npt.NDArray[np.float64]
_Bools: TypeAlias = npt.NDArray[np.bool_]
_Model: TypeAlias = Literal["independent", "correlated"]
_Scalar = TypeVar("_Scalar", bound=np.generic)
_STRETCH = 256
_TEMPLATES = 12


def _readonly(array: npt.NDArray[_Scalar]) -> npt.NDArray[_Scalar]:
    return np.frombuffer(array.tobytes(), dtype=array.dtype).reshape(array.shape)


def _probabilities(value: object, label: str) -> _Floats:
    items = np.asarray(value)
    if type(value) is int:
        valid = 0 <= value <= 1
    else:
        if items.dtype.kind not in "iuf" or _grch37.holds_bool(value):
            raise TypeError(f"{label} must contain real numbers")
        valid = bool(np.all(np.isfinite(items) & (items >= 0) & (items <= 1)))
    if not valid:
        raise ValueError(f"{label} must be finite and between 0 and 1")
    return np.asarray(items, dtype=np.float64)


@dataclass(frozen=True, init=False, eq=False, slots=True)
class Panel:
    """A synthetic biallelic marker panel for one chromosome.

    Parameters
    ----------
    loci : Loci
        The marker axis.
    alleles : array_like of int, shape (2, n)
        Distinct A, C, G and T codes (0 to 3) at each locus.
    frequency : array_like of float, shape (n,)
        The probability of the second allele at each locus.
    """

    loci: Loci
    alleles: _Bytes
    frequency: _Floats

    def __init__(
        self, loci: Loci, alleles: npt.ArrayLike, frequency: npt.ArrayLike
    ) -> None:
        if not isinstance(loci, Loci):
            raise TypeError("loci must be a Loci")
        n = len(loci.pos)
        codes = np.asarray(alleles)
        if codes.shape != (2, n):
            raise ValueError(f"alleles must have shape (2, {n})")
        if codes.dtype.kind not in "iuf" or _grch37.holds_bool(alleles):
            raise ValueError("alleles must contain numeric, non-boolean codes")
        checked = _grch37.check_whole_numbers(codes.ravel(), "alleles").reshape(2, n)
        if np.any((checked < 0) | (checked > 3)):
            raise ValueError("alleles must contain codes from 0 to 3")
        if np.any(checked[0] == checked[1]):
            raise ValueError("alleles must be distinct at each locus")
        frequencies = _probabilities(frequency, "frequency")
        if frequencies.shape != (n,):
            raise ValueError(f"frequency must have shape ({n},)")
        object.__setattr__(self, "loci", loci)
        object.__setattr__(self, "alleles", _readonly(checked.astype(np.int8)))
        object.__setattr__(self, "frequency", _readonly(frequencies))


@dataclass(frozen=True, slots=True)
class Observation:
    """Per-site assay and call probabilities.

    Parameters
    ----------
    error_rate : float, default 0
        Probability of reporting a different genotype.
    missing_rate : float, default 0
        Probability of a no-call at an assayed site.
    chip_keep : float, default 1
        Probability of assaying a site.
    """

    error_rate: float = 0.0
    missing_rate: float = 0.0
    chip_keep: float = 1.0

    def __post_init__(self) -> None:
        for name in ("error_rate", "missing_rate", "chip_keep"):
            value = getattr(self, name)
            if not isinstance(value, (int, float, np.integer, np.floating)):
                raise TypeError(f"{name} must be a real number")
            object.__setattr__(self, name, float(_probabilities(value, name)))


_DEFAULT_OBSERVATION = Observation()


@dataclass(frozen=True, eq=False, slots=True)
class _ObservedCalls:
    """Unphased allele calls and the assay mask on one panel axis."""

    calls: Haplotypes
    assayed: _Bools


@dataclass(frozen=True, eq=False, slots=True)
class _SimulatedArrays:
    """A simulation's panel, inherited copies, unphased calls and run parameters."""

    pedigree: Pedigree
    panel: Mapping[str, Panel]
    genetic_map: GeneticMap
    haplotypes: Mapping[str, Mapping[str, Haplotypes]]
    transmissions: Mapping[tuple[str, str], Mapping[str, _Bytes]]
    observations: Mapping[str, Mapping[str, _ObservedCalls]]
    parameters: Mapping[str, object]


def _integer(value: object, label: str, minimum: int) -> int:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)):
        raise TypeError(f"{label} must be an integer")
    checked = int(value)
    if checked < minimum:
        raise ValueError(f"{label} must be at least {minimum}")
    return checked


def _seed(seed: int | None) -> int:
    if seed is None:
        entropy = np.random.SeedSequence().entropy
        assert isinstance(entropy, int)
        return entropy
    return _integer(seed, "seed", 0)


class _Purpose(IntEnum):
    POSITIONS = 1
    ALLELES = 2
    FREQUENCY = 3
    TEMPLATES = 4
    FOUNDER = 5
    MEIOSIS = 6
    ERROR = 7
    ASSAY = 8
    MISSING = 9


def _rng(seed: int, purpose: _Purpose, *keys: str) -> np.random.Generator:
    words = [int(purpose)]
    for key in keys:
        digest = hashlib.sha256(key.encode("utf-8", errors="surrogatepass"))
        words.extend(map(int, np.frombuffer(digest.digest(), dtype="<u4")))
    return np.random.default_rng(np.random.SeedSequence(seed, spawn_key=tuple(words)))


def _check_region(region: str, genetic_map: GeneticMap) -> None:
    if not isinstance(region, str):
        raise TypeError("chromosomes must be str")
    if region not in _grch37.LENGTHS:
        raise ValueError(f"Unsupported chromosome {region!r}; use 1 to 22 or X")
    if region not in genetic_map.tables:
        raise ValueError(f"The genetic map has no region {region!r}")


def _positions(
    rng: np.random.Generator, region: str, first: int, last: int, n: int
) -> npt.NDArray[np.int64]:
    gaps = []
    if region == "X":
        for start, end in (_grch37.PAR1, _grch37.PAR2):
            begin, finish = max(start, first), min(end, last)
            if begin <= finish:
                gaps.append((begin, finish - begin + 1))
    capacity = last - first + 1 - sum(length for _, length in gaps)
    if n > capacity:
        raise ValueError(f"Chromosome {region!r} has room for {capacity} loci, not {n}")
    positions = np.sort(rng.choice(capacity, size=n, replace=False)) + first
    for start, length in gaps:
        positions[positions >= start] += length
    return positions


def synthetic_panel(
    genetic_map: GeneticMap,
    *,
    count: int,
    chromosomes: Sequence[str] | None = None,
    seed: int | None = None,
) -> Mapping[str, Panel]:
    """Draw a biallelic panel with evenly allocated markers in map spans.

    Parameters
    ----------
    genetic_map : GeneticMap
        The map whose regions and endpoints define the sampling spans.
    count : int
        Total number of markers. Extra markers go to the earliest chromosomes.
    chromosomes : sequence of str, optional
        Eligible chromosomes. Defaults to all mapped autosomes and non-PAR X.
    seed : int, optional
        Nonnegative seed for reproducible draws.
    """
    if not isinstance(genetic_map, GeneticMap):
        raise TypeError("genetic_map must be a GeneticMap")
    count = _integer(count, "count", 1)
    run_seed = _seed(seed)
    if chromosomes is None:
        chromosomes = tuple(c for c in genetic_map.regions if c in _grch37.LENGTHS)
    elif isinstance(chromosomes, str) or not isinstance(chromosomes, Sequence):
        raise TypeError("chromosomes must be a sequence of str")
    if not chromosomes:
        raise ValueError("Select at least one supported chromosome")
    seen: set[str] = set()
    for region in chromosomes:
        _check_region(region, genetic_map)
        if region in seen:
            raise ValueError(f"Chromosome {region!r} is repeated")
        seen.add(region)

    size, extra = divmod(count, len(seen))
    result: dict[str, Panel] = {}
    for index, region in enumerate(c for c in _grch37.LENGTHS if c in seen):
        n = size + (index < extra)
        if not n:
            continue
        pos = genetic_map.tables[region]["pos"]
        positions = _positions(
            _rng(run_seed, _Purpose.POSITIONS, region),
            region,
            int(pos.iloc[0]),
            int(pos.iloc[-1]),
            n,
        )
        alleles = _rng(run_seed, _Purpose.ALLELES, region)
        first = alleles.integers(4, size=n, dtype=np.int8)
        second = alleles.integers(3, size=n, dtype=np.int8)
        frequency = _rng(run_seed, _Purpose.FREQUENCY, region).uniform(
            0.05, 0.5, size=n
        )
        result[region] = Panel(
            Loci(region, positions),
            np.stack((first, second + (second >= first))),
            frequency,
        )
    return MappingProxyType(result)


def _select_tested(
    pedigree: Pedigree, tested: Collection[str] | None
) -> tuple[str, ...]:
    if tested is None:
        return pedigree.names
    if isinstance(tested, str) or not isinstance(tested, Collection):
        raise TypeError("tested must be a collection of str")
    selected: set[str] = set()
    for name in tested:
        if not isinstance(name, str):
            raise TypeError("tested must contain str names")
        try:
            pedigree[name]
        except KeyError as error:
            raise ValueError(error.args[0]) from None
        if name in selected:
            known = _known_names(pedigree.names)
            raise ValueError(f"Tested person {name!r} is repeated; {known}")
        selected.add(name)
    return tuple(name for name in pedigree.names if name in selected)


def _settings(
    pedigree: Pedigree, observation: Observation | Mapping[str, Observation]
) -> Mapping[str, Observation]:
    if isinstance(observation, Observation):
        return MappingProxyType(dict.fromkeys(pedigree.names, observation))
    if not isinstance(observation, Mapping):
        raise TypeError("observation must be an Observation or a mapping")
    for name, value in observation.items():
        if not isinstance(name, str):
            raise TypeError("observation keys must be str names")
        try:
            pedigree[name]
        except KeyError as error:
            raise ValueError(error.args[0]) from None
        if not isinstance(value, Observation):
            raise TypeError(f"The observation of {name!r} must be an Observation")
    return MappingProxyType(
        {name: observation.get(name, _DEFAULT_OBSERVATION) for name in pedigree.names}
    )


def _panels(panel: Mapping[str, Panel], genetic_map: GeneticMap) -> Mapping[str, Panel]:
    if not isinstance(panel, Mapping):
        raise TypeError("panel must be a mapping of chromosome to Panel")
    for region, block in panel.items():
        _check_region(region, genetic_map)
        if not isinstance(block, Panel):
            raise TypeError(f"The panel of region {region!r} must be a Panel")
        if block.loci.chrom != region:
            raise ValueError(f"The panel keyed {region!r} is on {block.loci.chrom!r}")
    return MappingProxyType({c: panel[c] for c in _grch37.LENGTHS if c in panel})


def _meiosis(
    rng: np.random.Generator, cm: _Floats, *, start_cm: float, end_cm: float
) -> _Bytes:
    """Draw the parental copy at each marker, switching at crossover events."""
    if not len(cm):
        return np.empty(0, dtype=np.int8)
    count = rng.poisson((end_cm - start_cm) / 100)
    events = np.sort(rng.uniform(start_cm, end_cm, size=count))
    initial = rng.integers(2)
    return ((initial + np.searchsorted(events, cm, side="right")) % 2).astype(np.int8)


def _templates(rng: np.random.Generator, frequency: _Floats) -> _Bools:
    """Draw shared allele choices for the correlated founder model."""
    return rng.random((_TEMPLATES, len(frequency))) < frequency


def _founder_copy(
    rng: np.random.Generator, panel: Panel, templates: _Bools | None
) -> _Bytes:
    n = len(panel.frequency)
    markers = np.arange(n)
    if templates is None:
        alleles = rng.random(n) < panel.frequency
    else:
        choice = rng.integers(len(templates), size=(n + _STRETCH - 1) // _STRETCH)
        alleles = templates[choice[markers // _STRETCH], markers]
    return panel.alleles[alleles.astype(np.intp), markers]


def _inherit(
    pedigree: Pedigree,
    region: str,
    panel: Panel,
    genetic_map: GeneticMap,
    seed: int,
    founders: _Model,
) -> tuple[dict[str, Haplotypes], dict[tuple[str, str], _Bytes]]:
    table = genetic_map.tables[region]["map"]
    start, end = float(table.iloc[0]), float(table.iloc[-1])
    cm = genetic_map.cm(region, panel.loci.pos)
    n = len(cm)
    markers = np.arange(n)
    shared = None
    if founders == "correlated":
        shared = _templates(_rng(seed, _Purpose.TEMPLATES, region), panel.frequency)

    truth: dict[str, Haplotypes] = {}
    transmissions: dict[tuple[str, str], _Bytes] = {}
    for name in pedigree.order:
        person = pedigree[name]
        copies = np.full((2, n), ABSENT, dtype=np.int8)
        for row, parent in enumerate((person.mother, person.father)):
            paternal_x = region == "X" and row == 1
            absent = paternal_x and person.sex == Sex.MALE
            if parent is None:
                if not absent:
                    copies[row] = _founder_copy(
                        _rng(seed, _Purpose.FOUNDER, region, name, str(row)),
                        panel,
                        shared,
                    )
                continue
            if paternal_x:
                inherited = np.full(n, ABSENT if absent else 0, dtype=np.int8)
            else:
                inherited = _meiosis(
                    _rng(seed, _Purpose.MEIOSIS, region, parent, name),
                    cm,
                    start_cm=start,
                    end_cm=end,
                )
            if not absent:
                copies[row] = truth[parent][inherited, markers]
            transmissions[parent, name] = _readonly(inherited)
        truth[name] = _readonly(copies)
    return truth, transmissions


def _observe(
    truth: Haplotypes,
    panel: Panel,
    setting: Observation,
    seed: int,
    name: str,
    *,
    haploid: bool,
) -> _ObservedCalls:
    region = panel.loci.chrom
    n = len(panel.frequency)
    error_rng = _rng(seed, _Purpose.ERROR, region, name)
    wrong = error_rng.random(n) < setting.error_rate
    assayed = _rng(seed, _Purpose.ASSAY, region, name).random(n) < setting.chip_keep
    missing = (
        _rng(seed, _Purpose.MISSING, region, name).random(n) < setting.missing_rate
    )

    first, second = panel.alleles
    if haploid:
        existing = truth[0] == second
        allele = np.where(existing ^ wrong, second, first)
        calls = np.stack((allele, truth[1]))
    else:
        dosage = np.count_nonzero(truth == second, axis=0)
        dosage = (dosage + wrong * error_rng.integers(1, 3, size=n)) % 3
        calls = np.sort(
            np.stack(
                (
                    np.where(dosage == 2, second, first),
                    np.where(dosage == 0, first, second),
                )
            ),
            axis=0,
        )
    calls[: 1 if haploid else 2, ~assayed | missing] = UNKNOWN
    return _ObservedCalls(_readonly(calls), _readonly(assayed))


def simulate_pedigree(
    pedigree: Pedigree,
    panel: Mapping[str, Panel],
    genetic_map: GeneticMap,
    *,
    tested: Collection[str] | None = None,
    observation: Observation | Mapping[str, Observation] = _DEFAULT_OBSERVATION,
    founders: _Model = "independent",
    seed: int | None = None,
) -> _SimulatedArrays:
    """Generate pedigree truth and unphased observations from a synthetic panel.

    Parameters
    ----------
    pedigree : Pedigree
        People and parent links to simulate.
    panel : mapping of str to Panel
        Synthetic markers on each simulated chromosome.
    genetic_map : GeneticMap
        Map for the simulated chromosomes.
    tested : collection of str, optional
        People whose test results to simulate. Defaults to everyone.
    observation : Observation or mapping of str to Observation, optional
        Shared or per-person settings. Omitted settings use the default rates.
    founders : {"independent", "correlated"}, default "independent"
        Draw founding copies per locus or from shared templates in stretches.
    seed : int, optional
        Nonnegative seed. The result records the seed used.
    """
    if not isinstance(pedigree, Pedigree):
        raise TypeError("pedigree must be a Pedigree")
    if not isinstance(genetic_map, GeneticMap):
        raise TypeError("genetic_map must be a GeneticMap")
    if not isinstance(founders, str):
        raise TypeError("founders must be a str")
    if founders not in ("independent", "correlated"):
        raise ValueError("founders must be 'independent' or 'correlated'")

    selected = _select_tested(pedigree, tested)
    settings = _settings(pedigree, observation)
    panel = _panels(panel, genetic_map)
    run_seed = _seed(seed)

    truth: dict[str, dict[str, Haplotypes]] = {name: {} for name in pedigree.names}
    transmissions: dict[tuple[str, str], dict[str, _Bytes]] = {
        (parent, person.name): {}
        for person in pedigree
        for parent in (person.mother, person.father)
        if parent is not None
    }
    for region, block in panel.items():
        copies, inherited = _inherit(
            pedigree, region, block, genetic_map, run_seed, founders
        )
        for name, value in copies.items():
            truth[name][region] = value
        for edge, value in inherited.items():
            transmissions[edge][region] = value

    observed = {
        name: MappingProxyType(
            {
                region: _observe(
                    truth[name][region],
                    block,
                    settings[name],
                    run_seed,
                    name,
                    haploid=region == "X" and pedigree[name].sex == Sex.MALE,
                )
                for region, block in panel.items()
            }
        )
        for name in selected
    }
    return _SimulatedArrays(
        pedigree,
        panel,
        genetic_map,
        MappingProxyType(
            {name: MappingProxyType(blocks) for name, blocks in truth.items()}
        ),
        MappingProxyType(
            {edge: MappingProxyType(blocks) for edge, blocks in transmissions.items()}
        ),
        MappingProxyType(observed),
        MappingProxyType(
            {
                "seed": run_seed,
                "model_version": "1",
                "founders": founders,
                "observation": settings,
                "tested": selected,
                "map_name": genetic_map.name,
                "map_sha256": genetic_map.sha256,
            }
        ),
    )
