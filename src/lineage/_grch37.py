"""Coordinates of the GRCh37 reference genome, and checks of positions against them.

Positions are 1-based and intervals include both ends.

References
----------
.. [1] Genome Reference Consortium, "Genome Reference Consortium Human Build 37
   (GRCh37)", 2009, RefSeq assembly GCF_000001405.13 (GenBank GCA_000001405.1),
   assembly report
   https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.13_GRCh37/GCF_000001405.13_GRCh37_assembly_report.txt
   (chromosome lengths).
.. [2] Genome Reference Consortium, "Human Genome Overview",
   https://www.ncbi.nlm.nih.gov/grc/human (GRCh37 pseudoautosomal regions).
"""

from __future__ import annotations

from collections.abc import Mapping
from types import MappingProxyType
from typing import Any, Final, TypeAlias

import numpy as np
import numpy.typing as npt

_Number: TypeAlias = np.integer[Any] | np.floating[Any]

# Neither type can be subclassed, so comparing exact types finds every bool.
_BOOLS: Final = frozenset((bool, np.bool_))

#: Lengths in bases of chromosomes 1 to 22 and X [1]_.
LENGTHS: Final[Mapping[str, int]] = MappingProxyType(
    {
        "1": 249_250_621,
        "2": 243_199_373,
        "3": 198_022_430,
        "4": 191_154_276,
        "5": 180_915_260,
        "6": 171_115_067,
        "7": 159_138_663,
        "8": 146_364_022,
        "9": 141_213_431,
        "10": 135_534_747,
        "11": 135_006_516,
        "12": 133_851_895,
        "13": 115_169_878,
        "14": 107_349_540,
        "15": 102_531_392,
        "16": 90_354_753,
        "17": 81_195_210,
        "18": 78_077_248,
        "19": 59_128_983,
        "20": 63_025_520,
        "21": 48_129_895,
        "22": 51_304_566,
        "X": 155_270_560,
    }
)

#: First and last base of the first pseudoautosomal region of X [2]_.
PAR1: Final = (60_001, 2_699_520)

#: First and last base of the second pseudoautosomal region of X [2]_.
PAR2: Final = (154_931_044, 155_260_560)

#: Regions of a genetic map in natural order: the autosomes, then X split at the
#: pseudoautosomal regions.
REGIONS: Final[tuple[str, ...]] = (
    *(str(i) for i in range(1, 23)),
    "X_par1",
    "X",
    "X_par2",
)

#: First and last base of each region. For "X" these are the ends of the chromosome;
#: :func:`check_positions` also excludes the pseudoautosomal regions.
BOUNDS: Final[Mapping[str, tuple[int, int]]] = MappingProxyType(
    {
        **{chrom: (1, length) for chrom, length in LENGTHS.items()},
        "X_par1": PAR1,
        "X_par2": PAR2,
    }
)


def check_positions(value: object, region: str, label: str) -> npt.NDArray[np.uint32]:
    """Return positions in `region` as a new read-only uint32 array."""
    if holds_bool(value):
        raise ValueError(f"{label} must be integers, not bool")
    pos = np.asarray(value)
    if pos.ndim != 1:
        raise ValueError(f"{label} must be a 1-D array, not {pos.ndim}-D")
    if pos.dtype.kind not in "iuf":
        raise ValueError(f"{label} must be integers, not {pos.dtype}")
    if pos.dtype.kind == "f":
        # Widen first: compared with a float32 or float16 array, a bound is rounded
        # to the array's precision.
        pos = pos.astype(np.promote_types(pos.dtype, np.float64))
        bad = ~np.isfinite(pos) | (np.trunc(pos) != pos)
        if bad.any():
            i = int(np.argmax(bad))
            raise ValueError(
                f"{label} must be finite whole numbers; {pos[i]} at index {i} is not"
            )
    outside = ~_in_region(region, pos)
    if outside.any():
        i = int(np.argmax(outside))
        raise ValueError(
            f"{label} must be {_describe(region)}; {pos[i]} at index {i} is not"
        )
    repeated = pos[1:] <= pos[:-1]
    if repeated.any():
        i = int(np.argmax(repeated)) + 1
        raise ValueError(
            f"{label} must be strictly increasing; {pos[i]} at index {i} does not "
            f"follow {pos[i - 1]}"
        )
    checked = pos.astype(np.uint32)
    checked.setflags(write=False)
    return checked


def holds_bool(value: object) -> bool:
    """Return whether plain Python input, such as a list or tuple, holds a bool."""
    if hasattr(value, "dtype"):
        return False
    items = np.asarray(value, dtype=object)
    return not _BOOLS.isdisjoint(map(type, items.flat))


def _in_region(region: str, pos: npt.NDArray[_Number]) -> npt.NDArray[np.bool_]:
    """Return which positions lie in a region."""
    first, last = BOUNDS[region]
    inside = (pos >= first) & (pos <= last)
    if region == "X":
        in_par1 = (pos >= PAR1[0]) & (pos <= PAR1[1])
        in_par2 = (pos >= PAR2[0]) & (pos <= PAR2[1])
        inside &= ~(in_par1 | in_par2)
    return inside


def _describe(region: str) -> str:
    """Describe the positions of a region, for an error message."""
    first, last = BOUNDS[region]
    text = f"{first:,} to {last:,}"
    if region == "X":
        text += (
            f", outside the pseudoautosomal regions {PAR1[0]:,} to {PAR1[1]:,} and "
            f"{PAR2[0]:,} to {PAR2[1]:,}"
        )
    return text
