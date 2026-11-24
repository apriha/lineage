"""Types shared across lineage."""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
from typing import Self, TypeAlias

import numpy as np
import numpy.typing as npt

from lineage import _grch37

__all__ = ["Loci", "Sex"]

_UInt32s: TypeAlias = npt.NDArray[np.uint32]
_Objects: TypeAlias = npt.NDArray[np.object_]


class Sex(StrEnum):
    """A person's sex."""

    MALE = "Male"
    FEMALE = "Female"


@dataclass(frozen=True, init=False, eq=False, slots=True)
class Loci:
    """The loci of one chromosome, in GRCh37 coordinates.

    The arrays are read-only copies of those given.

    Parameters
    ----------
    chrom : str
        The chromosome: "1" to "22" or "X".
    pos : array_like of int
        1-based positions within the chromosome, strictly increasing and, on X,
        outside the pseudoautosomal regions. Whole-number floats are accepted.
        Stored as uint32.
    rsid : array_like of str or None, optional
        Each locus's name, or None; several names are joined by ";". By default, no
        locus has a name.

    Raises
    ------
    TypeError
        If `chrom` is not a str.
    ValueError
        If `chrom`, `pos` or `rsid` is not valid.
    """

    chrom: str
    pos: npt.NDArray[np.uint32]
    rsid: npt.NDArray[np.object_]

    def __init__(
        self, chrom: str, pos: npt.ArrayLike, rsid: npt.ArrayLike | None = None
    ) -> None:
        chrom = _check_chrom(chrom)
        checked_pos = _grch37.check_positions(
            pos, chrom, f"Positions on chromosome {chrom!r}"
        )
        checked_rsid = _check_rsid(rsid, len(checked_pos))
        object.__setattr__(self, "chrom", chrom)
        object.__setattr__(self, "pos", checked_pos)
        object.__setattr__(self, "rsid", checked_rsid)

    def __reduce__(self) -> tuple[type[Self], tuple[str, _UInt32s, _Objects]]:
        return (type(self), (self.chrom, self.pos, self.rsid))


def _check_chrom(value: object) -> str:
    """Return `value` if it names a chromosome of :class:`Loci`."""
    if not isinstance(value, str):
        raise TypeError(f"A chromosome must be a str, not {type(value).__name__}")
    if value not in _grch37.LENGTHS:
        raise ValueError(f'A chromosome must be "1" to "22" or "X", not {value!r}')
    return value


def _check_rsid(value: object, n: int) -> npt.NDArray[np.object_]:
    """Return `n` locus names as a new read-only object array."""
    if value is None:
        rsid = np.full(n, None, dtype=object)
    else:
        rsid = np.array(value, dtype=object)
        if rsid.ndim != 1 or len(rsid) != n:
            raise ValueError(
                f"rsid must be a 1-D array of {n} names, one per position, not shape "
                f"{rsid.shape}"
            )
        valid = np.fromiter(map(_is_rsid, rsid), dtype=np.bool_, count=n)
        if not valid.all():
            i = int(np.argmin(valid))
            raise ValueError(
                f"Each rsid must be a non-empty str or None; {rsid[i]!r} at index {i} "
                "is not"
            )
    rsid.setflags(write=False)
    return rsid


def _is_rsid(value: object) -> bool:
    """Return whether `value` is a non-empty str or None."""
    return value is None or (isinstance(value, str) and value != "")
