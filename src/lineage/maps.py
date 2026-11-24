"""Genetic maps: :class:`GeneticMap`."""

from __future__ import annotations

import math
import re
from collections.abc import Iterator, Mapping, Sequence
from dataclasses import dataclass
from typing import Literal, Self

import numpy as np
import numpy.typing as npt
import pandas as pd

from lineage import _grch37

__all__ = ["GeneticMap"]

_COLUMNS = ("pos", "rate", "map")

_REGION_NAMES = '"1" to "22", "X_par1", "X" or "X_par2"'

_SHA256 = re.compile(r"[0-9A-Fa-f]{64}")


@dataclass(frozen=True, eq=False, slots=True)
class GeneticMap:
    """A GRCh37 genetic map, with a table per region.

    The regions, in natural order, are "1" to "22", "X_par1", "X" and "X_par2".
    "X_par1" and "X_par2" are the pseudoautosomal regions; "X" is the rest of X.

    Parameters
    ----------
    name : str
        A label, such as "HapMap2" or a 1000 Genomes population code.
    build : {37}
        The genome build.
    tables : mapping of str to pandas.DataFrame
        A table per region, with columns "pos" (bases), "rate" (cM/Mb) and "map"
        (cM). "pos" must be whole numbers inside the region and strictly increasing,
        "map" non-decreasing and "rate" non-negative. Other columns and the index are
        ignored.
    source : str, optional
        Where the map came from.
    sha256 : str, optional
        The source file's SHA-256 in hexadecimal; stored in lower case.

    Attributes
    ----------
    tables : mapping of str to pandas.DataFrame
        Read-only, in natural order. Each lookup returns a new copy with a RangeIndex
        and columns "pos" (uint32), "rate" and "map" (float64).

    Raises
    ------
    TypeError
        If an argument, or a value in `tables`, has the wrong type.
    ValueError
        If `name` or `tables` is empty, `build` is not 37, `sha256` is not valid, or
        `tables` names an unknown region or holds an invalid table.
    """

    name: str
    build: Literal[37]
    tables: Mapping[str, pd.DataFrame]
    source: str | None = None
    sha256: str | None = None

    def __post_init__(self) -> None:
        _check_name(self.name)
        _check_build(self.build)
        _check_source(self.source)
        object.__setattr__(self, "sha256", _check_sha256(self.sha256))
        object.__setattr__(self, "tables", _Tables(_check_tables(self.tables)))

    @classmethod
    def from_frames(
        cls,
        name: str,
        tables: Mapping[str, pd.DataFrame],
        *,
        source: str | None = None,
        sha256: str | None = None,
    ) -> Self:
        """Return a GRCh37 genetic map from a table per region.

        Parameters
        ----------
        name, tables, source, sha256
            As for :class:`GeneticMap`.

        Returns
        -------
        GeneticMap

        Raises
        ------
        TypeError, ValueError
            As for :class:`GeneticMap`.
        """
        return cls(name, 37, tables, source=source, sha256=sha256)

    @classmethod
    def uniform(
        cls, *, rate: float = 1.0, chromosomes: Sequence[str] | None = None
    ) -> Self:
        """Return a map with the same recombination rate everywhere.

        Each region's table has two rows, at the region's ends, and starts at 0 cM.
        The ends of "X" are the bases next to the pseudoautosomal regions.

        Parameters
        ----------
        rate : float, default 1.0
            The recombination rate in cM/Mb.
        chromosomes : sequence of str, optional
            The regions to include, in any order. By default, every region.

        Returns
        -------
        GeneticMap
            Named "uniform", with no source or sha256.

        Raises
        ------
        TypeError
            If `rate` is not an int or float, or `chromosomes` is a str or not a
            sequence.
        ValueError
            If `rate` is negative or not finite, or `chromosomes` is empty, repeats a
            region or names an unknown one.
        """
        value = _check_rate(rate)
        tables: dict[str, pd.DataFrame] = {}
        for region in _check_regions(chromosomes):
            first, last = _uniform_span(region)
            end = value * ((last - first) / 1e6)
            tables[region] = pd.DataFrame(
                {"pos": [first, last], "rate": [value, value], "map": [0.0, end]}
            )
        return cls.from_frames("uniform", tables)

    @property
    def regions(self) -> tuple[str, ...]:
        """The regions the map covers, in natural order."""
        return tuple(self.tables)

    def cm(self, region: str, positions: npt.ArrayLike) -> npt.NDArray[np.float64]:
        """Return the genetic positions of base-pair positions in a region.

        Interpolates linearly between the table's rows. Positions outside the table
        get the value at its nearest end.

        Parameters
        ----------
        region : str
            A region the map covers.
        positions : array_like of int or float
            Positions in bases.

        Returns
        -------
        numpy.ndarray of float64
            Genetic positions in cM, with the shape of `positions`.

        Raises
        ------
        KeyError
            If the map does not cover `region`.
        ValueError
            If a position is not a finite number.
        """
        if region not in self.tables:
            listed = ", ".join(repr(name) for name in self.tables)
            raise KeyError(
                f"The map has no region {region!r}; use one of its regions ({listed})"
            )
        table = _stored(self.tables)[region]
        x = _float_positions(positions)
        return np.asarray(
            np.interp(
                x,
                table["pos"].to_numpy(dtype=np.float64),
                table["map"].to_numpy(dtype=np.float64),
            ),
            dtype=np.float64,
        )


class _Tables(Mapping[str, pd.DataFrame]):
    """Read-only tables of a map, compared by identity; each lookup returns a copy."""

    __slots__ = ("_frames",)

    __eq__ = object.__eq__
    __hash__ = object.__hash__

    def __init__(self, frames: dict[str, pd.DataFrame]) -> None:
        self._frames = frames

    def __getitem__(self, region: str) -> pd.DataFrame:
        stored = self._frames[region]
        return pd.DataFrame(
            {name: column.to_numpy(copy=True) for name, column in stored.items()}
        )

    def __contains__(self, region: object) -> bool:
        return isinstance(region, str) and region in self._frames

    def __iter__(self) -> Iterator[str]:
        return iter(self._frames)

    def __len__(self) -> int:
        return len(self._frames)

    def __repr__(self) -> str:
        rows = ", ".join(f"{key!r}: {len(t)} rows" for key, t in self._frames.items())
        return f"<tables {{{rows}}}>"


def _stored(tables: Mapping[str, pd.DataFrame]) -> Mapping[str, pd.DataFrame]:
    """Return a map's stored tables without copying them, or `tables` itself."""
    return tables._frames if isinstance(tables, _Tables) else tables


def _check_name(value: object) -> None:
    """Check that a genetic map's name is a non-empty str."""
    if not isinstance(value, str):
        raise TypeError(
            f"A genetic map's name must be a str, not {type(value).__name__}"
        )
    if not value:
        raise ValueError("A genetic map's name is empty; give a name")


def _check_build(value: object) -> None:
    """Check that a genome build is the int 37."""
    if type(value) is not int or value != 37:
        raise ValueError(f"The build of a genetic map must be 37, not {value!r}")


def _check_source(value: object) -> None:
    """Check that a genetic map's source is a str or None."""
    if value is not None and not isinstance(value, str):
        raise TypeError(f"source must be a str or None, not {type(value).__name__}")


def _check_sha256(value: object) -> str | None:
    """Return a SHA-256 of 64 hexadecimal characters in lower case, or None."""
    if value is None:
        return None
    if not isinstance(value, str):
        raise TypeError(f"sha256 must be a str or None, not {type(value).__name__}")
    if _SHA256.fullmatch(value) is None:
        raise ValueError(f"sha256 must be 64 hexadecimal characters, not {value!r}")
    return value.lower()


def _check_tables(tables: object) -> dict[str, pd.DataFrame]:
    """Return new tables of a genetic map, in natural order of region."""
    if not isinstance(tables, Mapping):
        raise TypeError(
            f"tables must be a mapping of region to DataFrame, not "
            f"{type(tables).__name__}"
        )
    if not tables:
        raise ValueError("tables is empty; give a table for at least one region")
    checked = {}
    for region, frame in _stored(tables).items():
        if region not in _grch37.REGIONS:
            raise ValueError(
                f"{region!r} is not a region of a genetic map; use {_REGION_NAMES}"
            )
        if not isinstance(frame, pd.DataFrame):
            raise TypeError(
                f"The table of region {region!r} must be a DataFrame, not "
                f"{type(frame).__name__}"
            )
        checked[region] = _check_table(frame, region)
    return {region: checked[region] for region in _grch37.REGIONS if region in checked}


def _check_table(frame: pd.DataFrame, region: str) -> pd.DataFrame:
    """Return a new table of a region, with only the columns of a genetic map."""
    label = f"The table of region {region!r}"
    names = frame.columns.tolist()
    missing = [name for name in _COLUMNS if name not in names]
    if missing:
        raise ValueError(
            f"{label} has no {' or '.join(map(repr, missing))} column; it needs "
            f"columns {', '.join(map(repr, _COLUMNS))}"
        )
    repeated = [name for name in _COLUMNS if names.count(name) > 1]
    if repeated:
        raise ValueError(f"{label} has more than one {repeated[0]!r} column")
    if len(frame) == 0:
        raise ValueError(f"{label} is empty; give at least one row")
    pos = _grch37.check_positions(
        frame["pos"], region, f"Positions in the table of region {region!r}"
    )
    rate, cm = (
        _check_column(frame[name], f"Column {name!r} of the table of region {region!r}")
        for name in ("rate", "map")
    )
    _check_rows(
        frame["map"],
        np.r_[False, cm[1:] < cm[:-1]],
        f"{label} has genetic positions that decrease",
    )
    _check_rows(frame["rate"], rate < 0, f"{label} has negative rates")
    return pd.DataFrame({"pos": pos, "rate": rate, "map": cm})


def _check_column(column: pd.Series, label: str) -> npt.NDArray[np.float64]:
    """Return a column of finite numbers as a new float64 array."""
    dtype = column.dtype
    if (
        pd.api.types.is_bool_dtype(dtype)
        or pd.api.types.is_complex_dtype(dtype)
        or not pd.api.types.is_numeric_dtype(dtype)
    ):
        raise ValueError(f"{label} must be numeric, not {dtype}")
    if column.isna().any():
        raise ValueError(f"{label} has missing values")
    values = column.to_numpy(dtype=np.float64, copy=True)
    _check_rows(column, ~np.isfinite(values), f"{label} has values that are not finite")
    return values


def _check_rows(column: pd.Series, bad: npt.NDArray[np.bool_], message: str) -> None:
    """Raise a ValueError naming the first bad row, if there is one."""
    if bad.any():
        i = int(np.argmax(bad))
        raise ValueError(f"{message}; row {i} has {column.iloc[i]}")


def _check_rate(rate: object) -> float:
    """Return a finite, non-negative recombination rate as a float."""
    if isinstance(rate, bool) or not isinstance(rate, int | float):
        raise TypeError(f"rate must be an int or float, not {type(rate).__name__}")
    if not math.isfinite(rate) or rate < 0:
        raise ValueError(f"rate must be finite and at least 0, not {rate!r}")
    return float(rate)


def _check_regions(chromosomes: object) -> list[str]:
    """Return the regions of a uniform map, in natural order."""
    if chromosomes is None:
        return list(_grch37.REGIONS)
    if isinstance(chromosomes, str) or not isinstance(chromosomes, Sequence):
        raise TypeError(
            f"chromosomes must be a sequence of region names, not "
            f"{type(chromosomes).__name__}"
        )
    if not chromosomes:
        raise ValueError("chromosomes is empty; give at least one region")
    for region in chromosomes:
        if region not in _grch37.REGIONS:
            raise ValueError(
                f"{region!r} is not a region of a genetic map; use {_REGION_NAMES}"
            )
    if len(set(chromosomes)) != len(chromosomes):
        raise ValueError(f"chromosomes repeats a region: {list(chromosomes)!r}")
    return [region for region in _grch37.REGIONS if region in chromosomes]


def _uniform_span(region: str) -> tuple[int, int]:
    """Return the first and last base of a region in a uniform map."""
    if region == "X":
        return (_grch37.PAR1[1] + 1, _grch37.PAR2[0] - 1)
    return _grch37.BOUNDS[region]


def _float_positions(positions: npt.ArrayLike) -> npt.NDArray[np.float64]:
    """Return finite numeric positions as a new float64 array."""
    if _grch37.holds_bool(positions):
        raise ValueError("Positions must be numbers, not bool")
    values = np.asarray(positions)
    if values.dtype.kind not in "iuf":
        raise ValueError(f"Positions must be numbers, not {values.dtype}")
    x = values.astype(np.float64)
    if not np.isfinite(x).all():
        raise ValueError("Positions must be finite numbers")
    return x
