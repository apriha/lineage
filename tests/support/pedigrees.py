"""Invented pedigrees for reconstruction tests."""

from __future__ import annotations

from lineage.reconstruction.pedigree import Pedigree, Person
from lineage.types import Sex

_F = Sex.FEMALE
_M = Sex.MALE


def nuclear(n: int = 4) -> tuple[Pedigree, tuple[str, ...]]:
    """A mother, father and one to four children, with the mother and children tested."""
    if not 1 <= n <= 4:
        raise ValueError("n must be from 1 to 4")
    children = ("C1", "C2", "C3", "C4")[:n]
    pedigree = Pedigree(
        [
            Person("M1", _F),
            Person("F1", _M),
            *[
                Person(name, sex, "M1", "F1")
                for name, sex in zip(children, (_F, _M, _F, _M)[:n], strict=True)
            ],
        ]
    )
    return pedigree, ("M1", *children)


def collateral() -> tuple[Pedigree, tuple[str, ...]]:
    """An aunt and uncle tested alongside a father and two children."""
    pedigree = Pedigree(
        [
            Person("G1", _F),
            Person("G2", _M),
            Person("M1", _F, "G1", "G2"),
            Person("U1", _M, "G1", "G2"),
            Person("A1", _F, "G1", "G2"),
            Person("F1", _M),
            Person("C1", _F, "M1", "F1"),
            Person("C2", _M, "M1", "F1"),
        ]
    )
    return pedigree, ("U1", "A1", "F1", "C1", "C2")


def connected() -> tuple[Pedigree, tuple[str, ...]]:
    """Two disconnected families with twelve tested people."""
    pedigree = Pedigree(
        [
            Person("G1", _F),
            Person("G2", _M),
            Person("M2", _F, "G1", "G2"),
            Person("F1", _M, "G1", "G2"),
            Person("U1", _M, "G1", "G2"),
            Person("G3", _F),
            Person("G4", _M),
            Person("M1", _F, "G3", "G4"),
            Person("U2", _M, "G3", "G4"),
            Person("C1", _F, "M1", "F1"),
            Person("C2", _M, "M1", "F1"),
            Person("F2", _M),
            Person("C3", _F, "M2", "F2"),
            Person("G5", _F),
            Person("G6", _M),
            Person("M3", _F, "G5", "G6"),
            Person("G7", _F),
            Person("G8", _M),
            Person("F3", _M, "G7", "G8"),
            Person("C4", _M, "M3", "F3"),
            Person("C5", _F, "M3", "F3"),
        ]
    )
    return pedigree, (
        "M2",
        "U1",
        "G4",
        "M1",
        "U2",
        "C1",
        "C2",
        "F2",
        "C3",
        "M3",
        "C4",
        "C5",
    )


def maternal() -> tuple[Pedigree, tuple[str, ...]]:
    """A grandmother, mother, maternal uncle and four tested children."""
    children = ("Child1", "Child2", "Child3", "Child4")
    pedigree = Pedigree(
        [
            Person("Grandmother", _F),
            Person("Grandfather", _M),
            Person("Mother", _F, "Grandmother", "Grandfather"),
            Person("Uncle", _M, "Grandmother", "Grandfather"),
            Person("Father", _M),
            *[
                Person(name, sex, "Mother", "Father")
                for name, sex in zip(children, (_F, _M, _F, _M), strict=True)
            ],
        ]
    )
    return pedigree, ("Grandmother", "Mother", "Uncle", *children)
