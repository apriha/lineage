"""Exceptions raised by lineage."""

from __future__ import annotations

__all__ = ["LineageError", "PedigreeError"]


class LineageError(Exception):
    """Base class of the exceptions lineage raises."""


class PedigreeError(LineageError, ValueError):
    """A pedigree is invalid."""
