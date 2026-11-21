"""Types shared across lineage."""

from __future__ import annotations

from enum import StrEnum

__all__ = ["Sex"]


class Sex(StrEnum):
    """A person's sex."""

    MALE = "Male"
    FEMALE = "Female"
