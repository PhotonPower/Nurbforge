"""Helpers to apply transform patterns to geometry objects."""
from __future__ import annotations

from collections.abc import Iterable
from typing import TypeVar

from ._core import Transform

G = TypeVar("G")


def apply_pattern(geometry: G, transforms: Iterable[Transform]) -> list[G]:
    """Return one transformed copy of ``geometry`` per transform."""
    return [geometry.transformed(t) for t in transforms]  # type: ignore[attr-defined]
