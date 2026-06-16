"""Ports for persisted application reference data."""

from collections.abc import Mapping
from typing import Protocol, runtime_checkable


@runtime_checkable
class LhvRepository(Protocol):
    """Load lower-heating-value reference data by component identity."""

    def load_all(self) -> Mapping[str, float]: ...
