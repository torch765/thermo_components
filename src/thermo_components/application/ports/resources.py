"""Port for locating runtime resources."""

from pathlib import Path
from typing import Protocol, runtime_checkable


@runtime_checkable
class ResourceLocator(Protocol):
    """Resolve a resource path in source and packaged execution modes."""

    def resolve(self, relative_path: str | Path) -> Path: ...
