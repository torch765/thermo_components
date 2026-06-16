"""Resolve resources in source checkouts and PyInstaller bundles."""

from pathlib import Path
import sys


class RuntimeResourceLocator:
    """Locate bundled inputs without exposing PyInstaller details."""

    def __init__(
        self,
        source_root: str | Path | None = None,
        bundle_root: str | Path | None = None,
    ):
        detected_bundle_root = (
            bundle_root
            if bundle_root is not None
            else getattr(sys, "_MEIPASS", None)
        )
        self._base_path = Path(
            detected_bundle_root
            if detected_bundle_root is not None
            else source_root if source_root is not None else Path.cwd()
        )

    def resolve(self, relative_path: str | Path) -> Path:
        path = Path(relative_path)
        if path.is_absolute():
            return path
        return (self._base_path / path).resolve()
