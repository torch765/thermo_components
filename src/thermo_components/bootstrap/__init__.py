"""Application composition and startup."""

from .desktop import (
    DesktopDependencies,
    build_desktop_dependencies,
    load_lhv_data,
    resource_path,
)

__all__ = [
    "DesktopDependencies",
    "build_desktop_dependencies",
    "load_lhv_data",
    "resource_path",
]
