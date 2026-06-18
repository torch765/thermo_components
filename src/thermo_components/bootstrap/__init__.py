"""Application composition and startup."""

from importlib import import_module

__all__ = [
    "DesktopDependencies",
    "build_desktop_dependencies",
    "load_lhv_data",
    "resource_path",
]

_DESKTOP_EXPORTS = frozenset(__all__)


def __getattr__(name: str):
    """Load desktop composition only when a desktop export is requested."""

    if name not in _DESKTOP_EXPORTS:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    return getattr(import_module(f"{__name__}.desktop"), name)
