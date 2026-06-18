import ast
from pathlib import Path
import subprocess
import sys

import pytest


FORBIDDEN_IMPORT_ROOTS = {
    "PyQt6",
    "chemicals",
    "fastapi",
    "openpyxl",
    "sqlite3",
    "thermo",
    "uvicorn",
}
PACKAGE_ROOT = (
    Path(__file__).resolve().parents[1] / "src" / "thermo_components"
)
PROJECT_ROOT = Path(__file__).resolve().parents[1]


def imported_root_names(path: Path) -> set[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    roots: set[str] = set()

    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            roots.update(alias.name.split(".", 1)[0] for alias in node.names)
        elif (
            isinstance(node, ast.ImportFrom)
            and node.level == 0
            and node.module
        ):
            roots.add(node.module.split(".", 1)[0])

    return roots


@pytest.mark.parametrize("layer_name", ["domain", "application"])
def test_inner_layers_have_no_framework_or_infrastructure_imports(layer_name):
    layer_root = PACKAGE_ROOT / layer_name
    violations = {
        str(path.relative_to(layer_root)): sorted(
            imported_root_names(path) & FORBIDDEN_IMPORT_ROOTS
        )
        for path in layer_root.rglob("*.py")
        if imported_root_names(path) & FORBIDDEN_IMPORT_ROOTS
    }
    assert violations == {}


def test_thermo_dependencies_are_confined_to_the_thermo_adapter():
    source_paths = [PROJECT_ROOT / "density.py", *PACKAGE_ROOT.rglob("*.py")]
    violations = {}

    for path in source_paths:
        thermo_imports = imported_root_names(path) & {"chemicals", "thermo"}
        if not thermo_imports:
            continue
        relative_path = path.relative_to(PROJECT_ROOT)
        if relative_path.parts[:4] != (
            "src",
            "thermo_components",
            "adapters",
            "thermo",
        ):
            violations[str(relative_path)] = sorted(thermo_imports)

    assert violations == {}


def test_sqlite_dependency_is_confined_to_the_persistence_adapter():
    source_paths = [
        PROJECT_ROOT / "density.py",
        PROJECT_ROOT / "lhv_data.py",
        *PACKAGE_ROOT.rglob("*.py"),
    ]
    violations = {}

    for path in source_paths:
        sqlite_imports = imported_root_names(path) & {"sqlite3"}
        if not sqlite_imports:
            continue
        relative_path = path.relative_to(PROJECT_ROOT)
        if relative_path.parts[:4] != (
            "src",
            "thermo_components",
            "adapters",
            "persistence",
        ):
            violations[str(relative_path)] = sorted(sqlite_imports)

    assert violations == {}


def test_openpyxl_dependency_is_confined_to_the_reporting_adapter():
    source_paths = [PROJECT_ROOT / "density.py", *PACKAGE_ROOT.rglob("*.py")]
    violations = {}

    for path in source_paths:
        openpyxl_imports = imported_root_names(path) & {"openpyxl"}
        if not openpyxl_imports:
            continue
        relative_path = path.relative_to(PROJECT_ROOT)
        if relative_path.parts[:4] != (
            "src",
            "thermo_components",
            "adapters",
            "reporting",
        ):
            violations[str(relative_path)] = sorted(openpyxl_imports)

    assert violations == {}


def test_web_framework_dependencies_are_confined_to_the_web_adapter():
    source_paths = [PROJECT_ROOT / "density.py", *PACKAGE_ROOT.rglob("*.py")]
    violations = {}

    for path in source_paths:
        web_imports = imported_root_names(path) & {"fastapi", "uvicorn"}
        if not web_imports:
            continue
        relative_path = path.relative_to(PROJECT_ROOT)
        if relative_path.parts[:4] != (
            "src",
            "thermo_components",
            "adapters",
            "web",
        ):
            violations[str(relative_path)] = sorted(web_imports)

    assert violations == {}


def test_importing_web_app_does_not_load_pyqt():
    script = """
import sys
import thermo_components.adapters.web.app

pyqt_loaded = any(
    name == "PyQt6" or name.startswith("PyQt6.")
    for name in sys.modules
)
raise SystemExit(1 if pyqt_loaded else 0)
"""

    result = subprocess.run(
        [sys.executable, "-c", script],
        cwd=PROJECT_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
