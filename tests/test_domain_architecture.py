import ast
from pathlib import Path

import pytest


FORBIDDEN_IMPORT_ROOTS = {
    "PyQt6",
    "chemicals",
    "openpyxl",
    "sqlite3",
    "thermo",
}
PACKAGE_ROOT = (
    Path(__file__).resolve().parents[1] / "src" / "thermo_components"
)


def imported_root_names(path: Path) -> set[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    roots: set[str] = set()

    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            roots.update(alias.name.split(".", 1)[0] for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
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
