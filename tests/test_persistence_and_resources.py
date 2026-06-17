from pathlib import Path
import sys


PROJECT_ROOT = Path(__file__).resolve().parents[1]

from thermo_components.adapters.packaging import RuntimeResourceLocator
from thermo_components.adapters.persistence import SqliteLhvRepository
from thermo_components.application.ports import LhvRepository, ResourceLocator
from thermo_components.bootstrap import load_lhv_data, resource_path


def test_sqlite_lhv_repository_creates_and_loads_values(tmp_path):
    database_path = tmp_path / "lhv_data.db"
    repository = SqliteLhvRepository(database_path)

    assert isinstance(repository, LhvRepository)
    assert repository.upsert_all(
        {
            " Methane ": 35.8,
            "nitrogen": 0.0,
        }
    )
    assert repository.load_all() == {
        "methane": 35.8,
        "nitrogen": 0.0,
    }


def test_sqlite_lhv_repository_returns_empty_mapping_when_missing(
    tmp_path,
    capsys,
):
    repository = SqliteLhvRepository(tmp_path / "missing.db")

    assert repository.load_all() == {}
    assert "Database file not found" in capsys.readouterr().out


def test_runtime_resource_locator_uses_source_root(tmp_path):
    locator = RuntimeResourceLocator(source_root=tmp_path)

    assert isinstance(locator, ResourceLocator)
    assert locator.resolve("lhv_data.db") == (
        tmp_path / "lhv_data.db"
    ).resolve()


def test_runtime_resource_locator_prefers_bundle_root(tmp_path):
    source_root = tmp_path / "source"
    bundle_root = tmp_path / "bundle"
    locator = RuntimeResourceLocator(
        source_root=source_root,
        bundle_root=bundle_root,
    )

    assert locator.resolve("lhv_data.db") == (
        bundle_root / "lhv_data.db"
    ).resolve()


def test_runtime_resource_locator_detects_pyinstaller_root(
    monkeypatch,
    tmp_path,
):
    bundle_root = tmp_path / "bundle"
    monkeypatch.setattr(sys, "_MEIPASS", str(bundle_root), raising=False)

    locator = RuntimeResourceLocator(source_root=tmp_path / "source")

    assert locator.resolve("lhv_data.db") == (
        bundle_root / "lhv_data.db"
    ).resolve()


def test_lhv_loader_delegates_to_repository(tmp_path):
    database_path = tmp_path / "lhv_data.db"
    SqliteLhvRepository(database_path).upsert_all({"methane": 35.8})

    assert load_lhv_data(database_path) == {"methane": 35.8}


def test_resource_path_returns_absolute_string(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)

    resolved = resource_path("lhv_data.db")

    assert isinstance(resolved, str)
    assert Path(resolved) == (tmp_path / "lhv_data.db").resolve()


def test_pyinstaller_spec_uses_portable_package_collection():
    spec_text = (PROJECT_ROOT / "density.spec").read_text(
        encoding="utf-8"
    )
    forbidden_fragments = (
        "D:/",
        "D:\\",
        "OneDrive",
        ".venv/Lib/site-packages",
        ".venv\\Lib\\site-packages",
        "site-packages/chemicals",
        "site-packages\\chemicals",
        "site-packages/thermo",
        "site-packages\\thermo",
    )

    assert "collect_data_files" in spec_text
    assert "collect_submodules" in spec_text
    assert "collect_data_files('thermo')" in spec_text
    assert "collect_data_files('chemicals')" in spec_text
    assert "collect_submodules('thermo')" in spec_text
    assert "collect_submodules('chemicals')" in spec_text
    assert "('lhv_data.db', '.')" in spec_text
    assert "('gui.ui', '.')" in spec_text
    assert "hiddenimports=thermo_hiddenimports + chemicals_hiddenimports" in (
        spec_text
    )
    for fragment in forbidden_fragments:
        assert fragment not in spec_text
