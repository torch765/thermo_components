from pathlib import Path

import density
from thermo_components.adapters.thermo import ThermoGateway
from thermo_components.adapters.ui.qt_main_window import MainWindow
from thermo_components.bootstrap import resource_path


def test_density_launcher_exports_intentional_legacy_aliases(monkeypatch, tmp_path):
    expected_exports = {
        "MainWindow",
        "MixtureCalculator",
        "load_lhv_data",
        "main",
        "resource_path",
    }

    assert set(density.__all__) == expected_exports
    assert density.MainWindow is MainWindow
    assert density.MixtureCalculator is ThermoGateway

    monkeypatch.chdir(tmp_path)

    assert Path(density.resource_path("lhv_data.db")) == Path(
        resource_path("lhv_data.db")
    )
