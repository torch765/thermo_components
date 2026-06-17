"""Compatibility launcher for the Thermo Components desktop app."""

import sys
from pathlib import Path

from PyQt6.QtWidgets import QApplication, QMessageBox, QStyleFactory

from thermo_components.adapters.thermo import ThermoGateway
from thermo_components.adapters.ui.qt_main_window import MainWindow
from thermo_components.bootstrap import (
    build_desktop_dependencies,
    load_lhv_data as bootstrap_load_lhv_data,
    resource_path as bootstrap_resource_path,
)


MixtureCalculator = ThermoGateway

__all__ = [
    "MainWindow",
    "MixtureCalculator",
    "QMessageBox",
    "load_lhv_data",
    "main",
    "resource_path",
]


def load_lhv_data(db_path="lhv_data.db"):
    """Load LHV data through the SQLite persistence adapter."""
    return bootstrap_load_lhv_data(db_path)


def resource_path(relative_path):
    """Resolve a resource in source mode or a PyInstaller bundle."""
    return bootstrap_resource_path(relative_path)


def main():
    app = QApplication(sys.argv)
    app.setStyle(QStyleFactory.create("Windows"))
    app.setPalette(app.style().standardPalette())

    dependencies = build_desktop_dependencies(
        source_base_dir=Path(__file__).resolve().parent,
    )
    window = MainWindow(dependencies=dependencies)
    window.show()

    sys.exit(app.exec())


if __name__ == "__main__":
    main()
