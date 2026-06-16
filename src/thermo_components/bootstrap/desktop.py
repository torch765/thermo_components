"""Composition root for the PyQt desktop application."""

from collections.abc import Callable, Mapping
from dataclasses import dataclass
from pathlib import Path

from thermo_components.adapters.packaging import RuntimeResourceLocator
from thermo_components.adapters.persistence import SqliteLhvRepository
from thermo_components.adapters.reporting import OpenPyxlReportExporter
from thermo_components.adapters.thermo import ThermoGateway
from thermo_components.adapters.ui import QtReportExportController
from thermo_components.application.use_cases import (
    CalculatePropertiesUseCase,
    ConvertFlowUseCase,
    DeriveCompositionUseCase,
    NormalizeCompositionUseCase,
    PrepareReportUseCase,
)


@dataclass(frozen=True)
class DesktopDependencies:
    """Concrete services used by the desktop Qt adapter."""

    lhv_database: Mapping[str, float]
    calculator: ThermoGateway
    calculate_properties_use_case: CalculatePropertiesUseCase
    convert_flow_use_case: ConvertFlowUseCase
    derive_composition_use_case: DeriveCompositionUseCase
    normalize_composition_use_case: NormalizeCompositionUseCase
    prepare_report_use_case: PrepareReportUseCase
    report_exporter: OpenPyxlReportExporter
    report_export_controller_factory: Callable[
        [object],
        QtReportExportController,
    ]


def build_desktop_dependencies(
    *,
    lhv_data: Mapping[str, float] | None = None,
    lhv_db_path: str | Path | None = None,
    source_base_dir: str | Path | None = None,
) -> DesktopDependencies:
    """Build concrete adapters and application services for the desktop UI."""
    source_root = (
        Path(source_base_dir)
        if source_base_dir is not None
        else Path.cwd()
    )
    loaded_lhv_data = (
        dict(lhv_data)
        if lhv_data is not None
        else load_lhv_data(
            lhv_db_path
            if lhv_db_path is not None
            else RuntimeResourceLocator(source_root=source_root).resolve(
                "lhv_data.db"
            )
        )
    )

    calculator = ThermoGateway()
    report_exporter = OpenPyxlReportExporter()

    def build_report_controller(parent) -> QtReportExportController:
        return QtReportExportController(
            parent,
            report_exporter,
            source_base_dir=source_root,
        )

    return DesktopDependencies(
        lhv_database=loaded_lhv_data,
        calculator=calculator,
        calculate_properties_use_case=CalculatePropertiesUseCase(
            calculator,
            loaded_lhv_data,
        ),
        convert_flow_use_case=ConvertFlowUseCase(),
        derive_composition_use_case=DeriveCompositionUseCase(),
        normalize_composition_use_case=NormalizeCompositionUseCase(),
        prepare_report_use_case=PrepareReportUseCase(),
        report_exporter=report_exporter,
        report_export_controller_factory=build_report_controller,
    )


def load_lhv_data(db_path: str | Path = "lhv_data.db") -> dict[str, float]:
    """Load LHV data through the SQLite persistence adapter."""
    return dict(SqliteLhvRepository(db_path).load_all())


def resource_path(relative_path: str | Path) -> str:
    """Resolve a resource in source mode or a PyInstaller bundle."""
    return str(RuntimeResourceLocator().resolve(relative_path))
