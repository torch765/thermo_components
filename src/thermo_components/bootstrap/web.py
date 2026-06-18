"""Composition root for the FastAPI web application."""

from collections.abc import Callable, Mapping
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path

from thermo_components.adapters.packaging import RuntimeResourceLocator
from thermo_components.adapters.persistence import SqliteLhvRepository
from thermo_components.adapters.reporting import OpenPyxlReportExporter
from thermo_components.adapters.thermo import ThermoGateway
from thermo_components.application.ports import ReportExporter
from thermo_components.application.services import CalculationSessionService
from thermo_components.application.use_cases import (
    CalculatePropertiesUseCase,
    ConvertFlowUseCase,
    DeriveCompositionUseCase,
    NormalizeCompositionUseCase,
    PrepareReportUseCase,
)


@dataclass(frozen=True)
class WebDependencies:
    """Concrete services and factories used by the FastAPI adapter."""

    lhv_database: Mapping[str, float]
    derive_composition_use_case: DeriveCompositionUseCase
    calculation_session_factory: Callable[[str], CalculationSessionService]
    normalize_composition_use_case: NormalizeCompositionUseCase = field(
        default_factory=NormalizeCompositionUseCase
    )
    convert_flow_use_case: ConvertFlowUseCase = field(
        default_factory=ConvertFlowUseCase
    )
    report_exporter: ReportExporter = field(
        default_factory=OpenPyxlReportExporter
    )
    report_clock: Callable[[], datetime] = datetime.now


def build_web_dependencies(
    *,
    lhv_data: Mapping[str, float] | None = None,
    lhv_db_path: str | Path | None = None,
    source_base_dir: str | Path | None = None,
) -> WebDependencies:
    """Build web dependencies without importing the desktop Qt adapter."""

    source_root = (
        Path(source_base_dir)
        if source_base_dir is not None
        else Path.cwd()
    )
    loaded_lhv_data = (
        dict(lhv_data)
        if lhv_data is not None
        else dict(
            SqliteLhvRepository(
                lhv_db_path
                if lhv_db_path is not None
                else RuntimeResourceLocator(source_root=source_root).resolve(
                    "lhv_data.db"
                )
            ).load_all()
        )
    )

    def build_calculation_session(model: str) -> CalculationSessionService:
        # ThermoGateway is stateful, so each HTTP request receives a new one.
        calculator = ThermoGateway(eos=model)
        return CalculationSessionService(
            CalculatePropertiesUseCase(calculator, loaded_lhv_data),
            PrepareReportUseCase(),
        )

    return WebDependencies(
        lhv_database=loaded_lhv_data,
        derive_composition_use_case=DeriveCompositionUseCase(),
        normalize_composition_use_case=NormalizeCompositionUseCase(),
        convert_flow_use_case=ConvertFlowUseCase(),
        report_exporter=OpenPyxlReportExporter(),
        calculation_session_factory=build_calculation_session,
    )
