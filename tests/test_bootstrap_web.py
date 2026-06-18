from thermo_components.adapters.reporting import OpenPyxlReportExporter
from thermo_components.application.services import CalculationSessionService
from thermo_components.application.use_cases import (
    ConvertFlowUseCase,
    DeriveCompositionUseCase,
    NormalizeCompositionUseCase,
)
from thermo_components.bootstrap.web import build_web_dependencies


def test_build_web_dependencies_uses_supplied_lhv_data():
    dependencies = build_web_dependencies(
        lhv_data={"methane": 35.8},
    )

    assert dependencies.lhv_database == {"methane": 35.8}
    assert isinstance(
        dependencies.derive_composition_use_case,
        DeriveCompositionUseCase,
    )
    assert isinstance(
        dependencies.normalize_composition_use_case,
        NormalizeCompositionUseCase,
    )
    assert isinstance(
        dependencies.convert_flow_use_case,
        ConvertFlowUseCase,
    )
    assert isinstance(
        dependencies.report_exporter,
        OpenPyxlReportExporter,
    )


def test_web_dependencies_create_request_scoped_calculation_sessions():
    dependencies = build_web_dependencies(lhv_data={})

    first = dependencies.calculation_session_factory("PRMIX")
    second = dependencies.calculation_session_factory("PRMIX")

    assert isinstance(first, CalculationSessionService)
    assert isinstance(second, CalculationSessionService)
    assert first is not second
    assert (
        first.calculate_properties_use_case._calculator
        is not second.calculate_properties_use_case._calculator
    )
