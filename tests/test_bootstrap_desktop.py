from thermo_components.adapters.persistence import SqliteLhvRepository
from thermo_components.adapters.thermo import ThermoGateway
from thermo_components.adapters.ui import QtReportExportController
from thermo_components.application.use_cases import (
    CalculatePropertiesUseCase,
    ConvertFlowUseCase,
    DeriveCompositionUseCase,
    NormalizeCompositionUseCase,
    PrepareReportUseCase,
)
from thermo_components.bootstrap import build_desktop_dependencies


def test_build_desktop_dependencies_wires_services_from_mapping(tmp_path):
    dependencies = build_desktop_dependencies(
        lhv_data={"methane": 35.8},
        source_base_dir=tmp_path,
    )

    assert dependencies.lhv_database == {"methane": 35.8}
    assert isinstance(dependencies.calculator, ThermoGateway)
    assert isinstance(
        dependencies.calculate_properties_use_case,
        CalculatePropertiesUseCase,
    )
    assert isinstance(dependencies.convert_flow_use_case, ConvertFlowUseCase)
    assert isinstance(
        dependencies.derive_composition_use_case,
        DeriveCompositionUseCase,
    )
    assert isinstance(
        dependencies.normalize_composition_use_case,
        NormalizeCompositionUseCase,
    )
    assert isinstance(dependencies.prepare_report_use_case, PrepareReportUseCase)

    report_controller = dependencies.report_export_controller_factory(object())

    assert isinstance(report_controller, QtReportExportController)
    assert report_controller.report_exporter is dependencies.report_exporter
    assert report_controller.get_export_base_dir() == tmp_path


def test_build_desktop_dependencies_loads_lhv_database(tmp_path):
    database_path = tmp_path / "lhv_data.db"
    SqliteLhvRepository(database_path).upsert_all({"methane": 35.8})

    dependencies = build_desktop_dependencies(
        lhv_db_path=database_path,
        source_base_dir=tmp_path,
    )

    assert dependencies.lhv_database == {"methane": 35.8}
    assert isinstance(
        dependencies.calculate_properties_use_case,
        CalculatePropertiesUseCase,
    )
