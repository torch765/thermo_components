from datetime import datetime

import pytest

from density import MainWindow
from thermo_components.adapters.ui import FlowTabController
from thermo_components.bootstrap import build_desktop_dependencies
from thermo_components.domain.thermo_routes import PRMIX_DEFAULT_ROUTE


@pytest.fixture
def main_window(qt_app):
    window = MainWindow(lhv_data={"methane": 35.8})
    yield window
    window.close()


def _select_component(window, component_name):
    index = window.ui.comboBox_select_components.findText(component_name)
    assert index >= 0
    window.ui.comboBox_select_components.setCurrentIndex(index)


def _sample_result_mapping(**overrides):
    data = {
        "mw": 16.04,
        "comp_names": ["methane"],
        "mol_percents": [100.0],
        "wt_percents": [100.0],
        "thermo_route": PRMIX_DEFAULT_ROUTE,
        "model_display": "PRMIX",
        "density_result": 0.72,
        "phase": "Vapor",
        "density_error": None,
        "density_normal_kg_m3": 0.716,
        "density_normal_phase": "Vapor",
        "density_normal_error": None,
        "density_standard_kg_m3": 0.68,
        "density_standard_phase": "Vapor",
        "density_standard_error": None,
        "pressure_atm": 1.0,
        "bubble_point": -161.5,
        "bp_error": None,
        "mixture_lhv": 35.8,
        "missing_lhv": [],
        "basis": "Mol %",
        "eos": "PRMIX",
        "warnings": [],
    }
    data.update(overrides)
    return data


def test_main_window_accepts_bootstrap_dependencies(qt_app, tmp_path):
    dependencies = build_desktop_dependencies(
        lhv_data={"methane": 35.8},
        source_base_dir=tmp_path,
    )
    window = MainWindow(dependencies=dependencies)
    try:
        assert window.dependencies is dependencies
        assert window.lhv_database == {"methane": 35.8}
        assert window.calculator is dependencies.calculator
        assert (
            window.calculate_properties_use_case
            is dependencies.calculate_properties_use_case
        )
        assert isinstance(window.flow_tab_controller, FlowTabController)
        assert window.report_exporter is dependencies.report_exporter
        assert (
            window.report_export_controller.report_exporter
            is dependencies.report_exporter
        )
    finally:
        window.close()


def test_flow_tab_renders_mass_conversion(main_window):
    main_window.ui.lineEdit_enter_flow.setText("1")
    main_window.ui.comboBox_select_units.setCurrentText("kg/h")
    main_window.ui.comboBox_select_desired_units.setCurrentText("t/d")

    main_window.update_flow_conversion()

    assert main_window.ui.lineEdit_result.text() == "0.024"
    assert main_window.ui.lineEdit_in_desired_units.text() == "t/d"


def test_flow_tab_renders_missing_density_message(main_window):
    main_window.ui.lineEdit_enter_flow.setText("1")
    main_window.ui.comboBox_select_units.setCurrentText("kg/h")
    main_window.ui.comboBox_select_desired_units.setCurrentText("Nm3/h")

    main_window.update_flow_conversion()

    assert main_window.ui.lineEdit_result.text() == "Normal density required"
    assert main_window.ui.lineEdit_in_desired_units.text() == "Nm3/h"


def test_report_helpers_capture_conditions_from_ui_and_latest_result(
    main_window,
):
    main_window.ui.radioButton_wt_percent.setChecked(True)

    ui_conditions = dict(main_window.get_current_conditions())

    assert ui_conditions["Basis"] == "Wt %"
    assert ui_conditions["Model / EOS"] == "PRMIX"

    main_window.last_result_data = _sample_result_mapping(
        basis="Mol %",
        model_display="IAPWS-95",
    )
    result_conditions = dict(main_window.get_current_conditions())

    assert result_conditions["Basis"] == "Mol %"
    assert result_conditions["Model / EOS"] == "IAPWS-95"


def test_report_helpers_capture_composition_rows(main_window):
    _select_component(main_window, "methane")
    table = main_window.ui.tableWidget
    table.item(0, 1).setText("12.5")
    table.item(0, 2).setText("not numeric")

    rows = main_window.get_input_composition_rows()

    assert rows[0] == {
        "Component": "methane",
        "Mol %": 12.5,
        "Wt %": "not numeric",
    }
    assert rows[-1]["Component"] == "Total"
    assert rows[-1]["Mol %"] == 12.5


def test_report_export_request_builder_uses_current_helpers(
    main_window,
    tmp_path,
):
    _select_component(main_window, "methane")
    main_window.ui.tableWidget.item(0, 1).setText("100")
    main_window.last_result_data = _sample_result_mapping()

    exported_at = datetime(2026, 6, 17, 9, 15, 0)
    request = main_window._build_report_export_request(
        tmp_path / "report.xlsx",
        exported_at,
    )

    assert request.report_path == tmp_path / "report.xlsx"
    assert request.exported_at == exported_at
    assert request.conditions[0].setting == "Basis"
    assert request.conditions[0].value == "Mol %"
    assert request.input_rows[0].component == "methane"
    assert request.input_rows[-1].component == "Total"
    assert request.projection.result_rows


def test_calculate_and_display_invalid_input_characterizes_progress(
    main_window,
    monkeypatch,
):
    progress_calls = []
    monkeypatch.setattr(
        main_window,
        "animate_progress_to",
        lambda target, duration_ms=500: progress_calls.append(
            (target, duration_ms)
        ),
    )

    main_window.calculate_and_display()

    assert main_window.ui.results_list.item(0).text() == (
        "Error: No components selected."
    )
    assert progress_calls == [(80, 1000), (100, 500)]
    assert not hasattr(main_window, "worker_thread")


def test_calculation_error_clears_stale_result_and_records_error(
    main_window,
    monkeypatch,
):
    progress_calls = []
    monkeypatch.setattr(
        main_window,
        "animate_progress_to",
        lambda target, duration_ms=500: progress_calls.append(
            (target, duration_ms)
        ),
    )
    main_window.last_result_data = _sample_result_mapping()
    main_window.density_actual_kg_m3 = 1.0
    main_window.density_normal_kg_m3 = 2.0
    main_window.density_standard_kg_m3 = 3.0
    main_window.ui.results_list.addItem("old result")

    main_window.on_calculation_error("calculation failed")

    assert main_window.last_result_data is None
    assert main_window.density_actual_kg_m3 is None
    assert main_window.density_normal_kg_m3 is None
    assert main_window.density_standard_kg_m3 is None
    assert main_window.ui.results_list.count() == 1
    assert main_window.ui.results_list.item(0).text() == "calculation failed"
    assert progress_calls == [(100, 500)]
