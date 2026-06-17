from datetime import datetime

import pytest

from density import MainWindow
from thermo_components.adapters.ui import QtReportRequestBuilder
from thermo_components.domain.thermo_routes import PRMIX_DEFAULT_ROUTE
from thermo_components.domain.warnings import PRMIX_WATER_WARNING


@pytest.fixture
def main_window(qt_app):
    window = MainWindow(lhv_data={"methane": 35.8})
    yield window
    window.close()


def _select_component(window, component_name: str) -> None:
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


def test_report_request_builder_reads_ui_conditions_until_result_exists(
    main_window,
):
    builder = main_window.report_request_builder
    assert isinstance(builder, QtReportRequestBuilder)

    main_window.ui.radioButton_wt_percent.setChecked(True)
    ui_conditions = dict(builder.get_current_conditions())

    assert ui_conditions["Basis"] == "Wt %"
    assert ui_conditions["Model / EOS"] == "PRMIX"

    main_window.last_result_data = _sample_result_mapping(
        basis="Mol %",
        model_display="IAPWS-95",
    )
    result_conditions = dict(builder.get_current_conditions())

    assert result_conditions["Basis"] == "Mol %"
    assert result_conditions["Model / EOS"] == "IAPWS-95"


def test_report_request_builder_collects_composition_rows(main_window):
    _select_component(main_window, "methane")
    table = main_window.ui.tableWidget
    table.item(0, 1).setText("12.5")
    table.item(0, 2).setText("not numeric")

    rows = main_window.report_request_builder.get_input_composition_rows()

    assert rows[0] == {
        "Component": "methane",
        "Mol %": 12.5,
        "Wt %": "not numeric",
    }
    assert rows[-1]["Component"] == "Total"
    assert rows[-1]["Mol %"] == 12.5


def test_report_request_builder_builds_export_request(main_window, tmp_path):
    _select_component(main_window, "methane")
    main_window.ui.tableWidget.item(0, 1).setText("100")
    main_window.last_result_data = _sample_result_mapping()

    exported_at = datetime(2026, 6, 17, 9, 15, 0)
    request = main_window.report_request_builder.build_export_request(
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


def test_report_request_builder_builds_projection_rows(main_window):
    result_data = _sample_result_mapping(
        warnings=[PRMIX_WATER_WARNING],
        missing_lhv=["unknown"],
    )

    result_rows = main_window.report_request_builder.build_results_rows(
        result_data
    )
    warning_rows = main_window.report_request_builder.build_report_warning_rows(
        result_data
    )

    assert result_rows[0]["Property"] == "Average molecular weight"
    assert warning_rows == [
        {
            "Warning Type": "Thermo",
            "Details": PRMIX_WATER_WARNING,
        },
        {
            "Warning Type": "LHV",
            "Details": "LHV Warning: No data for: unknown",
        },
    ]
