from datetime import datetime

from density import MainWindow
from thermo_components.adapters.ui import QtReportExportController
from thermo_components.application.dto import (
    ReportExportRequest,
    ReportProjection,
)
from thermo_components.domain.thermo_routes import PRMIX_DEFAULT_ROUTE


class RecordingReportExporter:
    def __init__(self, exception=None):
        self.exception = exception
        self.requests = []

    def export(self, request):
        self.requests.append(request)
        if self.exception is not None:
            raise self.exception
        return request.report_path


def _build_request(report_path, exported_at):
    return ReportExportRequest(
        report_path=report_path,
        exported_at=exported_at,
        conditions=(),
        input_rows=(),
        projection=ReportProjection(result_rows=(), warning_rows=()),
    )


def _select_component(window, component_name):
    index = window.ui.comboBox_select_components.findText(component_name)
    assert index >= 0
    window.ui.comboBox_select_components.setCurrentIndex(index)


def _sample_result_mapping():
    return {
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


def test_report_controller_reports_no_results_without_exporting(
    tmp_path,
    monkeypatch,
):
    messages = []
    monkeypatch.setattr(
        "thermo_components.adapters.ui.report_controller.QMessageBox.information",
        lambda parent, title, message: messages.append((title, message)),
    )
    exporter = RecordingReportExporter()
    controller = QtReportExportController(
        object(),
        exporter,
        export_base_dir_provider=lambda: tmp_path,
    )

    result = controller.export_latest_report(
        has_results=False,
        build_request=_build_request,
    )

    assert result is None
    assert exporter.requests == []
    assert messages == [
        ("Print Results", "No calculation results available. Please click Go first.")
    ]


def test_report_controller_exports_and_opens_report(
    tmp_path,
    monkeypatch,
):
    messages = []
    opened_paths = []
    monkeypatch.setattr(
        "thermo_components.adapters.ui.report_controller.QMessageBox.information",
        lambda parent, title, message: messages.append((title, message)),
    )
    exporter = RecordingReportExporter()
    exported_at = datetime(2026, 6, 16, 10, 30, 0)
    controller = QtReportExportController(
        object(),
        exporter,
        export_base_dir_provider=lambda: tmp_path,
        clock=lambda: exported_at,
        file_opener=opened_paths.append,
    )

    result = controller.export_latest_report(
        has_results=True,
        build_request=_build_request,
    )

    expected_path = tmp_path / "Thermo_Report_20260616_103000.xlsx"
    assert result == str(expected_path)
    assert exporter.requests[0].report_path == expected_path
    assert exporter.requests[0].exported_at == exported_at
    assert messages == [
        ("Print Results", f"Excel report saved to:\n{expected_path}")
    ]
    assert opened_paths == [expected_path]


def test_report_controller_shows_openpyxl_install_hint(
    tmp_path,
    monkeypatch,
):
    warnings = []
    monkeypatch.setattr(
        "thermo_components.adapters.ui.report_controller.QMessageBox.warning",
        lambda parent, title, message: warnings.append((title, message)),
    )
    exporter = RecordingReportExporter(exception=ImportError("missing"))
    controller = QtReportExportController(
        object(),
        exporter,
        export_base_dir_provider=lambda: tmp_path,
        file_opener=lambda path: None,
    )

    result = controller.export_latest_report(
        has_results=True,
        build_request=_build_request,
    )

    assert result is None
    assert warnings == [
        (
            "Print Results",
            "Excel export requires the 'openpyxl' package.\n"
            "Install it in this environment and try again.\n\n"
            "Command:\npython -m pip install openpyxl",
        )
    ]


def test_report_controller_shows_export_failure(
    tmp_path,
    monkeypatch,
):
    critical_messages = []
    monkeypatch.setattr(
        "thermo_components.adapters.ui.report_controller.QMessageBox.critical",
        lambda parent, title, message: critical_messages.append(
            (title, message)
        ),
    )
    exporter = RecordingReportExporter(exception=RuntimeError("disk full"))
    controller = QtReportExportController(
        object(),
        exporter,
        export_base_dir_provider=lambda: tmp_path,
        file_opener=lambda path: None,
    )

    result = controller.export_latest_report(
        has_results=True,
        build_request=_build_request,
    )

    assert result is None
    assert critical_messages == [
        ("Print Results", "Failed to save Excel report:\ndisk full")
    ]


def test_main_window_report_export_uses_controller_and_request_builder(
    qt_app,
    tmp_path,
    monkeypatch,
):
    messages = []
    opened_paths = []
    monkeypatch.setattr(
        "thermo_components.adapters.ui.report_controller.QMessageBox.information",
        lambda parent, title, message: messages.append((title, message)),
    )
    exporter = RecordingReportExporter()
    exported_at = datetime(2026, 6, 16, 10, 30, 0)
    window = MainWindow(lhv_data={"methane": 35.8})
    try:
        _select_component(window, "methane")
        window.ui.tableWidget.item(0, 1).setText("100")
        window.last_result_data = _sample_result_mapping()
        window.report_export_controller.report_exporter = exporter
        window.report_export_controller.export_base_dir_provider = (
            lambda: tmp_path
        )
        window.report_export_controller.clock = lambda: exported_at
        window.report_export_controller.file_opener = opened_paths.append

        result = window.export_results_to_excel()

        expected_path = tmp_path / "Thermo_Report_20260616_103000.xlsx"
        request = exporter.requests[0]
        assert result == str(expected_path)
        assert request.report_path == expected_path
        assert request.exported_at == exported_at
        assert request.conditions[0].setting == "Basis"
        assert request.conditions[0].value == "Mol %"
        assert request.input_rows[0].component == "methane"
        assert request.projection.result_rows
        assert messages == [
            ("Print Results", f"Excel report saved to:\n{expected_path}")
        ]
        assert opened_paths == [expected_path]
    finally:
        window.close()
