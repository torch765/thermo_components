from datetime import datetime
from io import BytesIO

from fastapi.testclient import TestClient
from openpyxl import load_workbook

from thermo_components.adapters.reporting import OpenPyxlReportExporter
from thermo_components.adapters.web.app import create_app
from thermo_components.bootstrap.web import (
    WebDependencies,
    build_web_dependencies,
)


class RecordingReportExporter:
    def __init__(self):
        self.delegate = OpenPyxlReportExporter()
        self.request = None

    def export(self, request):
        self.request = request
        return self.delegate.export(request)


def methane_form(**overrides):
    form = {
        "component_name": "methane",
        "component_mole_percentage": "100",
        "component_weight_percentage": "100",
        "basis": "Mol %",
        "temperature_c": "25",
        "pressure_atm": "1",
        "model": "PRMIX",
    }
    form.update(overrides)
    return form


def build_report_client(exporter):
    base = build_web_dependencies(lhv_data={"methane": 35.8})
    dependencies = WebDependencies(
        lhv_database=base.lhv_database,
        derive_composition_use_case=base.derive_composition_use_case,
        calculation_session_factory=base.calculation_session_factory,
        normalize_composition_use_case=base.normalize_composition_use_case,
        convert_flow_use_case=base.convert_flow_use_case,
        report_exporter=exporter,
        report_clock=lambda: datetime(2026, 6, 18, 16, 45, 30),
    )
    return TestClient(create_app(dependencies))


def test_report_download_returns_temporary_excel_workbook():
    exporter = RecordingReportExporter()
    client = build_report_client(exporter)

    response = client.post(
        "/calculator/report",
        data=methane_form(),
    )

    assert response.status_code == 200
    assert response.headers["content-type"].startswith(
        "application/vnd.openxmlformats-officedocument"
    )
    assert (
        'filename="Thermo_Report_20260618_164530.xlsx"'
        in response.headers["content-disposition"]
    )
    workbook = load_workbook(BytesIO(response.content))
    worksheet = workbook["Results"]
    assert worksheet["A1"].value == "THERMO CALCULATOR REPORT"
    assert worksheet["A6"].value == "Basis"
    assert worksheet["B6"].value == "Mol %"
    assert worksheet["A7"].value == "Temperature"
    assert worksheet["B7"].value == "25 °C"
    assert worksheet["A13"].value == "methane"
    assert worksheet["B13"].value == 100.0
    assert worksheet["A14"].value == "Total"
    assert exporter.request.projection.result_rows
    assert not exporter.request.report_path.exists()


def test_report_download_rejects_invalid_calculation_input():
    exporter = RecordingReportExporter()
    client = build_report_client(exporter)

    response = client.post(
        "/calculator/report",
        data=methane_form(component_mole_percentage="90"),
    )

    assert response.status_code == 422
    assert response.headers["content-type"].startswith("text/html")
    assert "total must be 100" in response.text
    assert exporter.request is None
