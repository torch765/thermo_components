from datetime import datetime

from openpyxl import load_workbook

from thermo_components.adapters.reporting import OpenPyxlReportExporter
from thermo_components.application.dto import (
    ReportCompositionRow,
    ReportConditionRow,
    ReportExportRequest,
    ReportProjection,
    ReportResultRow,
    ReportWarningRow,
)
from thermo_components.application.ports import ReportExporter


def test_openpyxl_report_exporter_writes_prepared_report(tmp_path):
    report_path = tmp_path / "report.xlsx"
    exporter = OpenPyxlReportExporter()
    request = ReportExportRequest(
        report_path=report_path,
        exported_at=datetime(2026, 6, 16, 10, 30, 0),
        conditions=(ReportConditionRow("Basis", "Mol %"),),
        input_rows=(ReportCompositionRow("methane", 100.0, 100.0),),
        projection=ReportProjection(
            result_rows=(
                ReportResultRow(
                    "Density @ selected conditions",
                    0.6571965,
                    "kg/m³",
                    "Phase: Vapor",
                ),
            ),
            warning_rows=(
                ReportWarningRow("Thermo", "Example warning"),
            ),
        ),
    )

    exported_path = exporter.export(request)

    assert isinstance(exporter, ReportExporter)
    assert exported_path == report_path
    workbook = load_workbook(report_path)
    worksheet = workbook["Results"]

    assert worksheet["A1"].value == "THERMO CALCULATOR REPORT"
    assert worksheet["A2"].value == "Exported: 2026-06-16 10:30:00"
    assert worksheet["A4"].value == "Conditions"
    assert worksheet["A6"].value == "Basis"
    assert worksheet["B6"].value == "Mol %"
    assert worksheet["A8"].value == "Input Composition"
    assert worksheet["A10"].value == "methane"
    assert worksheet["B10"].value == 100.0
    assert worksheet["B10"].number_format == "0.0000"
    assert worksheet["A12"].value == "Warnings"
    assert worksheet["A14"].value == "Thermo"
    assert worksheet["B14"].value == "Example warning"
    assert worksheet["A16"].value == "Results"
    assert worksheet["A18"].value == "Density @ selected conditions"
    assert worksheet["B18"].value == 0.6571965
    assert worksheet["B18"].number_format == "0.000"
