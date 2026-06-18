"""Build and stream temporary Excel reports for the web adapter."""

import os
import tempfile
from collections.abc import Callable
from datetime import datetime
from pathlib import Path
from typing import Protocol

from fastapi.responses import FileResponse
from starlette.background import BackgroundTask

from thermo_components.application.dto import (
    ReportCompositionRow,
    ReportConditionRow,
    ReportExportRequest,
    ReportProjection,
    ReportResultRow,
    ReportWarningRow,
)

from .calculation import WebCalculationDependencies
from .schemas import CalculationRequestSchema, CalculationResponseSchema


EXCEL_MEDIA_TYPE = (
    "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
)


class _ReportExporter(Protocol):
    def export(self, request: ReportExportRequest) -> Path: ...


class WebReportDependencies(WebCalculationDependencies, Protocol):
    report_exporter: _ReportExporter
    report_clock: Callable[[], datetime]


def build_report_export_request(
    payload: CalculationRequestSchema,
    result: CalculationResponseSchema,
    *,
    report_path: Path,
    exported_at: datetime,
) -> ReportExportRequest:
    """Translate a web calculation into the shared report export DTO."""

    projection = result.report_projection
    if projection is None:
        raise ValueError("Report projection was not prepared.")

    calculation = result.calculation
    composition_rows = [
        ReportCompositionRow(
            component=component_name,
            mole_percent=calculation.mole_percents[index],
            weight_percent=calculation.weight_percents[index],
        )
        for index, component_name in enumerate(calculation.component_names)
    ]
    composition_rows.append(
        ReportCompositionRow(
            component="Total",
            mole_percent=sum(calculation.mole_percents),
            weight_percent=sum(calculation.weight_percents),
        )
    )

    return ReportExportRequest(
        report_path=report_path,
        exported_at=exported_at,
        conditions=(
            ReportConditionRow("Basis", calculation.basis),
            ReportConditionRow(
                "Temperature",
                f"{payload.temperature_c:g} °C",
            ),
            ReportConditionRow(
                "Pressure",
                f"{payload.pressure_atm:g} atm",
            ),
            ReportConditionRow("Model / EOS", calculation.model_display),
        ),
        input_rows=tuple(composition_rows),
        projection=ReportProjection(
            result_rows=tuple(
                ReportResultRow(
                    property_name=row.property_name,
                    value=row.value,
                    unit=row.unit,
                    notes=row.notes,
                )
                for row in projection.result_rows
            ),
            warning_rows=tuple(
                ReportWarningRow(
                    warning_type=row.warning_type,
                    details=row.details,
                )
                for row in projection.warning_rows
            ),
        ),
    )


def create_report_download(
    payload: CalculationRequestSchema,
    result: CalculationResponseSchema,
    dependencies: WebReportDependencies,
) -> FileResponse:
    """Export to a temporary file and delete it after the response."""

    exported_at = dependencies.report_clock()
    descriptor, temp_name = tempfile.mkstemp(
        prefix="thermo-components-",
        suffix=".xlsx",
    )
    os.close(descriptor)
    temp_path = Path(temp_name)

    try:
        request = build_report_export_request(
            payload,
            result,
            report_path=temp_path,
            exported_at=exported_at,
        )
        exported_path = Path(dependencies.report_exporter.export(request))
    except Exception:
        temp_path.unlink(missing_ok=True)
        raise

    if exported_path != temp_path:
        temp_path.unlink(missing_ok=True)

    filename = (
        f"Thermo_Report_{exported_at.strftime('%Y%m%d_%H%M%S')}.xlsx"
    )
    return FileResponse(
        path=exported_path,
        filename=filename,
        media_type=EXCEL_MEDIA_TYPE,
        background=BackgroundTask(_delete_report, exported_path),
    )


def _delete_report(path: Path) -> None:
    path.unlink(missing_ok=True)
