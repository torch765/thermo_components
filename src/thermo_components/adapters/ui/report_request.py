"""Qt report request builder for current UI state and latest calculation."""

from collections.abc import Callable

from thermo_components.application.dto import (
    ReportCompositionRow,
    ReportConditionRow,
    ReportExportRequest,
    ReportPreparationRequest,
    coerce_property_response,
)


class QtReportRequestBuilder:
    """Build report DTOs from Qt widgets and the latest calculation result."""

    def __init__(
        self,
        ui,
        prepare_report_use_case,
        *,
        lhv_data_available_provider: Callable[[], bool],
        latest_result_provider: Callable[[], object | None],
    ):
        self.ui = ui
        self.prepare_report_use_case = prepare_report_use_case
        self.lhv_data_available_provider = lhv_data_available_provider
        self.latest_result_provider = latest_result_provider

    def get_current_conditions(self) -> list[tuple[str, object]]:
        """Return the active calculation conditions for reporting."""
        latest_result = self.latest_result_provider()
        if latest_result:
            response = coerce_property_response(latest_result)
            basis = response.basis
            model_display = response.model_display
        else:
            basis = (
                "Mol %"
                if self.ui.radioButton_mol_percent.isChecked()
                else "Wt %"
            )
            model_display = self.ui.comboBox_select_EOS.currentText()

        return [
            ("Basis", basis),
            ("Temperature", self.ui.comboBox_select_temperature.currentText()),
            ("Pressure", self.ui.comboBox_select_pressure.currentText()),
            ("Model / EOS", model_display),
        ]

    def get_input_composition_rows(self) -> list[dict[str, object]]:
        """Return report-ready composition rows, including the Total row."""
        rows = []
        for row_index in range(self.ui.tableWidget.rowCount()):
            component_item = self.ui.tableWidget.item(row_index, 0)
            mol_item = self.ui.tableWidget.item(row_index, 1)
            wt_item = self.ui.tableWidget.item(row_index, 2)
            rows.append(
                {
                    "Component": (
                        component_item.text() if component_item else ""
                    ),
                    "Mol %": self.coerce_report_value(
                        mol_item.text() if mol_item else ""
                    ),
                    "Wt %": self.coerce_report_value(
                        wt_item.text() if wt_item else ""
                    ),
                }
            )
        return rows

    def build_report_projection(self, result_data):
        """Build the typed report projection used by display and export paths."""
        response = coerce_property_response(result_data)
        return self.prepare_report_use_case.execute(
            ReportPreparationRequest(
                calculation=response,
                lhv_data_available=self.lhv_data_available_provider(),
            )
        )

    def build_report_warning_rows(self, result_data) -> list[dict[str, object]]:
        """Build structured warning rows for the export report."""
        projection = self.build_report_projection(result_data)
        return [row.to_dict() for row in projection.warning_rows]

    def build_results_rows(self, result_data) -> list[dict[str, object]]:
        """Build structured result rows for report export without UI scraping."""
        projection = self.build_report_projection(result_data)
        return [row.to_dict() for row in projection.result_rows]

    def build_export_request(self, report_path, timestamp) -> ReportExportRequest:
        """Build the typed report export request for the latest calculation."""
        return ReportExportRequest(
            report_path=report_path,
            exported_at=timestamp,
            conditions=tuple(
                ReportConditionRow(setting, value)
                for setting, value in self.get_current_conditions()
            ),
            input_rows=tuple(
                ReportCompositionRow(
                    row["Component"],
                    row["Mol %"],
                    row["Wt %"],
                )
                for row in self.get_input_composition_rows()
            ),
            projection=self.build_report_projection(
                self.latest_result_provider()
            ),
        )

    @staticmethod
    def coerce_report_value(text: str):
        """Convert report cell text to float when possible."""
        cleaned = str(text).strip() if text is not None else ""
        if not cleaned:
            return ""
        try:
            return float(cleaned)
        except ValueError:
            return cleaned
