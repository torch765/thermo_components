"""Report projection workflow independent of the output file format."""

from thermo_components.application.dto import (
    ReportPreparationRequest,
    ReportProjection,
    ReportResultRow,
    ReportWarningRow,
)
from thermo_components.domain.lhv import build_lhv_display_values
from thermo_components.domain.results import build_density_note
from thermo_components.domain.thermo_routes import PURE_WATER_ROUTE


class PrepareReportUseCase:
    def execute(self, request: ReportPreparationRequest) -> ReportProjection:
        calculation = request.calculation
        warning_rows = [
            ReportWarningRow("Thermo", warning)
            for warning in calculation.warnings
        ]
        if calculation.missing_lhv:
            warning_rows.append(
                ReportWarningRow(
                    "LHV",
                    "LHV Warning: No data for: "
                    + ", ".join(calculation.missing_lhv),
                )
            )

        result_rows: list[ReportResultRow] = [
            ReportResultRow(
                "Average molecular weight",
                calculation.average_molecular_weight,
                "g/mol",
            )
        ]
        self._add_selected_density_rows(result_rows, calculation)
        self._add_reference_density_row(
            result_rows,
            "Density @ normal conditions",
            calculation.normal_density,
        )
        self._add_reference_density_row(
            result_rows,
            "Density @ standard conditions",
            calculation.standard_density,
        )
        self._add_bubble_point_row(result_rows, calculation)
        self._add_lhv_rows(
            result_rows,
            calculation,
            request.lhv_data_available,
        )

        return ReportProjection(
            result_rows=tuple(result_rows),
            warning_rows=tuple(warning_rows),
        )

    @staticmethod
    def _add_selected_density_rows(result_rows, calculation) -> None:
        density = calculation.selected_density
        if density.phase:
            result_rows.append(
                ReportResultRow(
                    "Phase @ selected conditions",
                    density.phase,
                )
            )
        else:
            result_rows.append(
                ReportResultRow(
                    "Phase @ selected conditions",
                    "Could not determine",
                    notes=density.error or "",
                )
            )

        if density.error:
            result_rows.append(
                ReportResultRow(
                    "Density @ selected conditions",
                    "N.A.",
                    "kg/m\u00b3",
                    density.error,
                )
            )
        elif (
            density.phase == "Two-Phase"
            and isinstance(density.value, tuple)
            and len(density.value) == 2
        ):
            liquid_density, vapor_density = density.value
            result_rows.extend(
                [
                    ReportResultRow(
                        "Density @ selected conditions (Liq)",
                        liquid_density
                        if liquid_density is not None
                        else "N.A.",
                        "kg/m\u00b3",
                    ),
                    ReportResultRow(
                        "Density @ selected conditions (Vap)",
                        vapor_density
                        if vapor_density is not None
                        else "N.A.",
                        "kg/m\u00b3",
                    ),
                ]
            )
        elif density.value is not None:
            result_rows.append(
                ReportResultRow(
                    "Density @ selected conditions",
                    density.value,
                    "kg/m\u00b3",
                )
            )
        else:
            result_rows.append(
                ReportResultRow(
                    "Density @ selected conditions",
                    "N.A.",
                    "kg/m\u00b3",
                )
            )

    @staticmethod
    def _add_reference_density_row(result_rows, property_name, density) -> None:
        result_rows.append(
            ReportResultRow(
                property_name,
                density.scalar_kg_m3
                if density.scalar_kg_m3 is not None
                else "N.A.",
                "kg/m\u00b3",
                build_density_note(density.phase, density.error),
            )
        )

    @staticmethod
    def _add_bubble_point_row(result_rows, calculation) -> None:
        label = (
            f"Bubble point @ {calculation.pressure_atm:g} atm"
            if calculation.pressure_atm is not None
            else "Bubble point"
        )
        note = (
            "IAPWS-95 saturation temperature for pure water."
            if calculation.thermo_route == PURE_WATER_ROUTE
            else ""
        )
        if calculation.bubble_point_error:
            result_rows.append(
                ReportResultRow(
                    label,
                    "N/A",
                    "\u00b0C",
                    calculation.bubble_point_error,
                )
            )
        elif calculation.bubble_point_c is not None:
            result_rows.append(
                ReportResultRow(
                    label,
                    calculation.bubble_point_c,
                    "\u00b0C",
                    note,
                )
            )
        else:
            result_rows.append(
                ReportResultRow(
                    label,
                    "Calculation Failed",
                    "\u00b0C",
                    note,
                )
            )

    @staticmethod
    def _add_lhv_rows(
        result_rows,
        calculation,
        lhv_data_available: bool,
    ) -> None:
        if not lhv_data_available:
            result_rows.append(
                ReportResultRow(
                    "LHV data",
                    "N/A",
                    notes="LHV database not loaded.",
                )
            )
            return

        values = build_lhv_display_values(
            calculation.mixture_lhv_mj_nm3,
            calculation.average_molecular_weight,
        )
        result_rows.append(
            ReportResultRow(
                "Mixture LHV",
                values["volumetric"]["MJ/Nm\u00b3"],
                "MJ/Nm\u00b3",
                "Base value",
            )
        )
        for unit in (
            "kcal/Nm\u00b3",
            "MMkcal/Nm\u00b3",
            "GJ/Nm\u00b3",
            "MMBtu/Nm\u00b3",
        ):
            result_rows.append(
                ReportResultRow(
                    "Mixture LHV",
                    values["volumetric"][unit],
                    unit,
                )
            )

        mass_note = (
            "Derived from mixture MW and 22.414 Nm\u00b3/kmol "
            "at 0 \u00b0C and 1 atm."
        )
        for index, unit in enumerate(
            (
                "MJ/kg",
                "MJ/t",
                "GJ/kg",
                "GJ/t",
                "kcal/kg",
                "kcal/t",
                "MMkcal/kg",
                "MMkcal/t",
                "MMBtu/kg",
                "MMBtu/t",
            )
        ):
            result_rows.append(
                ReportResultRow(
                    "Mixture LHV (Mass basis)",
                    values["mass_basis"][unit],
                    unit,
                    mass_note if index == 0 else "",
                )
            )
