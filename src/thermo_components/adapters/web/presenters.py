"""Translate application responses into web response schemas."""

from thermo_components.application.dto import (
    DensityCalculationResult,
    ReportProjection,
)
from thermo_components.application.services import CalculationSessionResponse
from thermo_components.domain.lhv import (
    build_lhv_display_values,
    format_lhv_display_value,
)

from .schemas import (
    CalculationResponseSchema,
    DensityResultSchema,
    FlowDensityStateSchema,
    LhvDisplaySchema,
    LhvValueSchema,
    PropertyCalculationSchema,
    ReportProjectionSchema,
    ReportResultRowSchema,
    ReportWarningRowSchema,
)

VOLUMETRIC_LHV_UNITS = (
    "MJ/Nm\u00b3",
    "kcal/Nm\u00b3",
    "MMkcal/Nm\u00b3",
    "GJ/Nm\u00b3",
    "MMBtu/Nm\u00b3",
)
MASS_BASIS_LHV_UNITS = (
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


def present_calculation_session(
    response: CalculationSessionResponse,
) -> CalculationResponseSchema:
    """Build the HTTP response without leaking application dataclasses."""

    calculation = response.calculation
    return CalculationResponseSchema(
        calculation=PropertyCalculationSchema(
            average_molecular_weight=calculation.average_molecular_weight,
            component_names=list(calculation.component_names),
            mole_percents=list(calculation.mole_percents),
            weight_percents=list(calculation.weight_percents),
            thermo_route=calculation.thermo_route,
            model_display=calculation.model_display,
            selected_density=_present_density(calculation.selected_density),
            normal_density=_present_density(calculation.normal_density),
            standard_density=_present_density(calculation.standard_density),
            bubble_point_c=calculation.bubble_point_c,
            bubble_point_error=calculation.bubble_point_error,
            mixture_lhv_mj_nm3=calculation.mixture_lhv_mj_nm3,
            missing_lhv=list(calculation.missing_lhv),
            basis=calculation.basis,
            eos=calculation.eos,
            pressure_atm=calculation.pressure_atm,
            warnings=list(calculation.warnings),
        ),
        flow_densities=FlowDensityStateSchema(
            selected_density_kg_m3=(
                response.flow_densities.selected_density_kg_m3
            ),
            normal_density_kg_m3=(
                response.flow_densities.normal_density_kg_m3
            ),
            standard_density_kg_m3=(
                response.flow_densities.standard_density_kg_m3
            ),
        ),
        lhv=_present_lhv(response),
        warnings=list(response.warning_messages),
        report_projection=_present_report_projection(
            response.report_projection
        ),
    )


def _present_lhv(
    response: CalculationSessionResponse,
) -> LhvDisplaySchema:
    calculation = response.calculation
    values = build_lhv_display_values(
        calculation.mixture_lhv_mj_nm3,
        calculation.average_molecular_weight,
    )
    available = response.lhv_data_available
    return LhvDisplaySchema(
        available=available,
        volumetric=[
            _present_lhv_value(
                unit,
                values["volumetric"][unit] if available else None,
            )
            for unit in VOLUMETRIC_LHV_UNITS
        ],
        mass_basis=[
            _present_lhv_value(
                unit,
                values["mass_basis"][unit] if available else None,
            )
            for unit in MASS_BASIS_LHV_UNITS
        ],
    )


def _present_lhv_value(
    unit: str,
    value: float | None,
) -> LhvValueSchema:
    return LhvValueSchema(
        unit=unit,
        value=value,
        display_value=format_lhv_display_value(value),
    )


def _present_density(
    density: DensityCalculationResult,
) -> DensityResultSchema:
    return DensityResultSchema(
        value=density.value,
        phase=density.phase,
        error=density.error,
        scalar_kg_m3=density.scalar_kg_m3,
    )


def _present_report_projection(
    projection: ReportProjection | None,
) -> ReportProjectionSchema | None:
    if projection is None:
        return None
    return ReportProjectionSchema(
        result_rows=[
            ReportResultRowSchema(
                property_name=row.property_name,
                value=row.value,
                unit=row.unit,
                notes=row.notes,
            )
            for row in projection.result_rows
        ],
        warning_rows=[
            ReportWarningRowSchema(
                warning_type=row.warning_type,
                details=row.details,
            )
            for row in projection.warning_rows
        ],
    )
