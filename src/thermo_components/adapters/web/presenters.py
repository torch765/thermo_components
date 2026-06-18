"""Translate application responses into web response schemas."""

from thermo_components.application.dto import (
    DensityCalculationResult,
    ReportProjection,
)
from thermo_components.application.services import CalculationSessionResponse

from .schemas import (
    CalculationResponseSchema,
    DensityResultSchema,
    FlowDensityStateSchema,
    PropertyCalculationSchema,
    ReportProjectionSchema,
    ReportResultRowSchema,
    ReportWarningRowSchema,
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
        warnings=list(response.warning_messages),
        report_projection=_present_report_projection(
            response.report_projection
        ),
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
