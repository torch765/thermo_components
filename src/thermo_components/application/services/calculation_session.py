"""UI-independent calculation workflow facade."""

from dataclasses import dataclass
from typing import Protocol

from thermo_components.application.dto import (
    PropertyCalculationRequest,
    PropertyCalculationResponse,
    ReportPreparationRequest,
    ReportProjection,
)


class _CalculatePropertiesWorkflow(Protocol):
    def execute(
        self,
        request: PropertyCalculationRequest,
    ) -> PropertyCalculationResponse: ...


class _PrepareReportWorkflow(Protocol):
    def execute(
        self,
        request: ReportPreparationRequest,
    ) -> ReportProjection: ...


@dataclass(frozen=True)
class CalculationSessionRequest:
    """Request for a full calculation session independent of any UI layer."""

    calculation_request: PropertyCalculationRequest
    lhv_data_available: bool
    include_report_projection: bool = True


@dataclass(frozen=True)
class FlowDensityState:
    """Density values needed by follow-on flow conversion workflows."""

    selected_density_kg_m3: float | None
    normal_density_kg_m3: float | None
    standard_density_kg_m3: float | None


@dataclass(frozen=True)
class CalculationSessionResponse:
    """Result of a UI-independent calculation session."""

    calculation: PropertyCalculationResponse
    flow_densities: FlowDensityState
    warning_messages: tuple[str, ...]
    report_projection: ReportProjection | None = None


class CalculationSessionService:
    """Compose calculation and report projection without UI dependencies."""

    def __init__(
        self,
        calculate_properties_use_case: _CalculatePropertiesWorkflow,
        prepare_report_use_case: _PrepareReportWorkflow,
    ):
        self.calculate_properties_use_case = calculate_properties_use_case
        self.prepare_report_use_case = prepare_report_use_case

    def execute(
        self,
        request: CalculationSessionRequest,
    ) -> CalculationSessionResponse:
        calculation = self.calculate_properties_use_case.execute(
            request.calculation_request
        )
        report_projection = None
        if request.include_report_projection:
            report_projection = self.prepare_report_use_case.execute(
                ReportPreparationRequest(
                    calculation=calculation,
                    lhv_data_available=request.lhv_data_available,
                )
            )

        return CalculationSessionResponse(
            calculation=calculation,
            flow_densities=FlowDensityState(
                selected_density_kg_m3=(
                    calculation.selected_density.scalar_kg_m3
                ),
                normal_density_kg_m3=(
                    calculation.normal_density.scalar_kg_m3
                ),
                standard_density_kg_m3=(
                    calculation.standard_density.scalar_kg_m3
                ),
            ),
            warning_messages=calculation.warnings,
            report_projection=report_projection,
        )
