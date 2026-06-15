"""Thermodynamic property calculation workflow."""

from collections.abc import Mapping

from thermo_components.application.dto import (
    DensityCalculationResult,
    PropertyCalculationRequest,
    PropertyCalculationResponse,
)
from thermo_components.application.ports import ThermoPropertyGateway
from thermo_components.domain.composition import (
    CompositionBasis,
    calculate_mixture_molecular_weight,
    percentages_to_mole_fractions,
)
from thermo_components.domain.conditions import (
    NORMAL_CONDITION,
    STANDARD_CONDITION,
)
from thermo_components.domain.lhv import calculate_mixture_lhv
from thermo_components.domain.thermo_routes import select_thermo_route
from thermo_components.domain.warnings import build_thermo_warning_messages


class CalculatePropertiesUseCase:
    def __init__(
        self,
        calculator: ThermoPropertyGateway,
        lhv_data: Mapping[str, float],
    ):
        self._calculator = calculator
        self._lhv_data = lhv_data

    def execute(
        self,
        request: PropertyCalculationRequest,
    ) -> PropertyCalculationResponse:
        self._validate_request(request)
        active_percentages = (
            request.mole_percents
            if request.basis == CompositionBasis.MOLE_PERCENT.value
            else request.weight_percents
        )
        mole_fractions = percentages_to_mole_fractions(
            request.component_names,
            active_percentages,
            request.basis,
        )
        self._calculator.set_components(mole_fractions)

        route = select_thermo_route(
            request.component_names,
            request.mole_percents,
            request.weight_percents,
            request.basis,
            self._calculator.eos,
        )
        route_id = route["route_id"]

        selected_density = DensityCalculationResult.from_raw(
            *self._calculator.calculate_density_for_route(
                request.temperature_k,
                request.pressure_pa,
                route_id,
            )
        )
        normal_density = DensityCalculationResult.from_raw(
            *self._calculator.calculate_density_for_route(
                NORMAL_CONDITION.temperature_k,
                NORMAL_CONDITION.pressure_pa,
                route_id,
            )
        )
        standard_density = DensityCalculationResult.from_raw(
            *self._calculator.calculate_density_for_route(
                STANDARD_CONDITION.temperature_k,
                STANDARD_CONDITION.pressure_pa,
                route_id,
            )
        )
        bubble_point_c, bubble_point_error = (
            self._calculator.calculate_bubble_point_for_route(
                request.pressure_pa,
                route_id,
            )
        )
        average_molecular_weight = calculate_mixture_molecular_weight(
            mole_fractions
        )
        mixture_lhv, missing_lhv = calculate_mixture_lhv(
            mole_fractions,
            self._lhv_data,
        )
        warning_data = {
            "phase": selected_density.phase,
            "density_normal_phase": normal_density.phase,
            "density_standard_phase": standard_density.phase,
        }
        warnings = build_thermo_warning_messages(
            request.component_names,
            request.mole_percents,
            request.weight_percents,
            request.basis,
            route_id,
            warning_data,
        )

        return PropertyCalculationResponse(
            average_molecular_weight=average_molecular_weight,
            component_names=request.component_names,
            mole_percents=request.mole_percents,
            weight_percents=request.weight_percents,
            thermo_route=route_id,
            model_display=route["model_display"],
            selected_density=selected_density,
            normal_density=normal_density,
            standard_density=standard_density,
            bubble_point_c=bubble_point_c,
            bubble_point_error=bubble_point_error,
            mixture_lhv_mj_nm3=mixture_lhv,
            missing_lhv=tuple(missing_lhv),
            basis=request.basis,
            eos=self._calculator.eos,
            pressure_atm=request.pressure_atm,
            warnings=tuple(warnings),
        )

    @staticmethod
    def _validate_request(request: PropertyCalculationRequest) -> None:
        if not request.component_names:
            raise ValueError("No components selected.")

        expected_length = len(request.component_names)
        if (
            len(request.mole_percents) != expected_length
            or len(request.weight_percents) != expected_length
        ):
            raise ValueError("Composition rows are incomplete.")

        if request.basis == CompositionBasis.MOLE_PERCENT.value:
            active_percentages = request.mole_percents
        elif request.basis == CompositionBasis.WEIGHT_PERCENT.value:
            active_percentages = request.weight_percents
        else:
            raise ValueError(f"Unsupported composition basis: {request.basis}")

        total_percent = sum(active_percentages)
        if abs(total_percent - 100.0) > 1e-4 or total_percent == 0:
            raise ValueError(f"Invalid {request.basis} input.")
        if request.temperature_k <= 0:
            raise ValueError("Temperature must be positive.")
        if request.pressure_pa <= 0 or request.pressure_atm <= 0:
            raise ValueError("Pressure must be positive.")
