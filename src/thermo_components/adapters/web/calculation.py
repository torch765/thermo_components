"""Shared calculation handling for JSON and HTML web adapters."""

from collections.abc import Callable, Mapping
from typing import Protocol

from thermo_components.application.dto import (
    DeriveCompositionRequest,
    DeriveCompositionResponse,
    PropertyCalculationRequest,
)
from thermo_components.application.services import (
    CalculationSessionRequest,
    CalculationSessionResponse,
)
from thermo_components.domain.composition import CompositionBasis
from thermo_components.domain.conditions import atm_to_pa, celsius_to_kelvin

from .presenters import present_calculation_session
from .schemas import CalculationRequestSchema, CalculationResponseSchema


class _DeriveCompositionWorkflow(Protocol):
    def execute(
        self,
        request: DeriveCompositionRequest,
    ) -> DeriveCompositionResponse: ...


class _CalculationSessionWorkflow(Protocol):
    def execute(
        self,
        request: CalculationSessionRequest,
    ) -> CalculationSessionResponse: ...


class WebCalculationDependencies(Protocol):
    lhv_database: Mapping[str, float]
    derive_composition_use_case: _DeriveCompositionWorkflow
    calculation_session_factory: Callable[
        [str],
        _CalculationSessionWorkflow,
    ]


def execute_web_calculation(
    payload: CalculationRequestSchema,
    dependencies: WebCalculationDependencies,
) -> CalculationResponseSchema:
    """Run the shared application workflow for a validated web request."""

    component_names = tuple(
        component.name for component in payload.components
    )
    active_percentages = tuple(
        component.percentage for component in payload.components
    )
    inactive_percentages = (
        dependencies.derive_composition_use_case.execute(
            DeriveCompositionRequest(
                component_names=component_names,
                active_percentages=active_percentages,
                basis=payload.basis,
            )
        ).percentages
    )
    if any(value is None for value in inactive_percentages):
        raise ValueError("The inactive composition basis could not be derived.")

    if payload.basis == CompositionBasis.MOLE_PERCENT.value:
        mole_percents = active_percentages
        weight_percents = tuple(
            float(value) for value in inactive_percentages
        )
    else:
        mole_percents = tuple(
            float(value) for value in inactive_percentages
        )
        weight_percents = active_percentages

    session = dependencies.calculation_session_factory(payload.model)
    response = session.execute(
        CalculationSessionRequest(
            calculation_request=PropertyCalculationRequest.from_sequences(
                component_names=component_names,
                mole_percents=mole_percents,
                weight_percents=weight_percents,
                basis=payload.basis,
                temperature_k=celsius_to_kelvin(payload.temperature_c),
                pressure_pa=atm_to_pa(payload.pressure_atm),
                pressure_atm=payload.pressure_atm,
            ),
            lhv_data_available=bool(dependencies.lhv_database),
            include_report_projection=payload.include_report_projection,
        )
    )
    return present_calculation_session(response)
