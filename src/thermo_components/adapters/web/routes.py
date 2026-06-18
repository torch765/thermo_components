"""FastAPI routes for calculation workflows."""

from typing import Protocol, cast

from fastapi import APIRouter, HTTPException, Request, status

from thermo_components.application.dto import (
    DeriveCompositionRequest,
    NormalizeCompositionRequest,
    NormalizeCompositionResponse,
)

from .calculation import (
    WebCalculationDependencies,
    execute_web_calculation,
)
from .schemas import (
    CalculationRequestSchema,
    CalculationResponseSchema,
    CompositionDeriveRequestSchema,
    CompositionDeriveResponseSchema,
    CompositionNormalizeRequestSchema,
    CompositionNormalizeResponseSchema,
)


router = APIRouter(prefix="/api", tags=["calculations"])


class _NormalizeCompositionWorkflow(Protocol):
    def execute(
        self,
        request: NormalizeCompositionRequest,
    ) -> NormalizeCompositionResponse: ...


class _CompositionDependencies(WebCalculationDependencies, Protocol):
    normalize_composition_use_case: _NormalizeCompositionWorkflow


@router.post(
    "/calculations",
    response_model=CalculationResponseSchema,
    status_code=status.HTTP_200_OK,
)
def calculate_properties(
    payload: CalculationRequestSchema,
    request: Request,
) -> CalculationResponseSchema:
    dependencies = cast(
        WebCalculationDependencies,
        request.app.state.dependencies,
    )

    try:
        return execute_web_calculation(payload, dependencies)
    except ValueError as exc:
        raise HTTPException(
            status_code=status.HTTP_422_UNPROCESSABLE_CONTENT,
            detail=str(exc),
        ) from exc


@router.post(
    "/compositions/derive",
    response_model=CompositionDeriveResponseSchema,
    tags=["composition"],
)
def derive_composition(
    payload: CompositionDeriveRequestSchema,
    request: Request,
) -> CompositionDeriveResponseSchema:
    dependencies = cast(
        WebCalculationDependencies,
        request.app.state.dependencies,
    )
    response = dependencies.derive_composition_use_case.execute(
        DeriveCompositionRequest(
            component_names=tuple(
                component.name for component in payload.components
            ),
            active_percentages=tuple(
                component.percentage for component in payload.components
            ),
            basis=payload.basis,
        )
    )
    return CompositionDeriveResponseSchema(
        percentages=list(response.percentages)
    )


@router.post(
    "/compositions/normalize",
    response_model=CompositionNormalizeResponseSchema,
    tags=["composition"],
)
def normalize_composition(
    payload: CompositionNormalizeRequestSchema,
    request: Request,
) -> CompositionNormalizeResponseSchema:
    dependencies = cast(
        _CompositionDependencies,
        request.app.state.dependencies,
    )
    try:
        response = dependencies.normalize_composition_use_case.execute(
            NormalizeCompositionRequest(
                percentages=tuple(payload.percentages),
                precision=payload.precision,
            )
        )
    except ValueError as exc:
        raise HTTPException(
            status_code=status.HTTP_422_UNPROCESSABLE_CONTENT,
            detail=str(exc),
        ) from exc
    return CompositionNormalizeResponseSchema(
        percentages=list(response.percentages)
    )
