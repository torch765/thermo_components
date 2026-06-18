"""FastAPI routes for calculation workflows."""

from typing import cast

from fastapi import APIRouter, HTTPException, Request, status

from .calculation import (
    WebCalculationDependencies,
    execute_web_calculation,
)
from .schemas import CalculationRequestSchema, CalculationResponseSchema


router = APIRouter(prefix="/api", tags=["calculations"])


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
