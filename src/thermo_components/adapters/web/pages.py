"""Server-rendered calculator pages."""

from itertools import zip_longest
from pathlib import Path
from typing import cast

from fastapi import APIRouter, Request, status
from fastapi.responses import HTMLResponse
from fastapi.templating import Jinja2Templates
from pydantic import ValidationError

from thermo_components.domain.composition import MOLECULAR_WEIGHTS
from thermo_components.domain.flow_units import FLOW_UNIT_ORDER

from .calculation import (
    WebCalculationDependencies,
    execute_web_calculation,
)
from .schemas import CalculationRequestSchema


WEB_ROOT = Path(__file__).resolve().parent
templates = Jinja2Templates(directory=WEB_ROOT / "templates")
page_router = APIRouter(include_in_schema=False)

COMPONENT_OPTIONS = tuple(
    {
        "value": component,
        "label": component.title(),
    }
    for component, molecular_weight in MOLECULAR_WEIGHTS.items()
    if component.strip() and molecular_weight > 0
)


@page_router.get("/", response_class=HTMLResponse, name="home")
@page_router.get(
    "/calculator",
    response_class=HTMLResponse,
    name="calculator",
)
def calculator_page(request: Request) -> HTMLResponse:
    return _render_calculator(
        request,
        form_state=_default_form_state(),
    )


@page_router.get(
    "/flow",
    response_class=HTMLResponse,
    name="flow",
)
def flow_page(request: Request) -> HTMLResponse:
    return templates.TemplateResponse(
        request=request,
        name="flow.html",
        context={
            "flow_units": FLOW_UNIT_ORDER,
        },
    )


@page_router.post(
    "/calculator",
    response_class=HTMLResponse,
    name="calculate_page",
)
async def submit_calculator(request: Request) -> HTMLResponse:
    form = await request.form()
    form_state = _form_state_from_submission(form)

    try:
        payload = _payload_from_form_state(form_state)
        dependencies = cast(
            WebCalculationDependencies,
            request.app.state.dependencies,
        )
        result = execute_web_calculation(payload, dependencies)
    except (ValidationError, ValueError) as exc:
        return _render_calculator(
            request,
            form_state=form_state,
            errors=_format_errors(exc),
            status_code=status.HTTP_422_UNPROCESSABLE_CONTENT,
        )

    return _render_calculator(
        request,
        form_state=_form_state_with_result(form_state, result),
        result=result,
    )


def _render_calculator(
    request: Request,
    *,
    form_state: dict,
    result=None,
    errors: tuple[str, ...] = (),
    status_code: int = status.HTTP_200_OK,
) -> HTMLResponse:
    return templates.TemplateResponse(
        request=request,
        name="calculator.html",
        context={
            "component_options": COMPONENT_OPTIONS,
            "errors": errors,
            "form_state": form_state,
            "result": result,
        },
        status_code=status_code,
    )


def _default_form_state() -> dict:
    return {
        "rows": [
            {
                "name": "methane",
                "mole_percentage": "100",
                "weight_percentage": "100",
            },
            {
                "name": "",
                "mole_percentage": "",
                "weight_percentage": "",
            },
        ],
        "basis": "Mol %",
        "temperature_c": "25",
        "pressure_atm": "1",
        "model": "PRMIX",
        "include_report_projection": False,
    }


def _form_state_from_submission(form) -> dict:
    component_names = form.getlist("component_name")
    mole_percentages = form.getlist("component_mole_percentage")
    weight_percentages = form.getlist("component_weight_percentage")
    rows = [
        {
            "name": str(name).strip(),
            "mole_percentage": str(mole_percentage).strip(),
            "weight_percentage": str(weight_percentage).strip(),
        }
        for name, mole_percentage, weight_percentage in zip_longest(
            component_names,
            mole_percentages,
            weight_percentages,
            fillvalue="",
        )
    ]
    if not rows:
        rows = [
            {
                "name": "",
                "mole_percentage": "",
                "weight_percentage": "",
            }
        ]

    return {
        "rows": rows,
        "basis": str(form.get("basis", "Mol %")),
        "temperature_c": str(form.get("temperature_c", "")),
        "pressure_atm": str(form.get("pressure_atm", "")),
        "model": "PRMIX",
        "include_report_projection": (
            form.get("include_report_projection") == "on"
        ),
    }


def _payload_from_form_state(form_state: dict) -> CalculationRequestSchema:
    components = []
    for row in form_state["rows"]:
        name = row["name"]
        percentage = (
            row["mole_percentage"]
            if form_state["basis"] == "Mol %"
            else row["weight_percentage"]
        )
        if not name and not percentage:
            continue
        if not name or not percentage:
            raise ValueError(
                "Each composition row needs both a component and percentage."
            )
        try:
            numeric_percentage = float(percentage)
        except ValueError as exc:
            raise ValueError(
                f"Invalid percentage for {name or 'component'}."
            ) from exc
        components.append(
            {
                "name": name,
                "percentage": numeric_percentage,
            }
        )

    return CalculationRequestSchema.model_validate(
        {
            "components": components,
            "basis": form_state["basis"],
            "temperature_c": form_state["temperature_c"],
            "pressure_atm": form_state["pressure_atm"],
            "model": form_state["model"],
            "include_report_projection": (
                form_state["include_report_projection"]
            ),
        }
    )


def _form_state_with_result(form_state: dict, result) -> dict:
    calculation = result.calculation
    rows = [
        {
            "name": component_name,
            "mole_percentage": _format_percentage(
                calculation.mole_percents[index]
            ),
            "weight_percentage": _format_percentage(
                calculation.weight_percents[index]
            ),
        }
        for index, component_name in enumerate(
            calculation.component_names
        )
    ]
    return {
        **form_state,
        "rows": rows,
    }


def _format_percentage(value: float) -> str:
    return f"{value:.4f}".rstrip("0").rstrip(".")


def _format_errors(exc: ValidationError | ValueError) -> tuple[str, ...]:
    if isinstance(exc, ValidationError):
        return tuple(error["msg"] for error in exc.errors())
    return (str(exc),)
