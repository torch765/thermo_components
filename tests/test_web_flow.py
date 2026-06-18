import re

import pytest
from fastapi.testclient import TestClient

from thermo_components.adapters.web.app import create_app
from thermo_components.application.dto import FlowConversionResponse
from thermo_components.bootstrap.web import (
    WebDependencies,
    build_web_dependencies,
)
from thermo_components.domain.flow_units import FLOW_UNIT_ORDER


class FakeFlowConversionUseCase:
    def __init__(self, response=None, error=None):
        self.response = response
        self.error = error
        self.requests = []

    def execute(self, request):
        self.requests.append(request)
        if self.error is not None:
            raise self.error
        return self.response


def build_dependencies(flow_use_case):
    base = build_web_dependencies(lhv_data={})
    return WebDependencies(
        lhv_database=base.lhv_database,
        derive_composition_use_case=base.derive_composition_use_case,
        calculation_session_factory=base.calculation_session_factory,
        normalize_composition_use_case=base.normalize_composition_use_case,
        convert_flow_use_case=flow_use_case,
    )


def test_flow_conversion_endpoint_maps_request_to_use_case():
    use_case = FakeFlowConversionUseCase(
        response=FlowConversionResponse(
            value=12.0,
            display_value="12",
        )
    )
    client = TestClient(create_app(build_dependencies(use_case)))

    response = client.post(
        "/api/flow-conversions",
        json={
            "input_text": "10",
            "from_unit": "Nm3/h",
            "to_unit": "Sm3/h",
            "normal_density_kg_m3": 1.2,
            "standard_density_kg_m3": 1.0,
        },
    )

    assert response.status_code == 200
    assert response.json() == {
        "value": 12.0,
        "display_value": "12",
    }
    request = use_case.requests[0]
    assert request.input_text == "10"
    assert request.from_unit == "Nm3/h"
    assert request.to_unit == "Sm3/h"
    assert request.normal_density_kg_m3 == pytest.approx(1.2)
    assert request.standard_density_kg_m3 == pytest.approx(1.0)


def test_flow_conversion_endpoint_runs_density_independent_conversion():
    client = TestClient(
        create_app(build_web_dependencies(lhv_data={}))
    )

    response = client.post(
        "/api/flow-conversions",
        json={
            "input_text": "1",
            "from_unit": "kg/h",
            "to_unit": "t/d",
        },
    )

    assert response.status_code == 200
    assert response.json()["value"] == pytest.approx(0.024)
    assert response.json()["display_value"] == "0.024"


def test_flow_conversion_endpoint_uses_both_reference_densities():
    client = TestClient(
        create_app(build_web_dependencies(lhv_data={}))
    )

    response = client.post(
        "/api/flow-conversions",
        json={
            "input_text": "10",
            "from_unit": "Nm3/h",
            "to_unit": "Sm3/h",
            "normal_density_kg_m3": 1.2,
            "standard_density_kg_m3": 1.0,
        },
    )

    assert response.status_code == 200
    assert response.json()["value"] == pytest.approx(12.0)
    assert response.json()["display_value"] == "12"


def test_flow_conversion_endpoint_returns_density_requirement():
    client = TestClient(
        create_app(build_web_dependencies(lhv_data={}))
    )

    response = client.post(
        "/api/flow-conversions",
        json={
            "input_text": "1",
            "from_unit": "kg/h",
            "to_unit": "Nm3/h",
        },
    )

    assert response.status_code == 422
    assert response.json() == {"detail": "Normal density required"}


def test_flow_conversion_endpoint_rejects_unsupported_unit():
    client = TestClient(
        create_app(build_web_dependencies(lhv_data={}))
    )

    response = client.post(
        "/api/flow-conversions",
        json={
            "input_text": "1",
            "from_unit": "kg/min",
            "to_unit": "t/d",
        },
    )

    assert response.status_code == 422
    assert "Unsupported flow unit" in response.text


def test_flow_page_renders_all_units_and_default_conversion():
    client = TestClient(
        create_app(build_web_dependencies(lhv_data={}))
    )

    response = client.get("/flow")

    assert response.status_code == 200
    assert "Flow conversion workspace" in response.text
    assert 'data-flow-input' in response.text
    assert 'data-normal-density' in response.text
    assert 'data-standard-density' in response.text
    assert 'src="http://testserver/static/flow.js"' in response.text
    assert ">0.024</output>" in response.text
    for unit in FLOW_UNIT_ORDER:
        assert f'value="{unit}"' in response.text
    assert re.search(
        r'<option value="kg/h" selected>',
        response.text,
    )
    assert re.search(
        r'<option value="t/d" selected>',
        response.text,
    )


def test_calculator_and_flow_pages_link_to_each_other():
    client = TestClient(
        create_app(build_web_dependencies(lhv_data={}))
    )

    calculator = client.get("/calculator")
    flow = client.get("/flow")

    assert 'href="http://testserver/flow"' in calculator.text
    assert 'href="http://testserver/calculator"' in flow.text
    assert "thermo-components.flow-densities" in client.get(
        "/static/calculator.js"
    ).text
    assert "thermo-components.flow-densities" in client.get(
        "/static/flow.js"
    ).text


def test_calculator_result_exposes_flow_densities_for_session_handoff():
    client = TestClient(
        create_app(
            build_web_dependencies(
                lhv_data={"methane": 35.8},
            )
        )
    )

    response = client.post(
        "/calculator",
        data={
            "component_name": "methane",
            "component_mole_percentage": "100",
            "component_weight_percentage": "100",
            "basis": "Mol %",
            "temperature_c": "25",
            "pressure_atm": "1",
            "model": "PRMIX",
        },
    )

    assert response.status_code == 200
    assert "data-flow-density-state" in response.text
    assert re.search(
        r'data-normal-density="[0-9.]+(?:e[+-]?[0-9]+)?"',
        response.text,
        re.IGNORECASE,
    )
    assert re.search(
        r'data-standard-density="[0-9.]+(?:e[+-]?[0-9]+)?"',
        response.text,
        re.IGNORECASE,
    )
