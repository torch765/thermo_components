import re
from urllib.parse import urlencode

from fastapi.testclient import TestClient

from thermo_components.adapters.web.app import create_app
from thermo_components.bootstrap.web import build_web_dependencies


def build_client() -> TestClient:
    return TestClient(
        create_app(
            build_web_dependencies(
                lhv_data={"methane": 35.8},
            )
        )
    )


def methane_form(**overrides):
    form = {
        "component_name": "methane",
        "component_mole_percentage": "100",
        "component_weight_percentage": "100",
        "basis": "Mol %",
        "temperature_c": "25",
        "pressure_atm": "1",
        "model": "PRMIX",
    }
    form.update(overrides)
    return form


def test_calculator_page_renders_default_form():
    response = build_client().get("/")

    assert response.status_code == 200
    assert response.headers["content-type"].startswith("text/html")
    assert "Mixture properties" in response.text
    assert 'action="http://testserver/calculator"' in response.text
    assert re.search(
        r'<option\s+value="methane"\s+selected',
        response.text,
    )
    assert 'name="temperature_c"' in response.text
    assert 'value="25"' in response.text
    assert 'name="component_mole_percentage"' in response.text
    assert 'name="component_weight_percentage"' in response.text
    assert "Normalize active basis" in response.text


def test_calculator_alias_and_static_assets_are_available():
    client = build_client()

    page_response = client.get("/calculator")
    css_response = client.get("/static/calculator.css")
    script_response = client.get("/static/calculator.js")

    assert page_response.status_code == 200
    assert css_response.status_code == 200
    assert "--teal:" in css_response.text
    assert script_response.status_code == 200
    assert "data-composition-list" in script_response.text


def test_calculator_form_renders_real_methane_results():
    response = build_client().post(
        "/calculator",
        data=methane_form(),
    )

    assert response.status_code == 200
    assert "Calculated properties" in response.text
    assert "0.6572" in response.text
    assert "Vapor" in response.text
    assert "16.040" in response.text
    assert "35.800" in response.text
    assert "Normal · 0 °C" in response.text


def test_calculator_form_accepts_multiple_component_rows():
    payload = urlencode(
        [
            ("component_name", "methane"),
            ("component_name", "ethane"),
            ("component_mole_percentage", "50"),
            ("component_mole_percentage", "50"),
            ("component_weight_percentage", "34.7864"),
            ("component_weight_percentage", "65.2136"),
            ("basis", "Mol %"),
            ("temperature_c", "25"),
            ("pressure_atm", "1"),
            ("model", "PRMIX"),
        ]
    )

    response = build_client().post(
        "/calculator",
        content=payload,
        headers={"content-type": "application/x-www-form-urlencoded"},
    )

    assert response.status_code == 200
    assert "Methane" in response.text
    assert "Ethane" in response.text
    assert len(re.findall(r"50\.000\s+%", response.text)) == 2


def test_calculator_form_uses_weight_basis_and_repopulates_both_columns():
    payload = urlencode(
        [
            ("component_name", "methane"),
            ("component_name", "ethane"),
            ("component_mole_percentage", ""),
            ("component_mole_percentage", ""),
            ("component_weight_percentage", "50"),
            ("component_weight_percentage", "50"),
            ("basis", "Wt %"),
            ("temperature_c", "25"),
            ("pressure_atm", "1"),
            ("model", "PRMIX"),
        ]
    )

    response = build_client().post(
        "/calculator",
        content=payload,
        headers={"content-type": "application/x-www-form-urlencoded"},
    )

    assert response.status_code == 200
    assert re.search(
        r'name="basis"\s+value="Wt %"\s+checked',
        response.text,
    )
    assert 'name="component_mole_percentage"' in response.text
    assert 'name="component_weight_percentage"' in response.text
    assert "Calculated properties" in response.text


def test_calculator_form_preserves_input_and_renders_validation_error():
    response = build_client().post(
        "/calculator",
        data=methane_form(component_mole_percentage="90"),
    )

    assert response.status_code == 422
    assert "Calculation input needs attention" in response.text
    assert "total must be 100" in response.text
    assert 'value="90"' in response.text
    assert re.search(
        r'<option\s+value="methane"\s+selected',
        response.text,
    )


def test_calculator_form_can_prepare_report_projection():
    response = build_client().post(
        "/calculator",
        data=methane_form(include_report_projection="on"),
    )

    assert response.status_code == 200
    assert "Report data prepared" in response.text
    assert "result rows are" in response.text
