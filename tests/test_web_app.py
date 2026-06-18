from fastapi import FastAPI
from fastapi.testclient import TestClient

from thermo_components.adapters.web.app import app, create_app


def test_create_app_returns_configured_fastapi_app():
    web_app = create_app()

    assert isinstance(web_app, FastAPI)
    assert web_app.title == "Thermo Components"
    assert web_app.version == "1.0.0"


def test_health_endpoint_reports_service_status():
    client = TestClient(create_app())

    response = client.get("/health")

    assert response.status_code == 200
    assert response.json() == {
        "service": "thermo-components",
        "status": "ok",
    }


def test_module_app_is_ready_for_uvicorn_import():
    with TestClient(app) as client:
        response = client.get("/health")

    assert response.status_code == 200


def test_openapi_document_includes_calculation_endpoint():
    response = TestClient(create_app()).get("/openapi.json")

    assert response.status_code == 200
    openapi = response.json()
    assert "/api/calculations" in openapi["paths"]
    assert "/api/flow-conversions" in openapi["paths"]
    assert openapi["components"]["schemas"]["CalculationRequestSchema"][
        "examples"
    ] == [
        {
            "components": [
                {
                    "name": "methane",
                    "percentage": 100.0,
                }
            ],
            "basis": "Mol %",
            "temperature_c": 25.0,
            "pressure_atm": 1.0,
            "model": "PRMIX",
            "include_report_projection": False,
        }
    ]
