import pytest
from fastapi.testclient import TestClient

from thermo_components.adapters.web.app import create_app
from thermo_components.bootstrap.web import build_web_dependencies


@pytest.fixture
def client():
    return TestClient(create_app(build_web_dependencies(lhv_data={})))


def test_derive_composition_returns_inactive_weight_basis(client):
    response = client.post(
        "/api/compositions/derive",
        json={
            "basis": "Mol %",
            "components": [
                {"name": "methane", "percentage": 50.0},
                {"name": "ethane", "percentage": 50.0},
            ],
        },
    )

    assert response.status_code == 200
    assert response.json()["percentages"] == pytest.approx(
        [
            16.04 / (16.04 + 30.07) * 100.0,
            30.07 / (16.04 + 30.07) * 100.0,
        ]
    )


def test_derive_composition_returns_inactive_mole_basis(client):
    response = client.post(
        "/api/compositions/derive",
        json={
            "basis": "Wt %",
            "components": [
                {"name": "methane", "percentage": 50.0},
                {"name": "ethane", "percentage": 50.0},
            ],
        },
    )

    assert response.status_code == 200
    assert sum(response.json()["percentages"]) == pytest.approx(100.0)
    assert (
        response.json()["percentages"][0]
        > response.json()["percentages"][1]
    )


def test_normalize_composition_scales_values_to_exactly_100(client):
    response = client.post(
        "/api/compositions/normalize",
        json={
            "percentages": [40.0, 40.0],
            "precision": 4,
        },
    )

    assert response.status_code == 200
    assert response.json() == {"percentages": [50.0, 50.0]}


def test_normalize_composition_rejects_zero_total(client):
    response = client.post(
        "/api/compositions/normalize",
        json={"percentages": [0.0, 0.0]},
    )

    assert response.status_code == 422
    assert response.json() == {
        "detail": "Cannot normalize: total is zero.",
    }


def test_composition_endpoints_are_in_openapi(client):
    openapi = client.get("/openapi.json").json()

    assert "/api/compositions/derive" in openapi["paths"]
    assert "/api/compositions/normalize" in openapi["paths"]
