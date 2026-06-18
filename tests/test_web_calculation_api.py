import pytest
from fastapi.testclient import TestClient

from thermo_components.adapters.web.app import create_app
from thermo_components.application.dto import (
    DensityCalculationResult,
    DeriveCompositionResponse,
    PropertyCalculationResponse,
    ReportProjection,
    ReportResultRow,
)
from thermo_components.application.services import (
    CalculationSessionResponse,
    FlowDensityState,
)
from thermo_components.bootstrap.web import (
    WebDependencies,
    build_web_dependencies,
)


class FakeDeriveCompositionUseCase:
    def __init__(self, percentages=(100.0,)):
        self.percentages = percentages
        self.requests = []

    def execute(self, request):
        self.requests.append(request)
        return DeriveCompositionResponse(percentages=self.percentages)


class FakeCalculationSession:
    def __init__(self, response=None, error=None):
        self.response = response
        self.error = error
        self.requests = []

    def execute(self, request):
        self.requests.append(request)
        if self.error is not None:
            raise self.error
        return self.response


class RecordingSessionFactory:
    def __init__(self, session):
        self.session = session
        self.models = []

    def __call__(self, model):
        self.models.append(model)
        return self.session


def methane_session_response(
    *,
    lhv_data_available=True,
) -> CalculationSessionResponse:
    calculation = PropertyCalculationResponse(
        average_molecular_weight=16.04,
        component_names=("methane",),
        mole_percents=(100.0,),
        weight_percents=(100.0,),
        thermo_route="prmix_default",
        model_display="PRMIX",
        selected_density=DensityCalculationResult.from_raw(
            0.657,
            "Vapor",
            None,
        ),
        normal_density=DensityCalculationResult.from_raw(
            0.718,
            "Vapor",
            None,
        ),
        standard_density=DensityCalculationResult.from_raw(
            0.679,
            "Vapor",
            None,
        ),
        bubble_point_c=-161.57,
        bubble_point_error=None,
        mixture_lhv_mj_nm3=35.8,
        missing_lhv=(),
        basis="Mol %",
        eos="PRMIX",
        pressure_atm=1.0,
        warnings=(),
    )
    return CalculationSessionResponse(
        calculation=calculation,
        flow_densities=FlowDensityState(
            selected_density_kg_m3=0.657,
            normal_density_kg_m3=0.718,
            standard_density_kg_m3=0.679,
        ),
        warning_messages=(),
        lhv_data_available=lhv_data_available,
        report_projection=ReportProjection(
            result_rows=(
                ReportResultRow(
                    "Average molecular weight",
                    16.04,
                    "g/mol",
                ),
            ),
            warning_rows=(),
        ),
    )


def valid_methane_payload(**overrides):
    payload = {
        "components": [
            {
                "name": " Methane ",
                "percentage": 100.0,
            }
        ],
        "basis": "Mol %",
        "temperature_c": 25.0,
        "pressure_atm": 1.0,
        "model": "PRMIX",
        "include_report_projection": True,
    }
    payload.update(overrides)
    return payload


def test_calculation_endpoint_maps_web_input_to_application_request():
    derive_use_case = FakeDeriveCompositionUseCase()
    session = FakeCalculationSession(methane_session_response())
    session_factory = RecordingSessionFactory(session)
    dependencies = WebDependencies(
        lhv_database={"methane": 35.8},
        derive_composition_use_case=derive_use_case,
        calculation_session_factory=session_factory,
    )

    response = TestClient(create_app(dependencies)).post(
        "/api/calculations",
        json=valid_methane_payload(),
    )

    assert response.status_code == 200
    assert response.json()["calculation"]["average_molecular_weight"] == 16.04
    assert response.json()["flow_densities"] == {
        "selected_density_kg_m3": 0.657,
        "normal_density_kg_m3": 0.718,
        "standard_density_kg_m3": 0.679,
    }
    assert len(response.json()["lhv"]["volumetric"]) == 5
    assert len(response.json()["lhv"]["mass_basis"]) == 10
    assert response.json()["lhv"]["volumetric"][0] == {
        "unit": "MJ/Nm\u00b3",
        "value": 35.8,
        "display_value": "35.80",
    }
    assert response.json()["lhv"]["mass_basis"][0]["unit"] == "MJ/kg"
    assert response.json()["lhv"]["mass_basis"][0]["value"] == (
        pytest.approx(35.8 / (16.04 / 22.414))
    )
    assert response.json()["report_projection"]["result_rows"][0] == {
        "property_name": "Average molecular weight",
        "value": 16.04,
        "unit": "g/mol",
        "notes": "",
    }
    assert session_factory.models == ["PRMIX"]
    assert derive_use_case.requests[0].component_names == ("methane",)
    assert derive_use_case.requests[0].active_percentages == (100.0,)

    session_request = session.requests[0]
    calculation_request = session_request.calculation_request
    assert calculation_request.component_names == ("methane",)
    assert calculation_request.mole_percents == (100.0,)
    assert calculation_request.weight_percents == (100.0,)
    assert calculation_request.temperature_k == pytest.approx(298.15)
    assert calculation_request.pressure_pa == pytest.approx(101325.0)
    assert calculation_request.pressure_atm == pytest.approx(1.0)
    assert session_request.lhv_data_available is True
    assert session_request.include_report_projection is True


def test_calculation_endpoint_derives_mole_percentages_for_weight_basis():
    derive_use_case = FakeDeriveCompositionUseCase(
        percentages=(70.0, 30.0),
    )
    session = FakeCalculationSession(methane_session_response())
    dependencies = WebDependencies(
        lhv_database={},
        derive_composition_use_case=derive_use_case,
        calculation_session_factory=RecordingSessionFactory(session),
    )
    payload = valid_methane_payload(
        components=[
            {"name": "methane", "percentage": 40.0},
            {"name": "ethane", "percentage": 60.0},
        ],
        basis="Wt %",
        include_report_projection=False,
    )

    response = TestClient(create_app(dependencies)).post(
        "/api/calculations",
        json=payload,
    )

    assert response.status_code == 200
    calculation_request = session.requests[0].calculation_request
    assert calculation_request.component_names == ("methane", "ethane")
    assert calculation_request.mole_percents == (70.0, 30.0)
    assert calculation_request.weight_percents == (40.0, 60.0)
    assert calculation_request.basis == "Wt %"


@pytest.mark.parametrize(
    ("payload_override", "message"),
    [
        (
            {
                "components": [
                    {"name": "methane", "percentage": 90.0},
                ]
            },
            "total must be 100",
        ),
        (
            {
                "components": [
                    {"name": "unobtainium", "percentage": 100.0},
                ]
            },
            "Unsupported component",
        ),
        (
            {
                "components": [
                    {"name": "methane", "percentage": 50.0},
                    {"name": "methane", "percentage": 50.0},
                ]
            },
            "Duplicate components",
        ),
    ],
)
def test_calculation_endpoint_rejects_invalid_web_input(
    payload_override,
    message,
):
    session_factory = RecordingSessionFactory(
        FakeCalculationSession(methane_session_response())
    )
    dependencies = WebDependencies(
        lhv_database={},
        derive_composition_use_case=FakeDeriveCompositionUseCase(),
        calculation_session_factory=session_factory,
    )

    response = TestClient(create_app(dependencies)).post(
        "/api/calculations",
        json=valid_methane_payload(**payload_override),
    )

    assert response.status_code == 422
    assert message in response.text
    assert session_factory.models == []


def test_calculation_endpoint_maps_application_validation_error_to_422():
    session = FakeCalculationSession(
        error=ValueError("Calculation request was rejected.")
    )
    dependencies = WebDependencies(
        lhv_database={},
        derive_composition_use_case=FakeDeriveCompositionUseCase(),
        calculation_session_factory=RecordingSessionFactory(session),
    )

    response = TestClient(create_app(dependencies)).post(
        "/api/calculations",
        json=valid_methane_payload(),
    )

    assert response.status_code == 422
    assert response.json() == {
        "detail": "Calculation request was rejected.",
    }


def test_calculation_endpoint_marks_lhv_values_unavailable_without_database():
    dependencies = WebDependencies(
        lhv_database={},
        derive_composition_use_case=FakeDeriveCompositionUseCase(),
        calculation_session_factory=RecordingSessionFactory(
            FakeCalculationSession(
                methane_session_response(lhv_data_available=False)
            )
        ),
    )

    response = TestClient(create_app(dependencies)).post(
        "/api/calculations",
        json=valid_methane_payload(),
    )

    assert response.status_code == 200
    lhv = response.json()["lhv"]
    assert lhv["available"] is False
    assert all(row["value"] is None for row in lhv["volumetric"])
    assert all(row["display_value"] == "N/A" for row in lhv["mass_basis"])


def test_calculation_endpoint_runs_real_methane_calculation():
    dependencies = build_web_dependencies(
        lhv_data={"methane": 35.8},
    )

    response = TestClient(create_app(dependencies)).post(
        "/api/calculations",
        json=valid_methane_payload(
            include_report_projection=False,
        ),
    )

    assert response.status_code == 200
    result = response.json()
    assert result["calculation"]["average_molecular_weight"] == pytest.approx(
        16.04
    )
    assert result["calculation"]["selected_density"]["scalar_kg_m3"] == (
        pytest.approx(0.657, rel=0.01)
    )
    assert result["calculation"]["mixture_lhv_mj_nm3"] == pytest.approx(35.8)
    assert result["flow_densities"]["normal_density_kg_m3"] is not None
    assert result["report_projection"] is None
