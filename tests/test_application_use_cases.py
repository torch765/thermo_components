import pytest

from thermo_components.application.dto import (
    DeriveCompositionRequest,
    FlowConversionRequest,
    NormalizeCompositionRequest,
    PropertyCalculationRequest,
    PropertyCalculationResponse,
    ReportPreparationRequest,
)
from thermo_components.application.use_cases import (
    CalculatePropertiesUseCase,
    ConvertFlowUseCase,
    DeriveCompositionUseCase,
    NormalizeCompositionUseCase,
    PrepareReportUseCase,
)
from thermo_components.domain.conditions import (
    NORMAL_CONDITION,
    STANDARD_CONDITION,
)
from thermo_components.domain.thermo_routes import (
    PRMIX_DEFAULT_ROUTE,
    PURE_WATER_ROUTE,
)


class FakePropertyCalculator:
    eos = "PRMIX"

    def __init__(self):
        self.components = {}
        self.density_calls = []
        self.bubble_point_calls = []

    def set_components(self, components):
        self.components = dict(components)

    def calculate_density_for_route(
        self,
        temperature_k,
        pressure_pa,
        route_id,
    ):
        self.density_calls.append((temperature_k, pressure_pa, route_id))
        return 1.25, "Vapor", None

    def calculate_bubble_point_for_route(self, pressure_pa, route_id):
        self.bubble_point_calls.append((pressure_pa, route_id))
        return -161.5, None


def methane_request() -> PropertyCalculationRequest:
    return PropertyCalculationRequest.from_sequences(
        component_names=["methane"],
        mole_percents=[100.0],
        weight_percents=[100.0],
        basis="Mol %",
        temperature_k=298.15,
        pressure_pa=101325.0,
        pressure_atm=1.0,
    )


def test_calculate_properties_orchestrates_all_reference_calculations():
    calculator = FakePropertyCalculator()
    use_case = CalculatePropertiesUseCase(
        calculator,
        {"methane": 35.8},
    )

    response = use_case.execute(methane_request())

    assert isinstance(response, PropertyCalculationResponse)
    assert response.average_molecular_weight == pytest.approx(16.04)
    assert response.mixture_lhv_mj_nm3 == pytest.approx(35.8)
    assert response.thermo_route == PRMIX_DEFAULT_ROUTE
    assert response.selected_density.scalar_kg_m3 == pytest.approx(1.25)
    assert response.normal_density.scalar_kg_m3 == pytest.approx(1.25)
    assert response.standard_density.scalar_kg_m3 == pytest.approx(1.25)
    assert calculator.components == {"methane": 1.0}
    assert calculator.density_calls == [
        (298.15, 101325.0, PRMIX_DEFAULT_ROUTE),
        (
            NORMAL_CONDITION.temperature_k,
            NORMAL_CONDITION.pressure_pa,
            PRMIX_DEFAULT_ROUTE,
        ),
        (
            STANDARD_CONDITION.temperature_k,
            STANDARD_CONDITION.pressure_pa,
            PRMIX_DEFAULT_ROUTE,
        ),
    ]
    assert calculator.bubble_point_calls == [
        (101325.0, PRMIX_DEFAULT_ROUTE)
    ]


def test_calculate_properties_selects_pure_water_route():
    calculator = FakePropertyCalculator()
    request = PropertyCalculationRequest.from_sequences(
        component_names=["water"],
        mole_percents=[100.0],
        weight_percents=[100.0],
        basis="Mol %",
        temperature_k=298.15,
        pressure_pa=101325.0,
        pressure_atm=1.0,
    )

    response = CalculatePropertiesUseCase(calculator, {"water": 0.0}).execute(
        request
    )

    assert response.thermo_route == PURE_WATER_ROUTE
    assert response.warnings == ()
    assert all(call[2] == PURE_WATER_ROUTE for call in calculator.density_calls)


def test_calculate_properties_rejects_invalid_active_total():
    calculator = FakePropertyCalculator()
    request = PropertyCalculationRequest.from_sequences(
        component_names=["methane"],
        mole_percents=[50.0],
        weight_percents=[100.0],
        basis="Mol %",
        temperature_k=298.15,
        pressure_pa=101325.0,
        pressure_atm=1.0,
    )

    with pytest.raises(ValueError, match="Invalid Mol % input"):
        CalculatePropertiesUseCase(calculator, {}).execute(request)


def test_flow_and_composition_use_cases_return_typed_responses():
    flow_response = ConvertFlowUseCase().execute(
        FlowConversionRequest(
            input_text="1,000",
            from_unit="kg/h",
            to_unit="t/d",
        )
    )
    normalized = NormalizeCompositionUseCase().execute(
        NormalizeCompositionRequest(percentages=(40.0, 40.0))
    )
    derived = DeriveCompositionUseCase().execute(
        DeriveCompositionRequest(
            component_names=("methane", "ethane"),
            active_percentages=(50.0, 50.0),
            basis="Mol %",
        )
    )

    assert flow_response.value == pytest.approx(24.0)
    assert flow_response.display_value == "24"
    assert normalized.percentages == (50.0, 50.0)
    assert sum(value for value in derived.percentages if value is not None) == (
        pytest.approx(100.0)
    )


def test_prepare_report_projects_typed_calculation_response():
    calculation = CalculatePropertiesUseCase(
        FakePropertyCalculator(),
        {"methane": 35.8},
    ).execute(methane_request())

    projection = PrepareReportUseCase().execute(
        ReportPreparationRequest(
            calculation=calculation,
            lhv_data_available=True,
        )
    )

    result_rows = [row.to_dict() for row in projection.result_rows]
    assert result_rows[0] == {
        "Property": "Average molecular weight",
        "Value": pytest.approx(16.04),
        "Unit": "g/mol",
        "Notes": "",
    }
    assert len(
        [
            row
            for row in result_rows
            if row["Property"] == "Mixture LHV (Mass basis)"
        ]
    ) == 10
