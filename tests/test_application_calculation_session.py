import pytest

from thermo_components.application.dto import (
    DensityCalculationResult,
    PropertyCalculationRequest,
    PropertyCalculationResponse,
    ReportProjection,
    ReportResultRow,
    ReportWarningRow,
)
from thermo_components.application.services import (
    CalculationSessionRequest,
    CalculationSessionResponse,
    CalculationSessionService,
    FlowDensityState,
)


class FakeCalculatePropertiesUseCase:
    def __init__(self, response=None, error=None):
        self.response = response
        self.error = error
        self.requests = []

    def execute(self, request):
        self.requests.append(request)
        if self.error is not None:
            raise self.error
        return self.response


class FakePrepareReportUseCase:
    def __init__(self, projection=None):
        self.projection = projection or ReportProjection(
            result_rows=(
                ReportResultRow("Average molecular weight", 16.04, "g/mol"),
            ),
            warning_rows=(
                ReportWarningRow("Thermo", "Check reference density"),
            ),
        )
        self.requests = []

    def execute(self, request):
        self.requests.append(request)
        return self.projection


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


def methane_response() -> PropertyCalculationResponse:
    return PropertyCalculationResponse(
        average_molecular_weight=16.04,
        component_names=("methane",),
        mole_percents=(100.0,),
        weight_percents=(100.0,),
        thermo_route="prmix",
        model_display="PRMIX",
        selected_density=DensityCalculationResult.from_raw(
            0.657,
            "Vapor",
            None,
        ),
        normal_density=DensityCalculationResult.from_raw(
            0.716,
            "Vapor",
            None,
        ),
        standard_density=DensityCalculationResult.from_raw(
            0.68,
            "Vapor",
            None,
        ),
        bubble_point_c=-161.5,
        bubble_point_error=None,
        mixture_lhv_mj_nm3=35.8,
        missing_lhv=(),
        basis="Mol %",
        eos="PRMIX",
        pressure_atm=1.0,
        warnings=("Check reference density",),
    )


def test_calculation_session_composes_calculation_and_report_projection():
    calculation_response = methane_response()
    calculate_use_case = FakeCalculatePropertiesUseCase(calculation_response)
    prepare_report_use_case = FakePrepareReportUseCase()
    service = CalculationSessionService(
        calculate_use_case,
        prepare_report_use_case,
    )
    request = CalculationSessionRequest(
        calculation_request=methane_request(),
        lhv_data_available=True,
    )

    response = service.execute(request)

    assert isinstance(response, CalculationSessionResponse)
    assert response.calculation is calculation_response
    assert response.flow_densities == FlowDensityState(
        selected_density_kg_m3=0.657,
        normal_density_kg_m3=0.716,
        standard_density_kg_m3=0.68,
    )
    assert response.warning_messages == ("Check reference density",)
    assert response.report_projection is prepare_report_use_case.projection
    assert calculate_use_case.requests == [request.calculation_request]
    assert prepare_report_use_case.requests[0].calculation is (
        calculation_response
    )
    assert prepare_report_use_case.requests[0].lhv_data_available is True


def test_calculation_session_can_skip_report_projection():
    calculate_use_case = FakeCalculatePropertiesUseCase(methane_response())
    prepare_report_use_case = FakePrepareReportUseCase()
    service = CalculationSessionService(
        calculate_use_case,
        prepare_report_use_case,
    )

    response = service.execute(
        CalculationSessionRequest(
            calculation_request=methane_request(),
            lhv_data_available=True,
            include_report_projection=False,
        )
    )

    assert response.report_projection is None
    assert prepare_report_use_case.requests == []


def test_calculation_session_propagates_validation_errors():
    calculate_use_case = FakeCalculatePropertiesUseCase(
        error=ValueError("Invalid Mol % input.")
    )
    prepare_report_use_case = FakePrepareReportUseCase()
    service = CalculationSessionService(
        calculate_use_case,
        prepare_report_use_case,
    )

    with pytest.raises(ValueError, match="Invalid Mol % input"):
        service.execute(
            CalculationSessionRequest(
                calculation_request=methane_request(),
                lhv_data_available=True,
            )
        )

    assert prepare_report_use_case.requests == []
