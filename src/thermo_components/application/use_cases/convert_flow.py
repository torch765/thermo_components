"""Flow conversion application workflow."""

from thermo_components.application.dto import (
    FlowConversionRequest,
    FlowConversionResponse,
)
from thermo_components.domain.flow_conversion import (
    convert_flow,
    format_flow_value,
    parse_flow_input,
)


class ConvertFlowUseCase:
    def execute(
        self,
        request: FlowConversionRequest,
    ) -> FlowConversionResponse:
        if not request.from_unit or not request.to_unit:
            raise ValueError("Select source and destination units")

        try:
            value = parse_flow_input(request.input_text)
        except ValueError as exc:
            raise ValueError("Invalid input") from exc

        if value is None:
            return FlowConversionResponse(value=None, display_value="")

        converted_value = convert_flow(
            value,
            request.from_unit,
            request.to_unit,
            density_normal_kg_m3=request.normal_density_kg_m3,
            density_standard_kg_m3=request.standard_density_kg_m3,
        )
        return FlowConversionResponse(
            value=converted_value,
            display_value=format_flow_value(converted_value),
        )
