"""HTTP request and response schemas for the web adapter."""

from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

from thermo_components.domain.composition import (
    MOLECULAR_WEIGHTS,
    normalize_component_identity,
)


CompositionBasisValue = Literal["Mol %", "Wt %"]
CalculationModel = Literal["PRMIX"]
DensityValue = float | tuple[float | None, float | None] | None


class WebSchema(BaseModel):
    model_config = ConfigDict(
        allow_inf_nan=False,
        extra="forbid",
    )


class ComponentInputSchema(WebSchema):
    name: str = Field(min_length=1)
    percentage: float = Field(ge=0.0, le=100.0)

    @field_validator("name")
    @classmethod
    def validate_component_name(cls, value: str) -> str:
        normalized = normalize_component_identity(value)
        if normalized not in MOLECULAR_WEIGHTS or not normalized:
            raise ValueError(f"Unsupported component: {value}")
        return normalized


class CalculationRequestSchema(WebSchema):
    components: list[ComponentInputSchema] = Field(min_length=1)
    basis: CompositionBasisValue
    temperature_c: float = Field(gt=-273.15)
    pressure_atm: float = Field(gt=0.0)
    model: CalculationModel = "PRMIX"
    include_report_projection: bool = False

    @model_validator(mode="after")
    def validate_composition(self) -> "CalculationRequestSchema":
        component_names = [component.name for component in self.components]
        if len(set(component_names)) != len(component_names):
            raise ValueError("Duplicate components are not allowed.")

        total = sum(component.percentage for component in self.components)
        if abs(total - 100.0) > 1e-4:
            raise ValueError(f"Invalid {self.basis} input: total must be 100.")
        return self


class DensityResultSchema(WebSchema):
    value: DensityValue
    phase: str | None
    error: str | None
    scalar_kg_m3: float | None


class PropertyCalculationSchema(WebSchema):
    average_molecular_weight: float
    component_names: list[str]
    mole_percents: list[float]
    weight_percents: list[float]
    thermo_route: str
    model_display: str
    selected_density: DensityResultSchema
    normal_density: DensityResultSchema
    standard_density: DensityResultSchema
    bubble_point_c: float | None
    bubble_point_error: str | None
    mixture_lhv_mj_nm3: float
    missing_lhv: list[str]
    basis: CompositionBasisValue
    eos: str
    pressure_atm: float | None
    warnings: list[str]


class FlowDensityStateSchema(WebSchema):
    selected_density_kg_m3: float | None
    normal_density_kg_m3: float | None
    standard_density_kg_m3: float | None


class ReportResultRowSchema(WebSchema):
    property_name: str
    value: Any
    unit: str
    notes: str


class ReportWarningRowSchema(WebSchema):
    warning_type: str
    details: str


class ReportProjectionSchema(WebSchema):
    result_rows: list[ReportResultRowSchema]
    warning_rows: list[ReportWarningRowSchema]


class CalculationResponseSchema(WebSchema):
    calculation: PropertyCalculationSchema
    flow_densities: FlowDensityStateSchema
    warnings: list[str]
    report_projection: ReportProjectionSchema | None
