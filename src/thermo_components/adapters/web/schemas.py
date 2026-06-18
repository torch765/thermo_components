"""HTTP request and response schemas for the web adapter."""

from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

from thermo_components.domain.composition import (
    MOLECULAR_WEIGHTS,
    normalize_component_identity,
)
from thermo_components.domain.flow_units import FLOW_UNITS


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
    model_config = ConfigDict(
        json_schema_extra={
            "examples": [
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
        }
    )

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


class CompositionDeriveRequestSchema(WebSchema):
    components: list[ComponentInputSchema] = Field(min_length=1)
    basis: CompositionBasisValue

    @model_validator(mode="after")
    def validate_unique_components(self) -> "CompositionDeriveRequestSchema":
        component_names = [component.name for component in self.components]
        if len(set(component_names)) != len(component_names):
            raise ValueError("Duplicate components are not allowed.")
        return self


class CompositionDeriveResponseSchema(WebSchema):
    percentages: list[float | None]


class CompositionNormalizeRequestSchema(WebSchema):
    percentages: list[float] = Field(min_length=1)
    precision: int = Field(default=4, ge=0, le=8)

    @field_validator("percentages")
    @classmethod
    def validate_non_negative_percentages(
        cls,
        values: list[float],
    ) -> list[float]:
        if any(value < 0 for value in values):
            raise ValueError("Percentages cannot be negative.")
        return values


class CompositionNormalizeResponseSchema(WebSchema):
    percentages: list[float]


class FlowConversionRequestSchema(WebSchema):
    input_text: str
    from_unit: str
    to_unit: str
    normal_density_kg_m3: float | None = Field(default=None, gt=0.0)
    standard_density_kg_m3: float | None = Field(default=None, gt=0.0)

    @field_validator("from_unit", "to_unit")
    @classmethod
    def validate_flow_unit(cls, value: str) -> str:
        if value not in FLOW_UNITS:
            raise ValueError(f"Unsupported flow unit: {value}")
        return value


class FlowConversionResponseSchema(WebSchema):
    value: float | None
    display_value: str


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


class LhvValueSchema(WebSchema):
    unit: str
    value: float | None
    display_value: str


class LhvDisplaySchema(WebSchema):
    available: bool
    volumetric: list[LhvValueSchema]
    mass_basis: list[LhvValueSchema]


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
    lhv: LhvDisplaySchema
    warnings: list[str]
    report_projection: ReportProjectionSchema | None
