"""Typed requests and responses used by application workflows."""

from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Any

from thermo_components.domain.results import extract_scalar_density_value


DensityValue = float | tuple[float | None, float | None] | None


@dataclass(frozen=True)
class PropertyCalculationRequest:
    component_names: tuple[str, ...]
    mole_percents: tuple[float, ...]
    weight_percents: tuple[float, ...]
    basis: str
    temperature_k: float
    pressure_pa: float
    pressure_atm: float

    @classmethod
    def from_sequences(
        cls,
        component_names,
        mole_percents,
        weight_percents,
        basis: str,
        temperature_k: float,
        pressure_pa: float,
        pressure_atm: float,
    ) -> "PropertyCalculationRequest":
        return cls(
            component_names=tuple(component_names),
            mole_percents=tuple(float(value) for value in mole_percents),
            weight_percents=tuple(float(value) for value in weight_percents),
            basis=basis,
            temperature_k=float(temperature_k),
            pressure_pa=float(pressure_pa),
            pressure_atm=float(pressure_atm),
        )


@dataclass(frozen=True)
class DensityCalculationResult:
    value: DensityValue
    phase: str | None
    error: str | None
    scalar_kg_m3: float | None

    @classmethod
    def from_raw(
        cls,
        value: DensityValue,
        phase: str | None,
        error: str | None,
    ) -> "DensityCalculationResult":
        return cls(
            value=value,
            phase=phase,
            error=error,
            scalar_kg_m3=extract_scalar_density_value(value, phase, error),
        )


@dataclass(frozen=True)
class PropertyCalculationResponse:
    average_molecular_weight: float
    component_names: tuple[str, ...]
    mole_percents: tuple[float, ...]
    weight_percents: tuple[float, ...]
    thermo_route: str
    model_display: str
    selected_density: DensityCalculationResult
    normal_density: DensityCalculationResult
    standard_density: DensityCalculationResult
    bubble_point_c: float | None
    bubble_point_error: str | None
    mixture_lhv_mj_nm3: float
    missing_lhv: tuple[str, ...]
    basis: str
    eos: str
    pressure_atm: float | None
    warnings: tuple[str, ...] = field(default_factory=tuple)

    @classmethod
    def from_mapping(
        cls,
        data: Mapping[str, Any],
    ) -> "PropertyCalculationResponse":
        normal_value = data.get(
            "density_normal_result",
            data.get("density_normal_kg_m3"),
        )
        normal_phase = data.get("density_normal_phase")
        normal_error = data.get("density_normal_error")
        normal_scalar = data.get("density_normal_kg_m3")
        if "density_normal_kg_m3" not in data:
            normal_scalar = extract_scalar_density_value(
                normal_value,
                normal_phase,
                normal_error,
            )

        standard_value = data.get(
            "density_standard_result",
            data.get("density_standard_kg_m3"),
        )
        standard_phase = data.get("density_standard_phase")
        standard_error = data.get("density_standard_error")
        standard_scalar = data.get("density_standard_kg_m3")
        if "density_standard_kg_m3" not in data:
            standard_scalar = extract_scalar_density_value(
                standard_value,
                standard_phase,
                standard_error,
            )

        pressure_atm = data.get("pressure_atm")
        return cls(
            average_molecular_weight=float(data.get("mw") or 0.0),
            component_names=tuple(data.get("comp_names") or ()),
            mole_percents=tuple(data.get("mol_percents") or ()),
            weight_percents=tuple(data.get("wt_percents") or ()),
            thermo_route=str(data.get("thermo_route") or ""),
            model_display=str(data.get("model_display") or data.get("eos") or ""),
            selected_density=DensityCalculationResult.from_raw(
                data.get("density_result"),
                data.get("phase"),
                data.get("density_error"),
            ),
            normal_density=DensityCalculationResult(
                value=normal_value,
                phase=normal_phase,
                error=normal_error,
                scalar_kg_m3=normal_scalar,
            ),
            standard_density=DensityCalculationResult(
                value=standard_value,
                phase=standard_phase,
                error=standard_error,
                scalar_kg_m3=standard_scalar,
            ),
            bubble_point_c=data.get("bubble_point"),
            bubble_point_error=data.get("bp_error"),
            mixture_lhv_mj_nm3=float(data.get("mixture_lhv") or 0.0),
            missing_lhv=tuple(data.get("missing_lhv") or ()),
            basis=str(data.get("basis") or ""),
            eos=str(data.get("eos") or ""),
            pressure_atm=(
                float(pressure_atm)
                if pressure_atm is not None
                else None
            ),
            warnings=tuple(data.get("warnings") or ()),
        )

    def to_legacy_dict(self) -> dict[str, Any]:
        """Expose the previous result shape during the UI migration."""
        return {
            "mw": self.average_molecular_weight,
            "comp_names": list(self.component_names),
            "mol_percents": list(self.mole_percents),
            "wt_percents": list(self.weight_percents),
            "thermo_route": self.thermo_route,
            "model_display": self.model_display,
            "density_result": self.selected_density.value,
            "phase": self.selected_density.phase,
            "density_error": self.selected_density.error,
            "density_actual_kg_m3": self.selected_density.scalar_kg_m3,
            "density_normal_result": self.normal_density.value,
            "density_normal_phase": self.normal_density.phase,
            "density_normal_error": self.normal_density.error,
            "density_normal_kg_m3": self.normal_density.scalar_kg_m3,
            "density_standard_result": self.standard_density.value,
            "density_standard_phase": self.standard_density.phase,
            "density_standard_error": self.standard_density.error,
            "density_standard_kg_m3": self.standard_density.scalar_kg_m3,
            "bubble_point": self.bubble_point_c,
            "bp_error": self.bubble_point_error,
            "mixture_lhv": self.mixture_lhv_mj_nm3,
            "missing_lhv": list(self.missing_lhv),
            "basis": self.basis,
            "eos": self.eos,
            "pressure_atm": self.pressure_atm,
            "warnings": list(self.warnings),
        }


def coerce_property_response(
    value: PropertyCalculationResponse | Mapping[str, Any],
) -> PropertyCalculationResponse:
    if isinstance(value, PropertyCalculationResponse):
        return value
    return PropertyCalculationResponse.from_mapping(value)


@dataclass(frozen=True)
class FlowConversionRequest:
    input_text: str
    from_unit: str
    to_unit: str
    normal_density_kg_m3: float | None = None
    standard_density_kg_m3: float | None = None


@dataclass(frozen=True)
class FlowConversionResponse:
    value: float | None
    display_value: str


@dataclass(frozen=True)
class NormalizeCompositionRequest:
    percentages: tuple[float, ...]
    precision: int = 4


@dataclass(frozen=True)
class NormalizeCompositionResponse:
    percentages: tuple[float, ...]


@dataclass(frozen=True)
class DeriveCompositionRequest:
    component_names: tuple[str, ...]
    active_percentages: tuple[float, ...]
    basis: str


@dataclass(frozen=True)
class DeriveCompositionResponse:
    percentages: tuple[float | None, ...]


@dataclass(frozen=True)
class ReportResultRow:
    property_name: str
    value: Any
    unit: str = ""
    notes: str = ""

    def to_dict(self) -> dict[str, Any]:
        return {
            "Property": self.property_name,
            "Value": self.value,
            "Unit": self.unit,
            "Notes": self.notes,
        }


@dataclass(frozen=True)
class ReportWarningRow:
    warning_type: str
    details: str

    def to_dict(self) -> dict[str, str]:
        return {
            "Warning Type": self.warning_type,
            "Details": self.details,
        }


@dataclass(frozen=True)
class ReportPreparationRequest:
    calculation: PropertyCalculationResponse
    lhv_data_available: bool


@dataclass(frozen=True)
class ReportProjection:
    result_rows: tuple[ReportResultRow, ...]
    warning_rows: tuple[ReportWarningRow, ...]
