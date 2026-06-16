"""Presentation helpers for the Qt results list."""

from thermo_components.application.dto import coerce_property_response
from thermo_components.domain.lhv import (
    build_lhv_display_values,
    format_lhv_display_value,
)


def build_result_list_items(
    calculation,
    *,
    lhv_data_loaded: bool,
    fallback_model_display: str = "",
) -> tuple[str, ...]:
    """Build the text rows rendered in the main results list."""
    response = coerce_property_response(calculation)
    result_data = response.to_legacy_dict()
    items: list[str] = [
        f"Calculating for {result_data['basis']} basis...",
        "Model / EOS: "
        f"{result_data.get('model_display', fallback_model_display)}",
        f"Avg. Molecular Wt: {result_data['mw']:.2f} g/mol",
    ]

    density_result = result_data["density_result"]
    phase = result_data["phase"]
    error = result_data["density_error"]
    if error:
        items.extend(
            [
                "Phase @ selected conditions: N.A.",
                "Density @ selected conditions: N.A.",
                f"Density/Phase Error: {error}",
            ]
        )
    elif phase == "Two-Phase" and density_result is not None:
        if isinstance(density_result, tuple) and len(density_result) == 2:
            density_liq, density_gas = density_result
            items.append(f"Phase @ selected conditions: {phase}")
            if density_liq is not None:
                items.append(
                    "  Density @ selected conditions (Liq): "
                    f"{density_liq:.3f} kg/m³"
                )
            else:
                items.append(
                    "  Density @ selected conditions (Liq): N.A."
                )
            if density_gas is not None:
                items.append(
                    "  Density @ selected conditions (Vap): "
                    f"{density_gas:.3f} kg/m³"
                )
            else:
                items.append(
                    "  Density @ selected conditions (Vap): N.A."
                )
        else:
            items.append(
                f"Phase @ selected conditions: {phase} "
                "(Result format unexpected)"
            )
    elif density_result is not None:
        items.append(
            "Phase @ selected conditions: "
            f"{phase or 'Could not determine.'}"
        )
        items.append(
            f"Density @ selected conditions: {density_result:.3f} kg/m³"
        )
    elif phase:
        items.append(f"Phase @ selected conditions: {phase}")
        items.append("Density @ selected conditions: N.A.")
    else:
        items.append("Phase @ selected conditions: Could not determine.")
        items.append("Density @ selected conditions: N.A.")

    items.append(
        _build_reference_density_line(
            "Density @ normal conditions",
            result_data.get("density_normal_kg_m3"),
            result_data.get("density_normal_phase"),
            result_data.get("density_normal_error"),
        )
    )
    items.append(
        _build_reference_density_line(
            "Density @ standard conditions",
            result_data.get("density_standard_kg_m3"),
            result_data.get("density_standard_phase"),
            result_data.get("density_standard_error"),
        )
    )

    bubble_point = result_data["bubble_point"]
    bp_error = result_data["bp_error"]
    pressure_atm = result_data["pressure_atm"]
    if bp_error:
        items.append(f"Bubble Point Error: {bp_error}")
    elif bubble_point is not None:
        items.append(
            f"Bubble Point @ {pressure_atm} atm: {bubble_point:.2f} °C"
        )
    else:
        items.append("Bubble Point: Calculation failed.")

    if lhv_data_loaded:
        _add_lhv_items(items, result_data)
    else:
        items.append("-" * 40)
        items.append("Mixture LHV = N/A (DB not loaded)")

    return tuple(items)


def _build_reference_density_line(
    label,
    density_value,
    phase_value,
    error_value,
) -> str:
    if density_value is not None:
        return f"{label}: {density_value:.3f} kg/m³"

    detail = ""
    if error_value:
        detail = f" ({error_value})"
    elif phase_value == "Two-Phase":
        detail = " (Two-Phase)"
    return f"{label}: N.A.{detail}"


def _add_lhv_items(items: list[str], result_data) -> None:
    mixture_lhv = result_data["mixture_lhv"]
    missing_lhv = result_data["missing_lhv"]
    if missing_lhv:
        items.append(f"LHV Warning: No data for: {', '.join(missing_lhv)}")
    items.append("-" * 40)
    lhv_display_values = build_lhv_display_values(
        mixture_lhv,
        result_data["mw"],
    )
    volumetric_lines = [
        ("Mixture LHV = ", "MJ/Nm³"),
        ("= ", "kcal/Nm³"),
        ("= ", "MMkcal/Nm³"),
        ("= ", "GJ/Nm³"),
        ("= ", "MMBtu/Nm³"),
    ]
    for prefix, unit in volumetric_lines:
        items.append(
            f"{prefix}"
            f"{format_lhv_display_value(lhv_display_values['volumetric'][unit])} "
            f"{unit}"
        )

    items.append("")
    items.append("Mass basis:")
    for unit in [
        "MJ/kg",
        "MJ/t",
        "GJ/kg",
        "GJ/t",
        "kcal/kg",
        "kcal/t",
        "MMkcal/kg",
        "MMkcal/t",
        "MMBtu/kg",
        "MMBtu/t",
    ]:
        items.append(
            f"= "
            f"{format_lhv_display_value(lhv_display_values['mass_basis'][unit])} "
            f"{unit}"
        )
