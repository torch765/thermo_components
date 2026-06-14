"""Thermodynamic warning policies."""

from .composition import has_water_component
from .thermo_routes import PRMIX_DEFAULT_ROUTE


PRMIX_WATER_WARNING = (
    "Warning: Water is present. PRMIX results may be unreliable for aqueous "
    "or polar behavior. Use with caution."
)
PRMIX_TWO_PHASE_WATER_WARNING = (
    "Warning: Water-containing two-phase behavior may be unreliable with "
    "PRMIX. Validate density and phase behavior with a water-capable method."
)


def phase_indicates_two_phase(phase) -> bool:
    """Return whether phase text represents two-phase behavior."""
    return "two-phase" in str(phase or "").strip().lower()


def build_thermo_warning_messages(
    component_names,
    mole_percents,
    weight_percents,
    basis: str,
    thermo_route: str,
    result_data: dict,
) -> list[str]:
    """Build persistent warnings for a completed property calculation."""
    if thermo_route != PRMIX_DEFAULT_ROUTE:
        return []

    if not has_water_component(
        component_names,
        mole_percents,
        weight_percents,
        basis,
    ):
        return []

    messages = [PRMIX_WATER_WARNING]
    has_two_phase_warning_condition = any(
        phase_indicates_two_phase(result_data.get(key))
        for key in (
            "phase",
            "density_normal_phase",
            "density_standard_phase",
        )
    )
    if has_two_phase_warning_condition:
        messages.append(PRMIX_TWO_PHASE_WATER_WARNING)

    return messages
