"""Thermodynamic route selection policy."""

from enum import Enum

from .composition import is_effectively_pure_water


class ThermoRoute(str, Enum):
    PURE_WATER_IAPWS95 = "iapws95_pure_water"
    PRMIX_DEFAULT = "prmix_default"


PURE_WATER_ROUTE = ThermoRoute.PURE_WATER_IAPWS95.value
PRMIX_DEFAULT_ROUTE = ThermoRoute.PRMIX_DEFAULT.value
IAPWS95_MODEL_DISPLAY = "IAPWS-95"
IAPWS95_TWO_PHASE_REL_TOL = 1e-4


def select_thermo_route(
    component_names,
    mole_percents,
    weight_percents,
    basis: str,
    default_eos: str,
) -> dict:
    """Select the calculation route for the current composition."""
    if is_effectively_pure_water(
        component_names,
        mole_percents,
        weight_percents,
        basis,
    ):
        return {
            "route_id": PURE_WATER_ROUTE,
            "model_display": IAPWS95_MODEL_DISPLAY,
        }
    return {
        "route_id": PRMIX_DEFAULT_ROUTE,
        "model_display": str(default_eos or "PRMIX").strip() or "PRMIX",
    }
