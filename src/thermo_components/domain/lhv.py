"""Lower-heating-value mixture and display rules."""

from collections.abc import Mapping


KCAL_PER_MJ = 238.845896627
MJ_PER_MMBTU = 1055.05585262
NORMAL_MOLAR_VOLUME_NM3_PER_KMOL = 22.414
CUBIC_METRE = "\u00b3"


def calculate_mixture_lhv(
    components: Mapping[str, float],
    lhv_data: Mapping[str, float],
) -> tuple[float, list[str]]:
    """Calculate mixture LHV in MJ/Nm3 from mole fractions."""
    if not components:
        return 0.0, []

    total_fraction = sum(components.values())
    if total_fraction == 0:
        return 0.0, []

    normalized_lhv_data = {
        str(component_name).strip().casefold(): float(value)
        for component_name, value in lhv_data.items()
    }
    mixture_lhv = 0.0
    missing_components: list[str] = []
    for component_name, fraction in components.items():
        mole_fraction = fraction / total_fraction
        component_lhv = normalized_lhv_data.get(
            str(component_name).strip().casefold()
        )
        if component_lhv is not None:
            mixture_lhv += mole_fraction * component_lhv
        elif component_name.strip():
            missing_components.append(component_name)

    return mixture_lhv, missing_components


def build_lhv_display_values(lhv_mj_nm3: float, mw_g_mol: float) -> dict:
    """Build display-ready LHV values from the base MJ/Nm3 result."""
    try:
        lhv_mj_nm3 = float(lhv_mj_nm3)
    except (TypeError, ValueError):
        lhv_mj_nm3 = 0.0

    try:
        mw_g_mol = float(mw_g_mol)
    except (TypeError, ValueError):
        mw_g_mol = 0.0

    kcal_per_nm3 = lhv_mj_nm3 * KCAL_PER_MJ
    volumetric_units = {
        f"MJ/Nm{CUBIC_METRE}": lhv_mj_nm3,
        f"kcal/Nm{CUBIC_METRE}": kcal_per_nm3,
        f"MMkcal/Nm{CUBIC_METRE}": kcal_per_nm3 / 1_000_000.0,
        f"GJ/Nm{CUBIC_METRE}": lhv_mj_nm3 / 1000.0,
        f"MMBtu/Nm{CUBIC_METRE}": lhv_mj_nm3 / MJ_PER_MMBTU,
    }
    mass_basis = {
        "MJ/kg": None,
        "MJ/t": None,
        "GJ/kg": None,
        "GJ/t": None,
        "kcal/kg": None,
        "kcal/t": None,
        "MMkcal/kg": None,
        "MMkcal/t": None,
        "MMBtu/kg": None,
        "MMBtu/t": None,
    }

    if mw_g_mol > 0:
        kg_per_nm3 = mw_g_mol / NORMAL_MOLAR_VOLUME_NM3_PER_KMOL
        if kg_per_nm3 > 0:
            mj_per_kg = lhv_mj_nm3 / kg_per_nm3
            mj_per_t = mj_per_kg * 1000.0
            kcal_per_kg = mj_per_kg * KCAL_PER_MJ
            kcal_per_t = kcal_per_kg * 1000.0
            mass_basis = {
                "MJ/kg": mj_per_kg,
                "MJ/t": mj_per_t,
                "GJ/kg": mj_per_kg / 1000.0,
                "GJ/t": mj_per_t / 1000.0,
                "kcal/kg": kcal_per_kg,
                "kcal/t": kcal_per_t,
                "MMkcal/kg": kcal_per_kg / 1_000_000.0,
                "MMkcal/t": kcal_per_t / 1_000_000.0,
                "MMBtu/kg": mj_per_kg / MJ_PER_MMBTU,
                "MMBtu/t": mj_per_t / MJ_PER_MMBTU,
            }

    return {
        "volumetric": volumetric_units,
        "mass_basis": mass_basis,
    }


def format_lhv_display_value(value: float | None) -> str:
    """Format an LHV value with the existing display precision."""
    if value is None:
        return "N/A"
    decimals = 2 if abs(value) >= 1 else 4
    return f"{value:.{decimals}f}"
