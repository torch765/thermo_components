"""Mass and reference-volume flow conversion rules."""

from .flow_units import FLOW_UNITS, FlowUnitDefinition


def parse_flow_input(text: str) -> float | None:
    """Parse a flow input, tolerating spaces and thousands separators."""
    cleaned = str(text).strip().replace(",", "").replace(" ", "")
    if not cleaned:
        return None
    return float(cleaned)


def format_flow_value(value: float) -> str:
    """Format flow results without scientific notation for typical values."""
    formatted = f"{value:,.6f}"
    if "." in formatted:
        formatted = formatted.rstrip("0").rstrip(".")
    return formatted


def _coerce_positive_density(value) -> float | None:
    try:
        density = float(value)
    except (TypeError, ValueError):
        return None
    if density <= 0:
        return None
    return density


def _required_density_bases(
    from_unit: FlowUnitDefinition,
    to_unit: FlowUnitDefinition,
) -> set[str]:
    if from_unit.dimension == "mass" and to_unit.dimension == "mass":
        return set()
    if from_unit.dimension == "mass":
        return {to_unit.basis}
    if to_unit.dimension == "mass":
        return {from_unit.basis}
    if from_unit.basis == to_unit.basis:
        return set()
    return {from_unit.basis, to_unit.basis}


def _missing_density_message(missing_bases: set[str]) -> str:
    if missing_bases == {"normal"}:
        return "Normal density required"
    if missing_bases == {"standard"}:
        return "Standard density required"
    return "Normal and standard densities required"


def convert_flow(
    value,
    from_unit,
    to_unit,
    density_normal_kg_m3=None,
    density_standard_kg_m3=None,
):
    """Convert between mass and reference-volume flow units."""
    from_definition = FLOW_UNITS.get(from_unit)
    to_definition = FLOW_UNITS.get(to_unit)
    if from_definition is None or to_definition is None:
        raise ValueError("Select source and destination units")

    value = float(value)
    required_bases = _required_density_bases(from_definition, to_definition)
    density_map = {
        "normal": _coerce_positive_density(density_normal_kg_m3),
        "standard": _coerce_positive_density(density_standard_kg_m3),
    }
    missing_bases = {
        basis
        for basis in required_bases
        if density_map.get(basis) is None
    }
    if missing_bases:
        raise ValueError(_missing_density_message(missing_bases))

    from_base_per_day = (
        value
        * from_definition.amount_to_base
        * from_definition.time_to_day
    )

    if from_definition.dimension == "mass":
        mass_kg_day = from_base_per_day
        if to_definition.dimension == "mass":
            return mass_kg_day / (
                to_definition.amount_to_base * to_definition.time_to_day
            )
        target_volume_m3_day = mass_kg_day / density_map[to_definition.basis]
        return target_volume_m3_day / (
            to_definition.amount_to_base * to_definition.time_to_day
        )

    source_volume_m3_day = from_base_per_day
    if (
        to_definition.dimension == "ref_volume"
        and from_definition.basis == to_definition.basis
    ):
        return source_volume_m3_day / (
            to_definition.amount_to_base * to_definition.time_to_day
        )

    mass_kg_day = source_volume_m3_day * density_map[from_definition.basis]
    if to_definition.dimension == "mass":
        return mass_kg_day / (
            to_definition.amount_to_base * to_definition.time_to_day
        )

    target_volume_m3_day = mass_kg_day / density_map[to_definition.basis]
    return target_volume_m3_day / (
        to_definition.amount_to_base * to_definition.time_to_day
    )
