"""Composition identities, basis conversion, and normalization rules."""

from collections.abc import Mapping, Sequence
from enum import Enum


class CompositionBasis(str, Enum):
    MOLE_PERCENT = "Mol %"
    WEIGHT_PERCENT = "Wt %"


MOLECULAR_WEIGHTS = {
    "        ": 0.0,
    "hydrogen": 2.016,
    "carbon dioxide": 44.01,
    "carbon monoxide": 28.01,
    "hydrogen sulfide": 34.08,
    "ammonia": 17.03,
    "carbonyl sulfide": 60.08,
    "sulfur dioxide": 64.06,
    "carbon disulfide": 76.14,
    "methyl mercaptan": 48.11,
    "methane": 16.04,
    "ethane": 30.07,
    "ethylene": 28.05,
    "propane": 44.10,
    "propylene": 42.08,
    "isobutane": 58.12,
    "n-butane": 58.12,
    "isobutylene": 56.11,
    "1-butene": 56.11,
    "cis-2-butene": 56.11,
    "trans-2-butene": 56.11,
    "butadiene": 54.09,
    "isopentane": 72.15,
    "n-pentane": 72.15,
    "hexane": 86.18,
    "cyclohexane": 84.16,
    "heptane": 100.20,
    "methylcyclohexane": 98.19,
    "benzene": 78.11,
    "toluene": 92.14,
    "o-xylene": 106.17,
    "m-xylene": 106.17,
    "p-xylene": 106.17,
    "ethylbenzene": 106.17,
    "cumene": 120.19,
    "octane": 114.23,
    "nonane": 128.26,
    "decane": 142.28,
    "dodecane": 170.33,
    "methanol": 32.04,
    "MTBE": 88.15,
    "dimethyl ether": 46.07,
    "oxygen": 32.00,
    "nitrogen": 28.01,
    "argon": 39.95,
    "water": 18.015,
}

WATER_COMPONENT_ALIASES = {"water", "h2o"}
PURE_WATER_WARNING_FRACTION = 0.999


def normalize_component_identity(name: str) -> str:
    """Normalize a component identifier for domain comparisons."""
    cleaned = str(name).strip()
    lowercase = cleaned.lower()
    compact = "".join(
        character for character in lowercase if character.isalnum()
    )
    if compact == "h2o":
        return "water"
    return next(
        (
            component_name
            for component_name in MOLECULAR_WEIGHTS
            if component_name.lower() == lowercase
        ),
        lowercase,
    )


def is_water_component(name: str) -> bool:
    """Return whether a component identifier represents water."""
    return normalize_component_identity(name) in WATER_COMPONENT_ALIASES


def active_basis_amount_rows(
    component_names: Sequence[str],
    mole_percents: Sequence[float],
    weight_percents: Sequence[float],
    basis: str,
) -> tuple[list[tuple[str, float]], float]:
    """Return positive amounts on the selected input basis."""
    active_values = (
        mole_percents
        if basis == CompositionBasis.MOLE_PERCENT.value
        else weight_percents
    )
    rows: list[tuple[str, float]] = []
    total_active = 0.0

    for component_name, raw_value in zip(component_names, active_values):
        try:
            value = float(raw_value)
        except (TypeError, ValueError):
            continue

        if value <= 0:
            continue

        rows.append((normalize_component_identity(component_name), value))
        total_active += value

    return rows, total_active


def has_water_component(
    component_names: Sequence[str],
    mole_percents: Sequence[float],
    weight_percents: Sequence[float],
    basis: str,
) -> bool:
    """Return whether water has a positive amount on the selected basis."""
    rows, _ = active_basis_amount_rows(
        component_names,
        mole_percents,
        weight_percents,
        basis,
    )
    return any(is_water_component(name) for name, _ in rows)


def water_fraction_active_basis(
    component_names: Sequence[str],
    mole_percents: Sequence[float],
    weight_percents: Sequence[float],
    basis: str,
) -> float:
    """Return the water fraction on the selected input basis."""
    rows, total_active = active_basis_amount_rows(
        component_names,
        mole_percents,
        weight_percents,
        basis,
    )
    if total_active <= 0:
        return 0.0
    water_active = sum(value for name, value in rows if is_water_component(name))
    return water_active / total_active


def is_effectively_pure_water(
    component_names: Sequence[str],
    mole_percents: Sequence[float],
    weight_percents: Sequence[float],
    basis: str,
) -> bool:
    """Return whether the selected-basis composition meets the water threshold."""
    rows, total_active = active_basis_amount_rows(
        component_names,
        mole_percents,
        weight_percents,
        basis,
    )
    if total_active <= 0:
        return False

    water_active = sum(value for name, value in rows if is_water_component(name))
    other_active = total_active - water_active
    if water_active <= 0:
        return False

    water_fraction = water_active / total_active
    other_fraction = other_active / total_active
    return (
        water_fraction >= PURE_WATER_WARNING_FRACTION
        and other_fraction <= (1.0 - PURE_WATER_WARNING_FRACTION)
    )


def derive_inactive_percentages(
    component_names: Sequence[str],
    active_percentages: Sequence[float],
    basis: str,
    molecular_weights: Mapping[str, float] = MOLECULAR_WEIGHTS,
) -> list[float | None]:
    """Convert the active composition basis into the inactive display basis."""
    weights = [molecular_weights.get(name, 0.0) or 0.0 for name in component_names]

    if basis == CompositionBasis.MOLE_PERCENT.value:
        masses = [
            percentage * molecular_weight
            if molecular_weight > 0 and percentage != 0
            else 0.0
            for percentage, molecular_weight in zip(active_percentages, weights)
        ]
        total = sum(masses)
        if total <= 0:
            return [None] * len(component_names)
        return [100.0 * mass / total for mass in masses]

    moles = [
        percentage / molecular_weight
        if molecular_weight > 0 and percentage != 0
        else 0.0
        for percentage, molecular_weight in zip(active_percentages, weights)
    ]
    total = sum(moles)
    if total <= 0:
        return [None] * len(component_names)
    return [100.0 * mole_amount / total for mole_amount in moles]


def normalize_percentages(
    percentages: Sequence[float],
    precision: int = 4,
) -> list[float]:
    """Normalize percentages to exactly 100 at the requested precision."""
    values = [float(value) for value in percentages]
    total = sum(values)
    if total == 0:
        raise ValueError("Cannot normalize: total is zero.")

    normalized: list[float] = []
    for index, value in enumerate(values):
        if index < len(values) - 1:
            normalized.append(round((value / total) * 100.0, precision))
        else:
            normalized.append(round(100.0 - sum(normalized), precision))
    return normalized


def percentages_to_mole_fractions(
    component_names: Sequence[str],
    percentages: Sequence[float],
    basis: str,
    molecular_weights: Mapping[str, float] = MOLECULAR_WEIGHTS,
) -> dict[str, float]:
    """Convert percentages on either supported basis into mole fractions."""
    if basis == CompositionBasis.MOLE_PERCENT.value:
        return {
            name: percentage / 100.0
            for name, percentage in zip(component_names, percentages)
        }

    mole_amounts: dict[str, float] = {}
    total_moles = 0.0
    for name, percentage in zip(component_names, percentages):
        molecular_weight = molecular_weights.get(name)
        if molecular_weight is None or molecular_weight <= 0:
            raise ValueError(f"Missing or invalid MW for {name}.")
        mole_amounts[name] = percentage / molecular_weight
        total_moles += mole_amounts[name]

    if total_moles == 0:
        raise ValueError("Total moles is zero after Wt% conversion.")

    return {
        name: mole_amounts[name] / total_moles
        for name in component_names
    }


def calculate_mixture_molecular_weight(
    components: Mapping[str, float],
    molecular_weights: Mapping[str, float] = MOLECULAR_WEIGHTS,
) -> float:
    """Calculate normalized mixture molecular weight in g/mol."""
    if not components:
        return 0.0

    fraction_sum = sum(components.values())
    if fraction_sum == 0:
        return 0.0

    return sum(
        molecular_weights.get(component_name, 0.0) * fraction / fraction_sum
        for component_name, fraction in components.items()
    )
