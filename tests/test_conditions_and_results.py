import pytest

from thermo_components.domain.conditions import (
    NORMAL_CONDITION,
    STANDARD_CONDITION,
    atm_to_pa,
    celsius_to_kelvin,
)
from thermo_components.domain.results import (
    build_density_note,
    extract_scalar_density_value,
)


def test_reference_conditions_expose_kelvin_and_pascal_values():
    assert NORMAL_CONDITION.temperature_k == pytest.approx(273.15)
    assert NORMAL_CONDITION.pressure_pa == pytest.approx(101325.0)
    assert STANDARD_CONDITION.temperature_c == pytest.approx(15.5555555556)
    assert celsius_to_kelvin(25.0) == pytest.approx(298.15)
    assert atm_to_pa(5.0) == pytest.approx(506625.0)


def test_density_result_helpers_preserve_scalar_and_two_phase_semantics():
    assert extract_scalar_density_value(1.25, "Vapor", None) == pytest.approx(1.25)
    assert extract_scalar_density_value((500.0, 2.0), "Two-Phase", None) is None
    assert extract_scalar_density_value(1.25, "Vapor", "failed") is None
    assert build_density_note("Liquid", None) == "Phase: Liquid"
    assert build_density_note("Two-Phase", None) == "Two-Phase result"
    assert build_density_note(None, "failed") == "failed"
