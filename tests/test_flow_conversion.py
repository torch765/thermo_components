import pytest

from thermo_components.domain.flow_conversion import (
    convert_flow,
    format_flow_value,
    parse_flow_input,
)


def test_mass_flow_conversion_does_not_require_density():
    assert convert_flow(1.0, "kg/h", "t/d") == pytest.approx(0.024)


def test_same_reference_basis_conversion_does_not_require_density():
    assert convert_flow(1.0, "Nm3/h", "Nm3/d") == pytest.approx(24.0)
    assert convert_flow(1_000.0, "SCFH", "MSCFD") == pytest.approx(24.0)


def test_mass_to_reference_volume_uses_target_basis_density():
    assert convert_flow(
        10.0,
        "kg/h",
        "Nm3/h",
        density_normal_kg_m3=2.0,
    ) == pytest.approx(5.0)


def test_cross_basis_volume_conversion_bridges_through_mass():
    assert convert_flow(
        10.0,
        "Nm3/h",
        "Sm3/h",
        density_normal_kg_m3=1.2,
        density_standard_kg_m3=1.0,
    ) == pytest.approx(12.0)


@pytest.mark.parametrize(
    ("from_unit", "to_unit", "message"),
    [
        ("kg/h", "Nm3/h", "Normal density required"),
        ("kg/h", "Sm3/h", "Standard density required"),
        ("Nm3/h", "Sm3/h", "Normal and standard densities required"),
    ],
)
def test_conversion_reports_required_density_basis(from_unit, to_unit, message):
    with pytest.raises(ValueError, match=f"^{message}$"):
        convert_flow(1.0, from_unit, to_unit)


def test_flow_input_and_output_formatting_are_stable():
    assert parse_flow_input(" 1,234.50 ") == pytest.approx(1234.5)
    assert parse_flow_input("   ") is None
    assert format_flow_value(1234.5) == "1,234.5"
