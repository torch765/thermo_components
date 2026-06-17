import pytest

from thermo_components.adapters.thermo import ThermoGateway
from thermo_components.domain.lhv import (
    build_lhv_display_values,
    calculate_mixture_lhv,
)


VOLUMETRIC_MJ_UNIT = "MJ/Nm\u00b3"


def test_lhv_display_values_preserve_base_and_derive_mass_basis():
    values = build_lhv_display_values(35.8, 16.04)

    assert values["volumetric"][VOLUMETRIC_MJ_UNIT] == pytest.approx(35.8)
    assert values["mass_basis"]["MJ/kg"] == pytest.approx(
        35.8 / (16.04 / 22.414)
    )
    assert values["mass_basis"]["MJ/t"] == pytest.approx(
        values["mass_basis"]["MJ/kg"] * 1000.0
    )


def test_lhv_mass_basis_is_unavailable_without_positive_molecular_weight():
    values = build_lhv_display_values(35.8, 0.0)
    assert all(value is None for value in values["mass_basis"].values())


def test_mixture_lhv_normalizes_composition_and_reports_missing_data():
    lhv, missing = calculate_mixture_lhv(
        {"methane": 2.0, "unknown": 1.0},
        {"methane": 35.8},
    )

    assert lhv == pytest.approx(35.8 * 2.0 / 3.0)
    assert missing == ["unknown"]


def test_thermo_gateway_keeps_lhv_method():
    calculator = ThermoGateway({"methane": 1.0})
    assert calculator.calculate_lhv({"methane": 35.8}) == (35.8, [])
