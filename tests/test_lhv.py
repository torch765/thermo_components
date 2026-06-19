import pytest

from lhv_data import LHV_DATA_RAW
from thermo_components.adapters.thermo import ThermoGateway
from thermo_components.domain.composition import MOLECULAR_WEIGHTS
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


def test_lhv_lookup_is_case_insensitive_for_canonical_component_names():
    lhv, missing = calculate_mixture_lhv(
        {"MTBE": 1.0},
        {"mtbe": 139.9},
    )

    assert lhv == pytest.approx(139.9)
    assert missing == []


def test_lhv_seed_covers_every_selectable_component():
    selectable_components = {
        component_name
        for component_name, molecular_weight in MOLECULAR_WEIGHTS.items()
        if component_name.strip() and molecular_weight > 0
    }

    assert set(LHV_DATA_RAW) == selectable_components
    assert "n-butylene" not in LHV_DATA_RAW
    assert "2-butene" not in LHV_DATA_RAW


def test_new_component_lhv_values_match_the_documented_seed():
    expected_values = {
        "ammonia": 14.1,
        "carbonyl sulfide": 24.6,
        "sulfur dioxide": 0.0,
        "carbon disulfide": 49.3,
        "methyl mercaptan": 51.4,
        "methanol": 30.1,
        "MTBE": 139.9,
        "dimethyl ether": 59.3,
        "decane": 283.1,
        "dodecane": 338.0,
        "cyclohexane": 164.6,
        "methylcyclohexane": 191.5,
        "cis-2-butene": 113.0,
        "trans-2-butene": 112.9,
    }

    assert {
        component_name: LHV_DATA_RAW[component_name]
        for component_name in expected_values
    } == expected_values
