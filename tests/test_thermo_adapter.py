import pytest

from thermo_components.adapters.thermo import ThermoGateway
from thermo_components.application.ports import ThermoPropertyGateway
from thermo_components.domain.thermo_routes import (
    PRMIX_DEFAULT_ROUTE,
    PURE_WATER_ROUTE,
)


def test_thermo_gateway_implements_port():
    gateway = ThermoGateway()

    assert isinstance(gateway, ThermoPropertyGateway)


def test_thermo_gateway_calculates_methane_with_prmix():
    gateway = ThermoGateway({"methane": 1.0})

    density, phase, density_error = gateway.calculate_density_for_route(
        298.15,
        101325.0,
        PRMIX_DEFAULT_ROUTE,
    )
    bubble_point_c, bubble_point_error = (
        gateway.calculate_bubble_point_for_route(
            101325.0,
            PRMIX_DEFAULT_ROUTE,
        )
    )

    assert density == pytest.approx(0.6571965, rel=1e-5)
    assert phase == "Vapor"
    assert density_error is None
    assert bubble_point_c == pytest.approx(-161.5699, abs=1e-3)
    assert bubble_point_error is None


def test_thermo_gateway_routes_pure_water_to_iapws95():
    gateway = ThermoGateway({"water": 1.0})

    density, phase, density_error = gateway.calculate_density_for_route(
        298.15,
        101325.0,
        PURE_WATER_ROUTE,
    )
    saturation_c, saturation_error = (
        gateway.calculate_bubble_point_for_route(
            101325.0,
            PURE_WATER_ROUTE,
        )
    )

    assert density == pytest.approx(997.0476, rel=1e-5)
    assert phase == "Liquid"
    assert density_error is None
    assert saturation_c == pytest.approx(99.9743, abs=1e-3)
    assert saturation_error is None
