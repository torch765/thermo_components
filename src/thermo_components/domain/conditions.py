"""Operating and reference-condition primitives."""

from dataclasses import dataclass


ATM_TO_PA = 101325.0
R = 0.082057


@dataclass(frozen=True)
class ReferenceCondition:
    temperature_c: float
    pressure_atm: float

    @property
    def temperature_k(self) -> float:
        return celsius_to_kelvin(self.temperature_c)

    @property
    def pressure_pa(self) -> float:
        return atm_to_pa(self.pressure_atm)


NORMAL_CONDITION = ReferenceCondition(temperature_c=0.0, pressure_atm=1.0)
STANDARD_CONDITION = ReferenceCondition(
    temperature_c=(60.0 - 32.0) * 5.0 / 9.0,
    pressure_atm=1.0,
)

NORMAL_T_C = NORMAL_CONDITION.temperature_c
NORMAL_P_ATM = NORMAL_CONDITION.pressure_atm
STANDARD_T_F = 60.0
STANDARD_T_C = STANDARD_CONDITION.temperature_c
STANDARD_P_ATM = STANDARD_CONDITION.pressure_atm


def celsius_to_kelvin(temperature_c: float) -> float:
    return temperature_c + 273.15


def atm_to_pa(pressure_atm: float) -> float:
    return pressure_atm * ATM_TO_PA
