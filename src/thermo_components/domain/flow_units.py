"""Flow-unit definitions and reference bases."""

from dataclasses import asdict, dataclass


KG_PER_TONNE = 1000.0
KG_PER_LB = 0.45359237
LB_PER_KLB = 1000.0
HOURS_PER_DAY = 24.0
M3_PER_BBL = 0.158987294928
FT3_PER_M3 = 35.3146667
M3_PER_FT3 = 0.028316846592


@dataclass(frozen=True)
class FlowUnitDefinition:
    dimension: str
    basis: str
    amount_to_base: float
    time_to_day: float


FLOW_UNIT_ORDER = [
    "kg/h",
    "kg/d",
    "t/h",
    "t/d",
    "lb/h",
    "lb/d",
    "Klb/h",
    "Klb/d",
    "Nm3/h",
    "Nm3/d",
    "Sm3/h",
    "Sm3/d",
    "bbl/h",
    "bbl/d",
    "SCFH",
    "MSCFD",
    "MMSCFD",
]

FLOW_UNITS = {
    "kg/h": FlowUnitDefinition("mass", "none", 1.0, HOURS_PER_DAY),
    "kg/d": FlowUnitDefinition("mass", "none", 1.0, 1.0),
    "t/h": FlowUnitDefinition("mass", "none", KG_PER_TONNE, HOURS_PER_DAY),
    "t/d": FlowUnitDefinition("mass", "none", KG_PER_TONNE, 1.0),
    "lb/h": FlowUnitDefinition("mass", "none", KG_PER_LB, HOURS_PER_DAY),
    "lb/d": FlowUnitDefinition("mass", "none", KG_PER_LB, 1.0),
    "Klb/h": FlowUnitDefinition(
        "mass",
        "none",
        KG_PER_LB * LB_PER_KLB,
        HOURS_PER_DAY,
    ),
    "Klb/d": FlowUnitDefinition(
        "mass",
        "none",
        KG_PER_LB * LB_PER_KLB,
        1.0,
    ),
    "Nm3/h": FlowUnitDefinition("ref_volume", "normal", 1.0, HOURS_PER_DAY),
    "Nm3/d": FlowUnitDefinition("ref_volume", "normal", 1.0, 1.0),
    "Sm3/h": FlowUnitDefinition("ref_volume", "standard", 1.0, HOURS_PER_DAY),
    "Sm3/d": FlowUnitDefinition("ref_volume", "standard", 1.0, 1.0),
    "bbl/h": FlowUnitDefinition(
        "ref_volume",
        "standard",
        M3_PER_BBL,
        HOURS_PER_DAY,
    ),
    "bbl/d": FlowUnitDefinition("ref_volume", "standard", M3_PER_BBL, 1.0),
    "SCFH": FlowUnitDefinition(
        "ref_volume",
        "standard",
        M3_PER_FT3,
        HOURS_PER_DAY,
    ),
    "MSCFD": FlowUnitDefinition(
        "ref_volume",
        "standard",
        1000.0 * M3_PER_FT3,
        1.0,
    ),
    "MMSCFD": FlowUnitDefinition(
        "ref_volume",
        "standard",
        1_000_000.0 * M3_PER_FT3,
        1.0,
    ),
}

# Compatibility mapping retained while density.py remains the application entry point.
FLOW_UNIT_DEFINITIONS = {
    unit: asdict(definition)
    for unit, definition in FLOW_UNITS.items()
}
