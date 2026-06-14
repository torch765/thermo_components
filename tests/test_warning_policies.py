from thermo_components.domain.thermo_routes import (
    PRMIX_DEFAULT_ROUTE,
    PURE_WATER_ROUTE,
)
from thermo_components.domain.warnings import (
    PRMIX_TWO_PHASE_WATER_WARNING,
    PRMIX_WATER_WARNING,
    build_thermo_warning_messages,
)


def build_warnings(component_names, mol_percents, route, **result_data):
    return build_thermo_warning_messages(
        component_names,
        mol_percents,
        [0.0] * len(component_names),
        "Mol %",
        route,
        result_data,
    )


def test_prmix_water_mixture_produces_persistent_warning():
    warnings = build_warnings(
        ["water", "methane"],
        [10.0, 90.0],
        PRMIX_DEFAULT_ROUTE,
        phase="Vapor",
    )
    assert warnings == [PRMIX_WATER_WARNING]


def test_two_phase_reference_result_adds_specific_water_warning():
    warnings = build_warnings(
        ["water", "methane"],
        [10.0, 90.0],
        PRMIX_DEFAULT_ROUTE,
        phase="Vapor",
        density_normal_phase="Two-Phase",
    )
    assert warnings == [
        PRMIX_WATER_WARNING,
        PRMIX_TWO_PHASE_WATER_WARNING,
    ]


def test_non_water_and_iapws_routes_do_not_produce_prmix_warning():
    assert build_warnings(
        ["methane"],
        [100.0],
        PRMIX_DEFAULT_ROUTE,
        phase="Vapor",
    ) == []
    assert build_warnings(
        ["water"],
        [100.0],
        PURE_WATER_ROUTE,
        phase="Liquid",
    ) == []
