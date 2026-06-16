from thermo_components.adapters.ui import build_result_list_items
from thermo_components.domain.thermo_routes import PRMIX_DEFAULT_ROUTE


def methane_result_data():
    return {
        "mw": 16.04,
        "basis": "Mol %",
        "model_display": "PRMIX",
        "eos": "PRMIX",
        "phase": "Vapor",
        "density_result": 0.6571965,
        "density_error": None,
        "density_normal_kg_m3": 0.716,
        "density_normal_phase": "Vapor",
        "density_normal_error": None,
        "density_standard_kg_m3": 0.68,
        "density_standard_phase": "Vapor",
        "density_standard_error": None,
        "pressure_atm": 1.0,
        "bubble_point": -161.5,
        "bp_error": None,
        "thermo_route": PRMIX_DEFAULT_ROUTE,
        "mixture_lhv": 35.8,
        "missing_lhv": [],
        "warnings": [],
    }


def test_result_list_presenter_builds_existing_vapor_display_lines():
    items = build_result_list_items(
        methane_result_data(),
        lhv_data_loaded=True,
    )

    assert items[:8] == (
        "Calculating for Mol % basis...",
        "Model / EOS: PRMIX",
        "Avg. Molecular Wt: 16.04 g/mol",
        "Phase @ selected conditions: Vapor",
        "Density @ selected conditions: 0.657 kg/m³",
        "Density @ normal conditions: 0.716 kg/m³",
        "Density @ standard conditions: 0.680 kg/m³",
        "Bubble Point @ 1.0 atm: -161.50 °C",
    )
    assert "Mixture LHV = 35.80 MJ/Nm³" in items
    assert "Mass basis:" in items


def test_result_list_presenter_preserves_no_lhv_database_message():
    items = build_result_list_items(
        methane_result_data(),
        lhv_data_loaded=False,
    )

    assert items[-2:] == (
        "----------------------------------------",
        "Mixture LHV = N/A (DB not loaded)",
    )


def test_result_list_presenter_formats_two_phase_density_values():
    data = methane_result_data()
    data.update(
        {
            "phase": "Two-Phase",
            "density_result": (500.0, 2.5),
            "density_normal_kg_m3": None,
            "density_normal_phase": "Two-Phase",
            "density_standard_kg_m3": None,
            "density_standard_error": "flash failed",
        }
    )

    items = build_result_list_items(data, lhv_data_loaded=False)

    assert "Phase @ selected conditions: Two-Phase" in items
    assert "  Density @ selected conditions (Liq): 500.000 kg/m³" in items
    assert "  Density @ selected conditions (Vap): 2.500 kg/m³" in items
    assert "Density @ normal conditions: N.A. (Two-Phase)" in items
    assert "Density @ standard conditions: N.A. (flash failed)" in items


def test_main_window_renders_presented_result_items(qt_app):
    from density import MainWindow

    window = MainWindow(lhv_data={"methane": 35.8})

    window.on_calculation_result(methane_result_data())

    assert window.ui.results_list.item(0).text() == (
        "Calculating for Mol % basis..."
    )
    assert window.ui.results_list.item(4).text() == (
        "Density @ selected conditions: 0.657 kg/m³"
    )
    assert window.density_actual_kg_m3 == 0.6571965
    assert window.density_normal_kg_m3 == 0.716
    assert window.density_standard_kg_m3 == 0.68

    window.close()
