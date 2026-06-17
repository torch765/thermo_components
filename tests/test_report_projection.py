from thermo_components.adapters.ui.qt_main_window import MainWindow
from thermo_components.domain.thermo_routes import PRMIX_DEFAULT_ROUTE
from thermo_components.domain.warnings import PRMIX_WATER_WARNING


def result_by_property(rows, property_name):
    return [row for row in rows if row["Property"] == property_name]


def test_report_projection_builds_results_without_scraping_results_widget(qt_app):
    window = MainWindow(lhv_data={"methane": 35.8})
    result_data = {
        "mw": 16.04,
        "phase": "Vapor",
        "density_result": 0.72,
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
    }

    rows = window.build_results_rows(result_data)

    assert result_by_property(rows, "Average molecular weight")[0]["Value"] == 16.04
    assert result_by_property(rows, "Phase @ selected conditions")[0]["Value"] == "Vapor"
    assert result_by_property(rows, "Density @ normal conditions")[0]["Value"] == 0.716
    assert result_by_property(rows, "Density @ standard conditions")[0]["Value"] == 0.68
    assert len(result_by_property(rows, "Mixture LHV")) == 5
    assert len(result_by_property(rows, "Mixture LHV (Mass basis)")) == 10

    window.close()


def test_report_warning_projection_combines_thermo_and_lhv_warnings(qt_app):
    window = MainWindow(lhv_data={"methane": 35.8})

    rows = window.build_report_warning_rows(
        {
            "warnings": [PRMIX_WATER_WARNING],
            "missing_lhv": ["unknown"],
        }
    )

    assert rows == [
        {
            "Warning Type": "Thermo",
            "Details": PRMIX_WATER_WARNING,
        },
        {
            "Warning Type": "LHV",
            "Details": "LHV Warning: No data for: unknown",
        },
    ]

    window.close()
