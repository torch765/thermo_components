import pytest

from thermo_components.adapters.ui.qt_main_window import MainWindow
from thermo_components.domain.composition import (
    derive_inactive_percentages,
    is_effectively_pure_water,
    normalize_component_identity,
    normalize_percentages,
    percentages_to_mole_fractions,
    water_fraction_active_basis,
)
from thermo_components.domain.thermo_routes import (
    IAPWS95_MODEL_DISPLAY,
    PRMIX_DEFAULT_ROUTE,
    PURE_WATER_ROUTE,
    select_thermo_route,
)


def test_component_identity_normalizes_water_alias():
    assert normalize_component_identity(" H2O ") == "water"
    assert normalize_component_identity("Methane") == "methane"


def test_water_fraction_uses_the_active_input_basis():
    fraction = water_fraction_active_basis(
        ["water", "methane"],
        [25.0, 75.0],
        [80.0, 20.0],
        "Wt %",
    )
    assert fraction == pytest.approx(0.8)


def test_effectively_pure_water_uses_the_existing_threshold():
    assert is_effectively_pure_water(
        ["water", "methane"],
        [99.9, 0.1],
        [0.0, 0.0],
        "Mol %",
    )
    assert not is_effectively_pure_water(
        ["water", "methane"],
        [99.89, 0.11],
        [0.0, 0.0],
        "Mol %",
    )


def test_route_selection_preserves_current_model_policy():
    pure_water = select_thermo_route(
        ["water"],
        [100.0],
        [100.0],
        "Mol %",
        "PRMIX",
    )
    mixture = select_thermo_route(
        ["water", "methane"],
        [50.0, 50.0],
        [0.0, 0.0],
        "Mol %",
        "PRMIX",
    )

    assert pure_water == {
        "route_id": PURE_WATER_ROUTE,
        "model_display": IAPWS95_MODEL_DISPLAY,
    }
    assert mixture == {
        "route_id": PRMIX_DEFAULT_ROUTE,
        "model_display": "PRMIX",
    }


def test_basis_conversion_derives_weight_and_mole_percentages():
    weight_percentages = derive_inactive_percentages(
        ["methane", "ethane"],
        [50.0, 50.0],
        "Mol %",
    )
    assert weight_percentages == pytest.approx(
        [16.04 / (16.04 + 30.07) * 100.0, 30.07 / (16.04 + 30.07) * 100.0]
    )

    mole_percentages = derive_inactive_percentages(
        ["methane", "ethane"],
        weight_percentages,
        "Wt %",
    )
    assert mole_percentages == pytest.approx([50.0, 50.0])


def test_normalization_and_weight_to_mole_fraction_conversion():
    assert normalize_percentages([40.0, 40.0]) == [50.0, 50.0]
    assert percentages_to_mole_fractions(
        ["methane", "ethane"],
        [50.0, 50.0],
        "Mol %",
    ) == {
        "methane": 0.5,
        "ethane": 0.5,
    }

    mole_fractions = percentages_to_mole_fractions(
        ["methane", "ethane"],
        [50.0, 50.0],
        "Wt %",
    )
    assert sum(mole_fractions.values()) == pytest.approx(1.0)
    assert mole_fractions["methane"] > mole_fractions["ethane"]


def test_normalize_composition_scales_active_column_to_100(qt_app):
    window = MainWindow(lhv_data={})

    for component in ("methane", "ethane"):
        index = window.ui.comboBox_select_components.findText(component)
        window.ui.comboBox_select_components.setCurrentIndex(index)

    window.ui.tableWidget.item(0, 1).setText("40")
    window.ui.tableWidget.item(1, 1).setText("40")
    window.normalize_composition()

    assert float(window.ui.tableWidget.item(0, 1).text()) == pytest.approx(50.0)
    assert float(window.ui.tableWidget.item(1, 1).text()) == pytest.approx(50.0)
    assert float(window.ui.tableWidget.item(2, 1).text()) == pytest.approx(100.0)

    window.close()
