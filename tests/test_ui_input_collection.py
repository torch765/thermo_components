import pytest

from thermo_components.adapters.ui import (
    collect_property_calculation_request,
)


class FakeItem:
    def __init__(self, text):
        self._text = text

    def text(self):
        return self._text


class FakeTable:
    def __init__(self, rows):
        self._rows = rows

    def rowCount(self):
        return len(self._rows) + 1

    def item(self, row, column):
        try:
            value = self._rows[row][column]
        except IndexError:
            return None
        return None if value is None else FakeItem(value)


class FakeRadio:
    def __init__(self, checked):
        self._checked = checked

    def isChecked(self):
        return self._checked


class FakeCombo:
    def __init__(self, text):
        self._text = text

    def currentText(self):
        return self._text


class FakeUi:
    def __init__(
        self,
        rows,
        *,
        mol_basis=True,
        temperature="25 °C",
        pressure="1 atm",
    ):
        self.tableWidget = FakeTable(rows)
        self.radioButton_mol_percent = FakeRadio(mol_basis)
        self.radioButton_wt_percent = FakeRadio(not mol_basis)
        self.comboBox_select_temperature = FakeCombo(temperature)
        self.comboBox_select_pressure = FakeCombo(pressure)


def test_input_collector_builds_property_request_from_widgets():
    result = collect_property_calculation_request(
        FakeUi([["methane", "100", "100"]])
    )

    assert result.error_message is None
    assert result.request.component_names == ("methane",)
    assert result.request.mole_percents == (100.0,)
    assert result.request.weight_percents == (100.0,)
    assert result.request.basis == "Mol %"
    assert result.request.temperature_k == pytest.approx(298.15)
    assert result.request.pressure_pa == pytest.approx(101325.0)
    assert result.request.pressure_atm == pytest.approx(1.0)


def test_input_collector_reports_invalid_active_numeric_input():
    result = collect_property_calculation_request(
        FakeUi([["methane", "not-a-number", "100"]])
    )

    assert result.request is None
    assert result.error_message == (
        "Error: Invalid numeric input in composition table."
    )


def test_input_collector_preserves_inactive_numeric_tolerance():
    result = collect_property_calculation_request(
        FakeUi([["methane", "100", "not-a-number"]])
    )

    assert result.error_message is None
    assert result.request.weight_percents == (0.0,)


def test_input_collector_reports_no_components():
    result = collect_property_calculation_request(FakeUi([]))

    assert result.request is None
    assert result.error_message == "Error: No components selected."


def test_input_collector_reports_invalid_temperature_and_pressure():
    invalid_temperature = collect_property_calculation_request(
        FakeUi([["methane", "100", "100"]], temperature="bad")
    )
    invalid_pressure = collect_property_calculation_request(
        FakeUi([["methane", "100", "100"]], pressure="0 atm")
    )

    assert invalid_temperature.error_message == (
        "Error: Invalid Temperature selection."
    )
    assert invalid_pressure.error_message == (
        "Error: Invalid Pressure selection."
    )
