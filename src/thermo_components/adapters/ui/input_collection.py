"""Collect calculation requests from the Qt widget tree."""

from dataclasses import dataclass

from thermo_components.application.dto import PropertyCalculationRequest


@dataclass(frozen=True)
class CalculationInputCollection:
    request: PropertyCalculationRequest | None = None
    error_message: str | None = None


def collect_property_calculation_request(ui) -> CalculationInputCollection:
    """Translate the current Qt widget state into an application request."""
    component_names = []
    mole_percents = []
    weight_percents = []
    row_count = ui.tableWidget.rowCount() - 1
    valid_input = True
    mol_basis_active = ui.radioButton_mol_percent.isChecked()
    wt_basis_active = ui.radioButton_wt_percent.isChecked()

    for row in range(row_count):
        component_item = ui.tableWidget.item(row, 0)
        mole_item = ui.tableWidget.item(row, 1)
        weight_item = ui.tableWidget.item(row, 2)
        if not component_item:
            continue

        component_names.append(component_item.text())
        mole_value, mole_valid = _parse_cell_float(mole_item)
        weight_value, weight_valid = _parse_cell_float(weight_item)
        mole_percents.append(mole_value)
        weight_percents.append(weight_value)

        if mol_basis_active and not mole_valid:
            valid_input = False
        if wt_basis_active and not weight_valid:
            valid_input = False

    if not valid_input:
        return CalculationInputCollection(
            error_message="Error: Invalid numeric input in composition table."
        )
    if not component_names:
        return CalculationInputCollection(
            error_message="Error: No components selected."
        )

    basis = "Mol %" if mol_basis_active else "Wt %"
    try:
        temperature_text = ui.comboBox_select_temperature.currentText()
        temperature_c = float(temperature_text.split("°")[0])
        temperature_k = temperature_c + 273.15
    except ValueError:
        return CalculationInputCollection(
            error_message="Error: Invalid Temperature selection."
        )

    try:
        pressure_text = ui.comboBox_select_pressure.currentText()
        pressure_atm = float(pressure_text.split(" ")[0])
        pressure_pa = pressure_atm * 101325.0
        if pressure_pa <= 0:
            raise ValueError("Pressure must be positive")
    except ValueError:
        return CalculationInputCollection(
            error_message="Error: Invalid Pressure selection."
        )

    return CalculationInputCollection(
        request=PropertyCalculationRequest.from_sequences(
            component_names=component_names,
            mole_percents=mole_percents,
            weight_percents=weight_percents,
            basis=basis,
            temperature_k=temperature_k,
            pressure_pa=pressure_pa,
            pressure_atm=pressure_atm,
        )
    )


def _parse_cell_float(item) -> tuple[float, bool]:
    text = item.text() if item is not None else ""
    if not text or not str(text).strip():
        return 0.0, True
    try:
        return float(text), True
    except ValueError:
        return 0.0, False
