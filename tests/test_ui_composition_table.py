import pytest
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QColor

from density import MainWindow


@pytest.fixture
def main_window(qt_app):
    window = MainWindow(lhv_data={})
    yield window
    window.close()


def _add_component(window, component_name: str) -> None:
    index = window.ui.comboBox_select_components.findText(component_name)
    assert index >= 0
    window.ui.comboBox_select_components.setCurrentIndex(index)


def _is_editable(item) -> bool:
    return bool(item.flags() & Qt.ItemFlag.ItemIsEditable)


def _selected_component_names(window) -> list[str]:
    selected_components = window.ui.selected_components_list
    return [
        selected_components.item(index).text()
        for index in range(selected_components.count())
    ]


def _table_component_names(window) -> list[str]:
    table = window.ui.tableWidget
    return [
        table.item(row, 0).text()
        for row in range(table.rowCount() - 1)
    ]


def test_composition_table_setup_keeps_total_row(main_window):
    table = main_window.ui.tableWidget

    assert table.columnCount() == 3
    assert table.horizontalHeaderItem(0).text() == "Component"
    assert table.horizontalHeaderItem(1).text() == "Mol %"
    assert table.horizontalHeaderItem(2).text() == "Wt %"
    assert table.rowCount() == 1
    assert table.item(0, 0).text() == "Total"
    assert not _is_editable(table.item(0, 0))
    assert not _is_editable(table.item(0, 1))
    assert not _is_editable(table.item(0, 2))


def test_wt_basis_makes_weight_column_editable(main_window):
    _add_component(main_window, "methane")

    main_window.ui.radioButton_wt_percent.setChecked(True)

    table = main_window.ui.tableWidget
    mol_item = table.item(0, 1)
    wt_item = table.item(0, 2)
    total_mol_item = table.item(1, 1)
    total_wt_item = table.item(1, 2)

    assert not _is_editable(mol_item)
    assert _is_editable(wt_item)
    assert mol_item.background().color() == QColor("lightGray")
    assert wt_item.background().color() == QColor("white")
    assert not _is_editable(total_mol_item)
    assert not _is_editable(total_wt_item)


def test_total_validation_updates_go_button(main_window):
    _add_component(main_window, "methane")
    table = main_window.ui.tableWidget

    table.item(0, 1).setText("bad")
    invalid_state = main_window.update_table_total()

    assert invalid_state.total_percent == pytest.approx(0.0)
    assert not invalid_state.valid_input
    assert not invalid_state.is_total_ok
    assert table.item(1, 1).text() == "0.0000"
    assert not main_window.ui.go_button.isEnabled()

    table.item(0, 1).setText("100")
    valid_state = main_window.update_table_total()

    assert valid_state.total_percent == pytest.approx(100.0)
    assert valid_state.valid_input
    assert valid_state.is_total_ok
    assert table.item(1, 1).text() == "100.0000"
    assert table.item(1, 2).text() == ""
    assert main_window.ui.go_button.isEnabled()


def test_add_component_updates_list_and_table(main_window):
    _add_component(main_window, "methane")

    table = main_window.ui.tableWidget

    assert _selected_component_names(main_window) == ["methane"]
    assert _table_component_names(main_window) == ["methane"]
    assert table.rowCount() == 2
    assert table.item(0, 0).text() == "methane"
    assert not _is_editable(table.item(0, 0))
    assert _is_editable(table.item(0, 1))
    assert not _is_editable(table.item(0, 2))
    assert main_window.ui.comboBox_select_components.currentIndex() == 0


def test_duplicate_component_keeps_existing_rows_and_warns(
    main_window,
    monkeypatch,
):
    warnings = []
    monkeypatch.setattr(
        "density.QMessageBox.warning",
        lambda parent, title, message: warnings.append((title, message)),
    )

    _add_component(main_window, "methane")
    _add_component(main_window, "methane")

    assert _selected_component_names(main_window) == ["methane"]
    assert _table_component_names(main_window) == ["methane"]
    assert warnings == [
        ("Duplicate", "Component 'methane' is already selected.")
    ]
    assert main_window.ui.comboBox_select_components.currentIndex() == 0


def test_remove_selected_component_updates_list_and_table(main_window):
    _add_component(main_window, "methane")
    _add_component(main_window, "ethane")

    main_window.ui.selected_components_list.setCurrentRow(0)
    main_window.remove_component()

    assert _selected_component_names(main_window) == ["ethane"]
    assert _table_component_names(main_window) == ["ethane"]
    assert main_window.ui.tableWidget.rowCount() == 2


def test_remove_without_components_shows_information(
    main_window,
    monkeypatch,
):
    messages = []
    monkeypatch.setattr(
        "density.QMessageBox.information",
        lambda parent, title, message: messages.append((title, message)),
    )

    main_window.remove_component()

    assert messages == [("Remove", "No components to remove.")]
    assert main_window.ui.selected_components_list.count() == 0
    assert main_window.ui.tableWidget.rowCount() == 1


def test_clear_all_removes_component_rows_and_resets_selection(main_window):
    _add_component(main_window, "methane")
    _add_component(main_window, "ethane")
    main_window.ui.comboBox_select_temperature.setCurrentText("25 °C")
    main_window.ui.comboBox_select_pressure.setCurrentText("5 atm")
    main_window.ui.radioButton_wt_percent.setChecked(True)

    main_window.clear_all()

    assert main_window.ui.selected_components_list.count() == 0
    assert main_window.ui.tableWidget.rowCount() == 1
    assert main_window.ui.tableWidget.item(0, 0).text() == "Total"
    assert main_window.ui.comboBox_select_components.currentIndex() == 0
    assert main_window.ui.comboBox_select_temperature.currentText() == "0 °C"
    assert main_window.ui.comboBox_select_pressure.currentText() == "1 atm"
    assert main_window.ui.radioButton_mol_percent.isChecked()


def test_normalize_composition_updates_active_column_and_success_message(
    main_window,
):
    _add_component(main_window, "methane")
    _add_component(main_window, "ethane")
    table = main_window.ui.tableWidget

    table.item(0, 1).setText("40")
    table.item(1, 1).setText("40")

    main_window.normalize_composition()

    assert float(table.item(0, 1).text()) == pytest.approx(50.0)
    assert float(table.item(1, 1).text()) == pytest.approx(50.0)
    assert table.item(2, 1).text() == "100.0000"
    assert main_window.ui.results_list.item(
        main_window.ui.results_list.count() - 1
    ).text() == "Composition normalized to 100%."


def test_normalize_composition_zero_total_warns_without_success_message(
    main_window,
    monkeypatch,
):
    warnings = []
    monkeypatch.setattr(
        "density.QMessageBox.warning",
        lambda parent, title, message: warnings.append((title, message)),
    )
    _add_component(main_window, "methane")

    main_window.normalize_composition()

    assert warnings == [("Normalize", "Cannot normalize: total is zero.")]
    assert main_window.ui.results_list.count() == 0
