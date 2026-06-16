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
