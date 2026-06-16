"""Qt controller for composition-table rendering and total validation."""

from collections.abc import Iterable
from dataclasses import dataclass
from enum import Enum

from PyQt6.QtCore import Qt
from PyQt6.QtGui import QColor, QPalette
from PyQt6.QtWidgets import QHeaderView, QTableWidgetItem


COMPONENT_COLUMN = 0
MOLE_PERCENT_COLUMN = 1
WEIGHT_PERCENT_COLUMN = 2
NORMALIZATION_SUCCESS_MESSAGE = "Composition normalized to 100%."


@dataclass(frozen=True)
class CompositionTotalState:
    total_percent: float
    valid_input: bool
    is_total_ok: bool


class ComponentAddStatus(str, Enum):
    ADDED = "added"
    DUPLICATE = "duplicate"
    EMPTY_SELECTION = "empty_selection"


class ComponentRemoveStatus(str, Enum):
    REMOVED = "removed"
    NO_COMPONENTS = "no_components"


@dataclass(frozen=True)
class ComponentAddResult:
    status: ComponentAddStatus
    component_name: str = ""


@dataclass(frozen=True)
class ComponentRemoveResult:
    status: ComponentRemoveStatus
    component_names: tuple[str, ...] = ()


class CompositionTableController:
    """Manage composition-table setup, basis styling, and total validation."""

    def __init__(self, ui):
        self.ui = ui

    @staticmethod
    def parse_float_or_zero(text: str) -> float:
        """Parse float from a table cell; treat blanks/invalid values as zero."""
        try:
            return float(text) if text and str(text).strip() else 0.0
        except (TypeError, ValueError):
            return 0.0

    def setup_table(self) -> None:
        """Set up the composition table headers and initial total row."""
        table = self.ui.tableWidget
        table.setColumnCount(3)
        table.setHorizontalHeaderLabels(["Component", "Mol %", "Wt %"])
        table.setColumnWidth(COMPONENT_COLUMN, 150)
        table.setColumnWidth(MOLE_PERCENT_COLUMN, 80)
        table.setColumnWidth(WEIGHT_PERCENT_COLUMN, 80)
        header = table.horizontalHeader()
        header.setSectionResizeMode(
            COMPONENT_COLUMN,
            QHeaderView.ResizeMode.Interactive,
        )
        header.setSectionResizeMode(
            MOLE_PERCENT_COLUMN,
            QHeaderView.ResizeMode.Stretch,
        )
        header.setSectionResizeMode(
            WEIGHT_PERCENT_COLUMN,
            QHeaderView.ResizeMode.Stretch,
        )

        table.setRowCount(1)
        total_item = QTableWidgetItem("Total")
        total_item.setFlags(Qt.ItemFlag.ItemIsEnabled)
        table.setItem(0, COMPONENT_COLUMN, total_item)

        mole_total_item = QTableWidgetItem("0.00")
        mole_total_item.setFlags(Qt.ItemFlag.ItemIsEnabled)
        table.setItem(0, MOLE_PERCENT_COLUMN, mole_total_item)

        weight_total_item = QTableWidgetItem("0.00")
        weight_total_item.setFlags(Qt.ItemFlag.ItemIsEnabled)
        table.setItem(0, WEIGHT_PERCENT_COLUMN, weight_total_item)

    def active_inactive_columns(self) -> tuple[int, int]:
        if self.ui.radioButton_mol_percent.isChecked():
            return MOLE_PERCENT_COLUMN, WEIGHT_PERCENT_COLUMN
        return WEIGHT_PERCENT_COLUMN, MOLE_PERCENT_COLUMN

    def apply_basis_state(self) -> None:
        """Apply editability and colors for the active/inactive basis columns."""
        active_col, inactive_col = self.active_inactive_columns()
        table = self.ui.tableWidget

        table.blockSignals(True)
        try:
            row_count = table.rowCount()
            total_row = row_count - 1

            for row in range(row_count):
                active_item = table.item(row, active_col)
                inactive_item = table.item(row, inactive_col)

                if active_item is None:
                    active_item = QTableWidgetItem("")
                    table.setItem(row, active_col, active_item)
                if inactive_item is None:
                    inactive_item = QTableWidgetItem("")
                    table.setItem(row, inactive_col, inactive_item)

                active_item.setBackground(QColor("white"))
                inactive_item.setBackground(QColor("lightGray"))

                if row != total_row:
                    active_item.setFlags(
                        Qt.ItemFlag.ItemIsEnabled
                        | Qt.ItemFlag.ItemIsSelectable
                        | Qt.ItemFlag.ItemIsEditable
                    )
                    inactive_item.setFlags(
                        Qt.ItemFlag.ItemIsEnabled
                        | Qt.ItemFlag.ItemIsSelectable
                    )
                else:
                    active_item.setFlags(
                        Qt.ItemFlag.ItemIsEnabled
                        | Qt.ItemFlag.ItemIsSelectable
                    )
                    inactive_item.setFlags(
                        Qt.ItemFlag.ItemIsEnabled
                        | Qt.ItemFlag.ItemIsSelectable
                    )
        finally:
            table.blockSignals(False)

    def update_table_total(self) -> CompositionTotalState:
        """Sum the active composition column and update total validation state."""
        active_col, inactive_col = self.active_inactive_columns()
        table = self.ui.tableWidget
        row_count = table.rowCount()
        total_row_index = row_count - 1
        total_percent = 0.0
        valid_input = True

        for row in range(total_row_index):
            percent_item = table.item(row, active_col)
            if percent_item is None:
                continue

            text = percent_item.text()
            try:
                total_percent += float(text) if text.strip() else 0.0
            except ValueError:
                if text.strip():
                    valid_input = False

        total_item_active = table.item(total_row_index, active_col)
        if total_item_active is None:
            total_item_active = QTableWidgetItem()
            table.setItem(total_row_index, active_col, total_item_active)
        total_item_active.setText(f"{total_percent:.4f}")
        total_item_active.setFlags(
            Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable
        )

        total_item_inactive = table.item(total_row_index, inactive_col)
        if total_item_inactive is None:
            total_item_inactive = QTableWidgetItem()
            table.setItem(total_row_index, inactive_col, total_item_inactive)
        total_item_inactive.setText("")
        total_item_inactive.setFlags(
            Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable
        )

        is_total_ok = abs(total_percent - 100.0) < 1e-4
        if is_total_ok and valid_input:
            total_item_active.setForeground(QColor("black"))
            self.ui.go_button.setEnabled(True)
        else:
            total_item_active.setForeground(QColor("red"))
            self.ui.go_button.setEnabled(False)

        return CompositionTotalState(
            total_percent=total_percent,
            valid_input=valid_input,
            is_total_ok=is_total_ok,
        )

    def read_active_percentages_for_normalization(self) -> tuple[float, ...]:
        """Read component percentages from the active basis column."""
        active_col, _ = self.active_inactive_columns()
        table = self.ui.tableWidget
        total_row = table.rowCount() - 1
        values: list[float] = []

        for row in range(total_row):
            item = table.item(row, active_col)
            try:
                value = (
                    float(item.text())
                    if item and item.text().strip()
                    else 0.0
                )
            except ValueError:
                value = 0.0
            values.append(value)

        return tuple(values)

    def write_normalized_active_percentages(
        self,
        percentages: Iterable[float],
    ) -> None:
        """Write normalized percentages into the active basis column."""
        active_col, _ = self.active_inactive_columns()
        table = self.ui.tableWidget

        table.blockSignals(True)
        try:
            for row, value in enumerate(percentages):
                item = table.item(row, active_col)
                if item is not None:
                    item.setText(f"{value:.4f}")
        finally:
            table.blockSignals(False)

    def show_normalization_success(self) -> None:
        if hasattr(self.ui, "results_list"):
            self.ui.results_list.addItem(NORMALIZATION_SUCCESS_MESSAGE)

    def add_selected_component(self) -> ComponentAddResult:
        """Add the selected component to the list and table widgets."""
        component_name = self.ui.comboBox_select_components.currentText()
        if not component_name.strip():
            return ComponentAddResult(ComponentAddStatus.EMPTY_SELECTION)

        if self.has_component(component_name):
            self.mark_component_combo_duplicate()
            return ComponentAddResult(
                ComponentAddStatus.DUPLICATE,
                component_name,
            )

        self.reset_component_combo_palette()
        self.ui.selected_components_list.addItem(component_name)
        self.insert_component_row(component_name)
        return ComponentAddResult(ComponentAddStatus.ADDED, component_name)

    def remove_selected_components(self) -> ComponentRemoveResult:
        """Remove selected components, or the current/last component fallback."""
        selected_items = list(self.ui.selected_components_list.selectedItems())
        removed_component_names: list[str] = []

        if not selected_items:
            if self.ui.selected_components_list.count() == 0:
                return ComponentRemoveResult(ComponentRemoveStatus.NO_COMPONENTS)

            current_row = self.ui.selected_components_list.currentRow()
            if current_row < 0:
                current_row = self.ui.selected_components_list.count() - 1
            item_to_remove = self.ui.selected_components_list.takeItem(
                current_row
            )
            if item_to_remove is not None:
                removed_component_names.append(item_to_remove.text())
                self.remove_component_from_table(item_to_remove.text())
        else:
            for item in selected_items:
                row = self.ui.selected_components_list.row(item)
                component_name = item.text()
                self.ui.selected_components_list.takeItem(row)
                removed_component_names.append(component_name)
                self.remove_component_from_table(component_name)

        if self.ui.selected_components_list.count() == 0:
            self.reset_component_combo_palette()

        return ComponentRemoveResult(
            ComponentRemoveStatus.REMOVED,
            tuple(removed_component_names),
        )

    def clear_component_rows_and_selection(self) -> None:
        """Remove all component rows while preserving the total row."""
        self.ui.selected_components_list.clear()
        while self.ui.tableWidget.rowCount() > 1:
            self.ui.tableWidget.removeRow(0)

        self.reset_component_combo_palette()
        self.reset_component_selection()

    def has_component(self, component_name: str) -> bool:
        current_items = [
            self.ui.selected_components_list.item(index).text()
            for index in range(self.ui.selected_components_list.count())
        ]
        return component_name in current_items

    def insert_component_row(self, component_name: str) -> None:
        total_row_index = self.ui.tableWidget.rowCount() - 1
        self.ui.tableWidget.insertRow(total_row_index)

        name_item = QTableWidgetItem(component_name)
        name_item.setFlags(name_item.flags() & ~Qt.ItemFlag.ItemIsEditable)
        self.ui.tableWidget.setItem(
            total_row_index,
            COMPONENT_COLUMN,
            name_item,
        )

        self.ui.tableWidget.setItem(
            total_row_index,
            MOLE_PERCENT_COLUMN,
            QTableWidgetItem(""),
        )
        self.ui.tableWidget.setItem(
            total_row_index,
            WEIGHT_PERCENT_COLUMN,
            QTableWidgetItem(""),
        )

    def remove_component_from_table(self, component_name: str) -> None:
        for row_index in range(self.ui.tableWidget.rowCount() - 2, -1, -1):
            table_item = self.ui.tableWidget.item(row_index, COMPONENT_COLUMN)
            if table_item and table_item.text() == component_name:
                self.ui.tableWidget.removeRow(row_index)
                break

    def mark_component_combo_duplicate(self) -> None:
        self._set_component_combo_base_color(QColor("red"))

    def reset_component_combo_palette(self) -> None:
        self._set_component_combo_base_color(QColor("white"))

    def reset_component_selection(self) -> None:
        self.ui.comboBox_select_components.setCurrentIndex(0)

    def _set_component_combo_base_color(self, color: QColor) -> None:
        palette = self.ui.comboBox_select_components.palette()
        palette.setColor(QPalette.ColorRole.Base, color)
        self.ui.comboBox_select_components.setPalette(palette)
