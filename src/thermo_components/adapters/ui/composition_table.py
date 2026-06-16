"""Qt controller for composition-table rendering and total validation."""

from dataclasses import dataclass

from PyQt6.QtCore import Qt
from PyQt6.QtGui import QColor
from PyQt6.QtWidgets import QHeaderView, QTableWidgetItem


COMPONENT_COLUMN = 0
MOLE_PERCENT_COLUMN = 1
WEIGHT_PERCENT_COLUMN = 2


@dataclass(frozen=True)
class CompositionTotalState:
    total_percent: float
    valid_input: bool
    is_total_ok: bool


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
