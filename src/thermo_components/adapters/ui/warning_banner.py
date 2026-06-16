"""Qt controller for the persistent thermo warning banner."""

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QLabel


class ThermoWarningBannerController:
    """Manage the warning banner and related results-area geometry."""

    def __init__(self, ui):
        self.ui = ui
        self.label = QLabel(ui.tab)
        self.label.setObjectName("thermoWarningLabel")
        self.label.setWordWrap(True)
        self.label.setTextFormat(Qt.TextFormat.PlainText)
        self.label.setAlignment(
            Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop
        )
        self.label.setMargin(8)
        self.label.setStyleSheet(
            "QLabel#thermoWarningLabel {"
            "background-color: rgb(255, 244, 204);"
            "color: rgb(122, 63, 0);"
            "border: 1px solid rgb(230, 184, 0);"
            "border-radius: 4px;"
            "}"
        )
        self.label.hide()

        self._results_label_base_geometry = ui.results_label.geometry()
        self._results_list_base_geometry = ui.results_list.geometry()
        self._progress_bar_base_geometry = (
            ui.progressBar.geometry() if hasattr(ui, "progressBar") else None
        )
        self._warning_banner_y = ui.groupBox.geometry().bottom() + 8

    def set_messages(self, warning_messages) -> None:
        cleaned_messages = [
            str(message).strip()
            for message in warning_messages
            if str(message).strip()
        ]
        if not cleaned_messages:
            self.hide()
            return

        warning_text = "\n".join(cleaned_messages)
        banner_width = self._results_list_base_geometry.width()
        text_rect = self.label.fontMetrics().boundingRect(
            0,
            0,
            banner_width - 16,
            1000,
            int(Qt.TextFlag.TextWordWrap),
            warning_text,
        )
        banner_height = max(44, text_rect.height() + 16)
        self.label.setText(warning_text)
        self.label.setGeometry(
            self._results_list_base_geometry.x(),
            self._warning_banner_y,
            banner_width,
            banner_height,
        )
        self.label.show()
        self.label.raise_()

        results_label_y = self._warning_banner_y + banner_height + 8
        self.ui.results_label.setGeometry(
            self._results_label_base_geometry.x(),
            results_label_y,
            self._results_label_base_geometry.width(),
            self._results_label_base_geometry.height(),
        )

        list_y_offset = (
            self._results_list_base_geometry.y()
            - self._results_label_base_geometry.y()
        )
        results_list_y = results_label_y + list_y_offset
        base_results_bottom = (
            self._results_list_base_geometry.y()
            + self._results_list_base_geometry.height()
        )
        results_list_height = max(100, base_results_bottom - results_list_y)
        self.ui.results_list.setGeometry(
            self._results_list_base_geometry.x(),
            results_list_y,
            self._results_list_base_geometry.width(),
            results_list_height,
        )

        if self._progress_bar_base_geometry is not None:
            self.ui.progressBar.setGeometry(
                self._progress_bar_base_geometry.x(),
                base_results_bottom + 8,
                self._progress_bar_base_geometry.width(),
                self._progress_bar_base_geometry.height(),
            )

    def hide(self) -> None:
        self.label.hide()
        self.ui.results_label.setGeometry(self._results_label_base_geometry)
        self.ui.results_list.setGeometry(self._results_list_base_geometry)
        if self._progress_bar_base_geometry is not None:
            self.ui.progressBar.setGeometry(
                self._progress_bar_base_geometry
            )
