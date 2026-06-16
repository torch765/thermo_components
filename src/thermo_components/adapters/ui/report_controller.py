"""Qt controller for report export actions and user messaging."""

import os
import sys
from collections.abc import Callable
from datetime import datetime
from pathlib import Path

from PyQt6.QtWidgets import QMessageBox

from thermo_components.application.dto import ReportExportRequest


class QtReportExportController:
    """Coordinate report export paths, dialogs, and optional file opening."""

    def __init__(
        self,
        parent,
        report_exporter,
        *,
        source_base_dir: Path | str | None = None,
        export_base_dir_provider: Callable[[], Path | str] | None = None,
        clock: Callable[[], datetime] | None = None,
        file_opener: Callable[[Path], None] | None = None,
    ):
        self.parent = parent
        self.report_exporter = report_exporter
        self.source_base_dir = (
            Path(source_base_dir)
            if source_base_dir is not None
            else Path.cwd()
        )
        self.export_base_dir_provider = export_base_dir_provider
        self.clock = clock or datetime.now
        self.file_opener = file_opener

    def get_export_base_dir(self) -> Path:
        """Return the directory where generated reports should be saved."""
        if self.export_base_dir_provider is not None:
            return Path(self.export_base_dir_provider())
        if getattr(sys, "frozen", False):
            return Path(sys.executable).parent
        return self.source_base_dir

    def build_report_filename(
        self,
        timestamp: datetime | None = None,
    ) -> str:
        """Build the timestamped Excel report filename."""
        timestamp = timestamp or self.clock()
        return f"Thermo_Report_{timestamp.strftime('%Y%m%d_%H%M%S')}.xlsx"

    def export_latest_report(
        self,
        *,
        has_results: bool,
        build_request: Callable[[Path, datetime], ReportExportRequest],
    ) -> str | None:
        """Export the current report request and show Qt user feedback."""
        if not has_results:
            QMessageBox.information(
                self.parent,
                "Print Results",
                "No calculation results available. Please click Go first.",
            )
            return None

        timestamp = self.clock()
        report_path = self.get_export_base_dir() / self.build_report_filename(
            timestamp
        )
        export_request = build_request(report_path, timestamp)

        try:
            exported_path = self.report_exporter.export(export_request)
        except ImportError:
            QMessageBox.warning(
                self.parent,
                "Print Results",
                "Excel export requires the 'openpyxl' package.\n"
                "Install it in this environment and try again.\n\n"
                "Command:\npython -m pip install openpyxl",
            )
            return None
        except Exception as exc:
            QMessageBox.critical(
                self.parent,
                "Print Results",
                f"Failed to save Excel report:\n{exc}",
            )
            return None

        QMessageBox.information(
            self.parent,
            "Print Results",
            f"Excel report saved to:\n{exported_path}",
        )
        self.open_exported_file(Path(exported_path))
        return str(exported_path)

    def open_exported_file(self, exported_path: Path) -> None:
        if self.file_opener is not None:
            self.file_opener(exported_path)
            return

        if hasattr(os, "startfile"):
            try:
                os.startfile(exported_path)
            except OSError:
                pass
