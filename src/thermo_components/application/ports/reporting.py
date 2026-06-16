"""Ports for exporting prepared reports."""

from pathlib import Path
from typing import Protocol, runtime_checkable

from thermo_components.application.dto import ReportExportRequest


@runtime_checkable
class ReportExporter(Protocol):
    """Write a prepared report projection to an external representation."""

    def export(self, request: ReportExportRequest) -> Path: ...
