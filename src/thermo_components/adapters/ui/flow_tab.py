"""Qt controller for Flow-tab setup, signals, and conversion rendering."""

from collections.abc import Callable

from PyQt6.QtCore import Qt

from thermo_components.application.dto import FlowConversionRequest
from thermo_components.domain.flow_units import FLOW_UNIT_ORDER


DensityProvider = Callable[[], float | None]


class FlowTabController:
    """Manage Flow-tab widgets while delegating conversion rules to a use case."""

    def __init__(
        self,
        ui,
        convert_flow_use_case,
        *,
        normal_density_provider: DensityProvider,
        standard_density_provider: DensityProvider,
    ):
        self.ui = ui
        self.convert_flow_use_case = convert_flow_use_case
        self.normal_density_provider = normal_density_provider
        self.standard_density_provider = standard_density_provider

    def setup(self) -> None:
        """Initialize Flow-tab widgets from the current UI definition."""
        if not self.has_required_widgets():
            return

        self.ui.comboBox_select_units.clear()
        self.ui.comboBox_select_units.addItems(FLOW_UNIT_ORDER)
        self.ui.comboBox_select_desired_units.clear()
        self.ui.comboBox_select_desired_units.addItems(FLOW_UNIT_ORDER)
        self.ui.comboBox_select_units.setCurrentText("kg/h")
        self.ui.comboBox_select_desired_units.setCurrentText("t/d")
        self.configure_copyable_output(self.ui.lineEdit_result)
        self.configure_copyable_output(self.ui.lineEdit_in_desired_units)
        self.update_conversion()

    def connect_signals(self) -> None:
        """Connect Flow-tab input widgets to conversion rendering."""
        if not self.has_required_widgets():
            return

        self.ui.lineEdit_enter_flow.textChanged.connect(self.update_conversion)
        self.ui.comboBox_select_units.currentIndexChanged.connect(
            self.update_conversion
        )
        self.ui.comboBox_select_desired_units.currentIndexChanged.connect(
            self.update_conversion
        )

    def update_conversion(self) -> None:
        """Render the latest flow conversion result or validation message."""
        if not self.has_required_widgets():
            return

        from_unit = self.ui.comboBox_select_units.currentText().strip()
        to_unit = self.ui.comboBox_select_desired_units.currentText().strip()
        self.ui.lineEdit_in_desired_units.setText(to_unit)

        if not from_unit or not to_unit:
            self.ui.lineEdit_result.clear()
            return

        try:
            response = self.convert_flow_use_case.execute(
                FlowConversionRequest(
                    input_text=self.ui.lineEdit_enter_flow.text(),
                    from_unit=from_unit,
                    to_unit=to_unit,
                    normal_density_kg_m3=self.normal_density_provider(),
                    standard_density_kg_m3=self.standard_density_provider(),
                )
            )
        except ValueError as exc:
            self.ui.lineEdit_result.setText(str(exc))
            return

        self.ui.lineEdit_result.setText(response.display_value)

    def has_required_widgets(self) -> bool:
        """Return whether this UI contains the expected Flow-tab widgets."""
        return all(
            hasattr(self.ui, widget_name)
            for widget_name in (
                "comboBox_select_units",
                "comboBox_select_desired_units",
                "lineEdit_enter_flow",
                "lineEdit_result",
                "lineEdit_in_desired_units",
            )
        )

    @staticmethod
    def configure_copyable_output(line_edit) -> None:
        """Keep outputs read-only while preserving selection and copy behavior."""
        line_edit.setReadOnly(True)
        line_edit.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        line_edit.setContextMenuPolicy(Qt.ContextMenuPolicy.DefaultContextMenu)
        line_edit.setCursor(Qt.CursorShape.IBeamCursor)
        line_edit.setDragEnabled(True)
