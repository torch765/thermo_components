import pytest
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QMainWindow

from gui import Ui_Dialog
from thermo_components.adapters.ui import FlowTabController
from thermo_components.application.dto import (
    FlowConversionRequest,
    FlowConversionResponse,
)
from thermo_components.domain.flow_units import FLOW_UNIT_ORDER


class RecordingFlowUseCase:
    def __init__(self, response=None, error_message: str | None = None):
        self.response = response or FlowConversionResponse(
            value=None,
            display_value="",
        )
        self.error_message = error_message
        self.requests: list[FlowConversionRequest] = []

    def execute(self, request: FlowConversionRequest):
        self.requests.append(request)
        if self.error_message is not None:
            raise ValueError(self.error_message)
        return self.response


@pytest.fixture
def flow_ui(qt_app):
    window = QMainWindow()
    ui = Ui_Dialog()
    ui.setupUi(window)
    yield ui
    window.close()


def _build_controller(
    ui,
    use_case,
    *,
    normal_density=0.716,
    standard_density=0.68,
):
    return FlowTabController(
        ui,
        use_case,
        normal_density_provider=lambda: normal_density,
        standard_density_provider=lambda: standard_density,
    )


def test_flow_tab_controller_setup_populates_units_and_copyable_outputs(flow_ui):
    use_case = RecordingFlowUseCase()
    controller = _build_controller(flow_ui, use_case)

    controller.setup()

    assert [
        flow_ui.comboBox_select_units.itemText(index)
        for index in range(flow_ui.comboBox_select_units.count())
    ] == FLOW_UNIT_ORDER
    assert flow_ui.comboBox_select_units.currentText() == "kg/h"
    assert flow_ui.comboBox_select_desired_units.currentText() == "t/d"
    assert flow_ui.lineEdit_in_desired_units.text() == "t/d"
    assert flow_ui.lineEdit_result.isReadOnly()
    assert flow_ui.lineEdit_result.focusPolicy() == Qt.FocusPolicy.StrongFocus
    assert flow_ui.lineEdit_result.dragEnabled()


def test_flow_tab_controller_builds_request_and_renders_response(flow_ui):
    use_case = RecordingFlowUseCase(
        FlowConversionResponse(value=0.024, display_value="0.024")
    )
    controller = _build_controller(
        flow_ui,
        use_case,
        normal_density=1.2,
        standard_density=1.0,
    )
    controller.setup()
    use_case.requests.clear()

    flow_ui.lineEdit_enter_flow.setText("1")
    controller.update_conversion()

    assert use_case.requests == [
        FlowConversionRequest(
            input_text="1",
            from_unit="kg/h",
            to_unit="t/d",
            normal_density_kg_m3=1.2,
            standard_density_kg_m3=1.0,
        )
    ]
    assert flow_ui.lineEdit_result.text() == "0.024"
    assert flow_ui.lineEdit_in_desired_units.text() == "t/d"


def test_flow_tab_controller_renders_use_case_errors(flow_ui):
    use_case = RecordingFlowUseCase(error_message="Normal density required")
    controller = _build_controller(flow_ui, use_case)
    controller.setup()

    flow_ui.comboBox_select_units.setCurrentText("kg/h")
    flow_ui.comboBox_select_desired_units.setCurrentText("Nm3/h")
    flow_ui.lineEdit_enter_flow.setText("1")
    controller.update_conversion()

    assert flow_ui.lineEdit_result.text() == "Normal density required"


def test_flow_tab_controller_connects_widget_signals(flow_ui):
    use_case = RecordingFlowUseCase(
        FlowConversionResponse(value=24.0, display_value="24")
    )
    controller = _build_controller(flow_ui, use_case)
    controller.setup()
    controller.connect_signals()
    use_case.requests.clear()

    flow_ui.lineEdit_enter_flow.setText("2")

    assert use_case.requests[-1].input_text == "2"
    assert flow_ui.lineEdit_result.text() == "24"


def test_flow_tab_controller_ignores_missing_flow_widgets():
    use_case = RecordingFlowUseCase()
    controller = _build_controller(object(), use_case)

    controller.setup()
    controller.connect_signals()
    controller.update_conversion()

    assert use_case.requests == []
