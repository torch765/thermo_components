import pytest

from thermo_components.adapters.ui import QtCalculationWorkflowController
from thermo_components.adapters.ui.qt_main_window import MainWindow
from thermo_components.domain.thermo_routes import PRMIX_DEFAULT_ROUTE


class FakeSignal:
    def __init__(self):
        self.callbacks = []

    def connect(self, callback):
        self.callbacks.append(callback)

    def emit(self, *args):
        for callback in list(self.callbacks):
            callback(*args)


class FakeThread:
    def __init__(self):
        self.started = FakeSignal()
        self.finished = FakeSignal()
        self.start_called = False
        self.quit_called = False
        self.delete_later_called = False

    def start(self):
        self.start_called = True

    def quit(self):
        self.quit_called = True

    def deleteLater(self):
        self.delete_later_called = True


class FakeWorker:
    def __init__(self, use_case, request):
        self.use_case = use_case
        self.request = request
        self.result = FakeSignal()
        self.error = FakeSignal()
        self.finished = FakeSignal()
        self.thread = None
        self.delete_later_called = False

    def moveToThread(self, thread):
        self.thread = thread

    def run(self):
        pass

    def deleteLater(self):
        self.delete_later_called = True


@pytest.fixture
def main_window(qt_app):
    window = MainWindow(lhv_data={"methane": 35.8})
    yield window
    window.close()


def _select_component(window, component_name: str) -> None:
    index = window.ui.comboBox_select_components.findText(component_name)
    assert index >= 0
    window.ui.comboBox_select_components.setCurrentIndex(index)


def _sample_result_mapping(**overrides):
    data = {
        "mw": 16.04,
        "basis": "Mol %",
        "model_display": "PRMIX",
        "eos": "PRMIX",
        "phase": "Vapor",
        "density_result": 0.6571965,
        "density_error": None,
        "density_normal_kg_m3": 0.716,
        "density_normal_phase": "Vapor",
        "density_normal_error": None,
        "density_standard_kg_m3": 0.68,
        "density_standard_phase": "Vapor",
        "density_standard_error": None,
        "pressure_atm": 1.0,
        "bubble_point": -161.5,
        "bp_error": None,
        "thermo_route": PRMIX_DEFAULT_ROUTE,
        "mixture_lhv": 35.8,
        "missing_lhv": [],
        "warnings": [],
    }
    data.update(overrides)
    return data


def test_calculation_workflow_controller_handles_invalid_input(
    main_window,
    monkeypatch,
):
    progress_calls = []
    monkeypatch.setattr(
        main_window,
        "animate_progress_to",
        lambda target, duration_ms=500: progress_calls.append(
            (target, duration_ms)
        ),
    )

    main_window.calculation_workflow_controller.calculate_and_display()

    assert main_window.ui.results_list.item(0).text() == (
        "Error: No components selected."
    )
    assert progress_calls == [(80, 1000), (100, 500)]
    assert not hasattr(main_window.calculation_workflow_controller, "worker")
    assert not hasattr(
        main_window.calculation_workflow_controller,
        "worker_thread",
    )


def test_calculation_workflow_controller_starts_worker_for_valid_input(
    main_window,
    monkeypatch,
):
    progress_calls = []
    fake_threads = []
    fake_workers = []

    def build_thread():
        thread = FakeThread()
        fake_threads.append(thread)
        return thread

    def build_worker(use_case, request):
        worker = FakeWorker(use_case, request)
        fake_workers.append(worker)
        return worker

    controller = main_window.calculation_workflow_controller
    controller.thread_factory = build_thread
    controller.worker_factory = build_worker
    monkeypatch.setattr(
        main_window,
        "animate_progress_to",
        lambda target, duration_ms=500: progress_calls.append(
            (target, duration_ms)
        ),
    )

    _select_component(main_window, "methane")
    main_window.ui.tableWidget.item(0, 1).setText("100")

    controller.calculate_and_display()

    assert progress_calls == [(80, 1000)]
    assert not main_window.ui.go_button.isEnabled()
    assert not main_window.ui.printResultsButton.isEnabled()
    assert fake_threads[0].start_called
    assert fake_workers[0].request.component_names == ("methane",)
    assert fake_workers[0].thread is fake_threads[0]
    assert fake_threads[0].started.callbacks == [fake_workers[0].run]
    assert controller.worker is fake_workers[0]
    assert controller.worker_thread is fake_threads[0]


def test_calculation_workflow_controller_renders_result(
    main_window,
    monkeypatch,
):
    progress_calls = []
    flow_updates = []
    monkeypatch.setattr(
        main_window,
        "animate_progress_to",
        lambda target, duration_ms=500: progress_calls.append(
            (target, duration_ms)
        ),
    )
    monkeypatch.setattr(
        main_window.flow_tab_controller,
        "update_conversion",
        lambda: flow_updates.append(True),
    )

    main_window.calculation_workflow_controller.on_calculation_result(
        _sample_result_mapping()
    )

    assert main_window.ui.results_list.item(0).text() == (
        "Calculating for Mol % basis..."
    )
    assert main_window.density_actual_kg_m3 == 0.6571965
    assert main_window.density_normal_kg_m3 == 0.716
    assert main_window.density_standard_kg_m3 == 0.68
    assert flow_updates == [True]
    assert progress_calls == [(100, 500)]


def test_calculation_workflow_controller_renders_error_and_finish_state(
    main_window,
    monkeypatch,
):
    progress_calls = []
    monkeypatch.setattr(
        main_window,
        "animate_progress_to",
        lambda target, duration_ms=500: progress_calls.append(
            (target, duration_ms)
        ),
    )
    main_window.last_result_data = _sample_result_mapping()
    main_window.density_actual_kg_m3 = 1.0
    main_window.density_normal_kg_m3 = 2.0
    main_window.density_standard_kg_m3 = 3.0
    main_window.ui.results_list.addItem("old result")
    main_window.ui.go_button.setEnabled(False)
    main_window.ui.printResultsButton.setEnabled(False)

    controller = main_window.calculation_workflow_controller
    controller.on_calculation_error("calculation failed")
    controller.on_calculation_finished()

    assert main_window.last_result_data is None
    assert main_window.density_actual_kg_m3 is None
    assert main_window.density_normal_kg_m3 is None
    assert main_window.density_standard_kg_m3 is None
    assert main_window.ui.results_list.count() == 1
    assert main_window.ui.results_list.item(0).text() == "calculation failed"
    assert progress_calls == [(100, 500)]
    assert main_window.ui.go_button.isEnabled()
    assert main_window.ui.printResultsButton.isEnabled()


def test_main_window_accepts_calculation_workflow_controller(main_window):
    assert isinstance(
        main_window.calculation_workflow_controller,
        QtCalculationWorkflowController,
    )
