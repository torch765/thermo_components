"""Qt controller for property calculation workflow and result lifecycle."""

from collections.abc import Callable

from PyQt6.QtCore import QThread

from thermo_components.adapters.ui.input_collection import (
    collect_property_calculation_request,
)
from thermo_components.adapters.ui.presenters import build_result_list_items
from thermo_components.adapters.ui.qt_worker import CalculationWorker
from thermo_components.application.dto import coerce_property_response


class QtCalculationWorkflowController:
    """Coordinate calculation input, worker threading, and result rendering."""

    def __init__(
        self,
        ui,
        calculate_properties_use_case,
        *,
        lhv_data_loaded_provider: Callable[[], bool],
        fallback_model_provider: Callable[[], str],
        invalidate_results: Callable[..., None],
        set_result_state: Callable[[object, dict], None],
        set_warning_messages: Callable[[list[str]], None],
        update_flow_conversion: Callable[[], None],
        animate_progress_to: Callable[[int, int], None],
        thread_factory: Callable[[], object] = QThread,
        worker_factory: Callable[[object, object], object] | None = None,
    ):
        self.ui = ui
        self.calculate_properties_use_case = calculate_properties_use_case
        self.lhv_data_loaded_provider = lhv_data_loaded_provider
        self.fallback_model_provider = fallback_model_provider
        self.invalidate_results = invalidate_results
        self.set_result_state = set_result_state
        self.set_warning_messages = set_warning_messages
        self.update_flow_conversion = update_flow_conversion
        self.animate_progress_to = animate_progress_to
        self.thread_factory = thread_factory
        self.worker_factory = worker_factory or self._build_worker

    def calculate_and_display(self) -> None:
        """Gather inputs, run calculations, and display results."""
        self.reset_progress()
        self.animate_progress_to(80, 1000)
        self.invalidate_results()

        collected_input = collect_property_calculation_request(self.ui)
        if collected_input.error_message:
            self.ui.results_list.addItem(collected_input.error_message)
            self.animate_progress_to(100, 500)
            return

        self.ui.go_button.setEnabled(False)
        if hasattr(self.ui, "printResultsButton"):
            self.ui.printResultsButton.setEnabled(False)

        self.worker_thread = self.thread_factory()
        self.worker = self.worker_factory(
            self.calculate_properties_use_case,
            collected_input.request,
        )
        self.worker.moveToThread(self.worker_thread)
        self.worker.result.connect(self.on_calculation_result)
        self.worker.error.connect(self.on_calculation_error)
        self.worker.finished.connect(self.on_calculation_finished)
        self.worker_thread.started.connect(self.worker.run)
        self.worker.finished.connect(self.worker_thread.quit)
        self.worker.finished.connect(self.worker.deleteLater)
        self.worker_thread.finished.connect(self.worker_thread.deleteLater)
        self.worker_thread.start()

    def on_calculation_result(self, result_data) -> None:
        """Render a successful calculation result into the Qt widgets."""
        response = coerce_property_response(result_data)
        legacy_result = response.to_legacy_dict()

        self.set_result_state(response, legacy_result)
        self.set_warning_messages(legacy_result.get("warnings") or [])
        self.update_flow_conversion()

        self.ui.results_list.clear()
        for item in build_result_list_items(
            response,
            lhv_data_loaded=self.lhv_data_loaded_provider(),
            fallback_model_display=self.fallback_model_provider(),
        ):
            self.ui.results_list.addItem(item)
        self.animate_progress_to(100, 500)

    def on_calculation_error(self, error_msg: str) -> None:
        """Clear stale results and render a calculation error."""
        self.invalidate_results(clear_visible_results=False)
        self.ui.results_list.clear()
        self.ui.results_list.addItem(error_msg)
        self.animate_progress_to(100, 500)

    def on_calculation_finished(self) -> None:
        """Restore buttons after worker completion."""
        self.ui.go_button.setEnabled(True)
        if hasattr(self.ui, "printResultsButton"):
            self.ui.printResultsButton.setEnabled(True)

    def reset_progress(self) -> None:
        if hasattr(self.ui, "progressBar"):
            self.ui.progressBar.setValue(0)

    @staticmethod
    def _build_worker(use_case, request):
        return CalculationWorker(use_case=use_case, request=request)
