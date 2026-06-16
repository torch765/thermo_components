"""Qt worker bridge for executing application use cases off the UI thread."""

from PyQt6.QtCore import QObject, pyqtSignal


class CalculationWorker(QObject):
    """Execute one application request and emit Qt-friendly signals."""

    result = pyqtSignal(object)
    error = pyqtSignal(str)
    finished = pyqtSignal()

    def __init__(self, use_case, request):
        super().__init__()
        self.use_case = use_case
        self.request = request

    def run(self):
        try:
            self.result.emit(self.use_case.execute(self.request))
        except ValueError as exc:
            self.error.emit(f"Error: {exc}")
        except Exception as exc:
            import traceback

            self.error.emit(
                f"Worker Exception: {exc}\n{traceback.format_exc()}"
            )
        finally:
            self.finished.emit()
