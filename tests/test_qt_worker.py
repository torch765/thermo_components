from thermo_components.adapters.ui import CalculationWorker


class StubUseCase:
    def __init__(self, result=None, error=None):
        self.result = result
        self.error = error
        self.requests = []

    def execute(self, request):
        self.requests.append(request)
        if self.error is not None:
            raise self.error
        return self.result


def test_calculation_worker_emits_use_case_response_and_finished(qt_app):
    use_case = StubUseCase(result=object())
    request = object()
    results = []
    errors = []
    finished = []
    worker = CalculationWorker(use_case, request)
    worker.result.connect(results.append)
    worker.error.connect(errors.append)
    worker.finished.connect(lambda: finished.append(True))

    worker.run()

    assert use_case.requests == [request]
    assert results == [use_case.result]
    assert errors == []
    assert finished == [True]


def test_calculation_worker_maps_validation_error_and_finishes(qt_app):
    use_case = StubUseCase(error=ValueError("Invalid Mol % input."))
    errors = []
    finished = []
    worker = CalculationWorker(use_case, object())
    worker.error.connect(errors.append)
    worker.finished.connect(lambda: finished.append(True))

    worker.run()

    assert errors == ["Error: Invalid Mol % input."]
    assert finished == [True]
