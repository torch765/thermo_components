"""Composition normalization and inactive-basis workflows."""

from thermo_components.application.dto import (
    DeriveCompositionRequest,
    DeriveCompositionResponse,
    NormalizeCompositionRequest,
    NormalizeCompositionResponse,
)
from thermo_components.domain.composition import (
    derive_inactive_percentages,
    normalize_percentages,
)


class NormalizeCompositionUseCase:
    def execute(
        self,
        request: NormalizeCompositionRequest,
    ) -> NormalizeCompositionResponse:
        values = normalize_percentages(
            request.percentages,
            precision=request.precision,
        )
        return NormalizeCompositionResponse(percentages=tuple(values))


class DeriveCompositionUseCase:
    def execute(
        self,
        request: DeriveCompositionRequest,
    ) -> DeriveCompositionResponse:
        values = derive_inactive_percentages(
            request.component_names,
            request.active_percentages,
            request.basis,
        )
        return DeriveCompositionResponse(percentages=tuple(values))
