"""Port for thermodynamic property calculations."""

from collections.abc import Mapping
from typing import Protocol, runtime_checkable


DensityValue = float | tuple[float, float] | None
DensityCalculation = tuple[DensityValue, str | None, str | None]
BubblePointCalculation = tuple[float | None, str | None]


@runtime_checkable
class ThermoPropertyGateway(Protocol):
    """Application-facing contract for a thermodynamic calculation engine."""

    eos: str

    def set_components(self, components: Mapping[str, float]) -> None: ...

    def set_eos(self, eos: str) -> None: ...

    def calculate_density_for_route(
        self,
        temperature_k: float,
        pressure_pa: float,
        route_id: str,
    ) -> DensityCalculation: ...

    def calculate_bubble_point_for_route(
        self,
        pressure_pa: float,
        route_id: str,
    ) -> BubblePointCalculation: ...
