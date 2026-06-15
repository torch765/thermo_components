"""Adapter between the application thermo port and the ``thermo`` package."""

from collections.abc import Mapping

from chemicals import iapws as chemicals_iapws
from thermo import (
    CEOSGas,
    CEOSLiquid,
    ChemicalConstantsPackage,
    FlashPureVLS,
    FlashVL,
    IAPWS95,
    PRMIX,
)

from thermo_components.application.ports import (
    BubblePointCalculation,
    DensityCalculation,
)
from thermo_components.domain.composition import (
    calculate_mixture_molecular_weight,
)
from thermo_components.domain.lhv import calculate_mixture_lhv
from thermo_components.domain.thermo_routes import (
    IAPWS95_TWO_PHASE_REL_TOL,
    PURE_WATER_ROUTE,
)


class ThermoGateway:
    """Provide route-aware thermodynamic calculations through ``thermo``."""

    def __init__(
        self,
        components: Mapping[str, float] | None = None,
        eos: str = "PRMIX",
    ):
        self.components = dict(components) if components is not None else {}
        self.eos = eos
        self.constants = None
        self.properties = None

    def set_components(self, components: Mapping[str, float]) -> None:
        self.components = dict(components)
        self.constants = None
        self.properties = None

    def set_eos(self, eos: str) -> None:
        self.eos = eos
        print(f"EOS set to: {self.eos}")

    def calculate_molecular_weight(self) -> float:
        return calculate_mixture_molecular_weight(self.components)

    def calculate_density_for_route(
        self,
        temperature_k: float,
        pressure_pa: float,
        route_id: str,
    ) -> DensityCalculation:
        """Dispatch density calculation to the selected thermo route."""
        if route_id == PURE_WATER_ROUTE:
            return self.calculate_pure_water_density_iapws95(
                temperature_k,
                pressure_pa,
            )
        return self.calculate_density(temperature_k, pressure_pa)

    def calculate_bubble_point_for_route(
        self,
        pressure_pa: float,
        route_id: str,
    ) -> BubblePointCalculation:
        """Dispatch saturation calculation to the selected thermo route."""
        if route_id == PURE_WATER_ROUTE:
            return self.calculate_pure_water_bubble_point_iapws95(pressure_pa)
        return self.calculate_bubble_point(pressure_pa)

    def calculate_pure_water_density_iapws95(
        self,
        temperature_k: float,
        pressure_pa: float,
    ) -> DensityCalculation:
        """Calculate pure-water density using the local IAPWS-95 implementation."""
        if temperature_k is None or pressure_pa is None:
            return None, None, "Missing temperature or pressure."
        if temperature_k <= 0 or pressure_pa <= 0:
            return None, None, "Temperature and pressure must be positive."

        try:
            phase_label = None
            critical_t = IAPWS95.Tc
            critical_p = IAPWS95.Pc

            if 235.0 <= temperature_k < critical_t and pressure_pa < critical_p:
                psat = chemicals_iapws.iapws95_Psat(temperature_k)
                relative_difference = abs(pressure_pa - psat) / max(psat, 1.0)
                if relative_difference <= IAPWS95_TWO_PHASE_REL_TOL:
                    density_liq = chemicals_iapws.iapws95_rhol_sat(
                        temperature_k
                    )
                    density_gas = chemicals_iapws.iapws95_rhog_sat(
                        temperature_k
                    )
                    return (density_liq, density_gas), "Two-Phase", None
                phase_label = "Liquid" if pressure_pa > psat else "Vapor"
            elif temperature_k >= critical_t and pressure_pa >= critical_p:
                phase_label = "Supercritical"
            elif temperature_k >= critical_t:
                phase_label = "Vapor"
            elif pressure_pa >= critical_p:
                phase_label = "Liquid"

            water_state = IAPWS95(T=temperature_k, P=pressure_pa, zs=[1.0])
            density_mass = water_state.rho_mass()

            if phase_label is None:
                phase_label = "Liquid" if density_mass >= 200.0 else "Vapor"

            return density_mass, phase_label, None
        except Exception as exc:
            return None, None, f"IAPWS-95 pure-water density failed: {exc}"

    def calculate_pure_water_bubble_point_iapws95(
        self,
        pressure_pa: float,
    ) -> BubblePointCalculation:
        """Calculate pure-water saturation temperature using IAPWS-95."""
        if pressure_pa is None or pressure_pa <= 0:
            return None, "Pressure must be positive."

        try:
            saturation_temperature_k = chemicals_iapws.iapws95_Tsat(
                pressure_pa
            )
        except Exception as exc:
            return None, f"IAPWS-95 saturation calculation failed: {exc}"

        return saturation_temperature_k - 273.15, None

    def calculate_density(
        self,
        temperature_k: float,
        pressure_pa: float,
    ) -> DensityCalculation:
        """Calculate density while preserving the existing phase logic."""
        if not self.components:
            return None, None, "No components selected."

        total = sum(self.components.values())
        if total == 0:
            return None, None, "Total fraction is zero."

        zs = [value / total for value in self.components.values()]
        component_names = list(self.components.keys())
        print(
            "Calculating density (original logic) for: "
            f"{component_names} at T={temperature_k:.2f} K, "
            f"P={pressure_pa:.1f} Pa with EOS={self.eos}"
        )

        try:
            if self.constants is None or self.properties is None:
                print("Fetching chemical constants...")
                self.constants, self.properties = (
                    ChemicalConstantsPackage.from_IDs(component_names)
                )
                print("Constants fetched.")

            if (
                None in self.constants.Tcs
                or None in self.constants.Pcs
                or None in self.constants.omegas
            ):
                missing_components = [
                    component_names[index]
                    for index, critical_temperature in enumerate(
                        self.constants.Tcs
                    )
                    if critical_temperature is None
                ]
                return (
                    None,
                    None,
                    "Missing critical properties for: "
                    + ", ".join(missing_components),
                )

            eos_kwargs = {
                "Pcs": self.constants.Pcs,
                "Tcs": self.constants.Tcs,
                "omegas": self.constants.omegas,
            }
            heat_capacity_gases = (
                self.properties.HeatCapacityGases
                if hasattr(self.properties, "HeatCapacityGases")
                else None
            )

            if len(component_names) == 1:
                gas = CEOSGas(
                    PRMIX,
                    eos_kwargs=eos_kwargs,
                    HeatCapacityGases=heat_capacity_gases,
                )
                liquid = CEOSLiquid(
                    PRMIX,
                    eos_kwargs=eos_kwargs,
                    HeatCapacityGases=heat_capacity_gases,
                )
                flasher = FlashPureVLS(
                    self.constants,
                    self.properties,
                    gas=gas,
                    liquids=[liquid],
                    solids=[],
                )
            else:
                gas = CEOSGas(
                    PRMIX,
                    HeatCapacityGases=heat_capacity_gases,
                    eos_kwargs=eos_kwargs,
                    zs=zs,
                )
                liquid = CEOSLiquid(
                    PRMIX,
                    HeatCapacityGases=heat_capacity_gases,
                    eos_kwargs=eos_kwargs,
                    zs=zs,
                )
                flasher = FlashVL(
                    self.constants,
                    self.properties,
                    liquid=liquid,
                    gas=gas,
                )

            print(f"Performing flash T={temperature_k}, P={pressure_pa}, zs={zs}")
            result = flasher.flash(T=temperature_k, P=pressure_pa, zs=zs)
            print(
                f"Flash result: Phases present = {len(result.phases)}, "
                f"PurePhase='{getattr(result, 'phase', 'N/A')}'"
            )
            if hasattr(result, "phases"):
                for index, phase in enumerate(result.phases):
                    print(
                        f"  Phase {index}: Liquid={phase.is_liquid}, "
                        f"Gas={phase.is_gas}, V={phase.V()}"
                    )

            molecular_weights = (
                self.constants.MWs if self.constants is not None else []
            )
            if not molecular_weights or None in molecular_weights:
                return (
                    None,
                    None,
                    "Could not retrieve molecular weight for all components.",
                )
            molecular_weight = (
                sum(
                    fraction * component_mw
                    for fraction, component_mw in zip(zs, molecular_weights)
                )
                / 1000.0
            )
            if molecular_weight <= 0:
                return (
                    None,
                    None,
                    "Calculated average molecular weight is zero or negative.",
                )

            density_liquid = None
            density_gas = None

            if len(component_names) == 1 and isinstance(
                flasher,
                FlashPureVLS,
            ):
                pure_phase = getattr(result, "phase", "").lower()
                if pure_phase in ["l", "liquid"]:
                    if (
                        result.liquid0
                        and result.liquid0.V() is not None
                        and result.liquid0.V() > 0
                    ):
                        density_liquid = molecular_weight / result.liquid0.V()
                elif pure_phase in ["g", "v", "gas", "vapor"]:
                    if (
                        result.gas
                        and result.gas.V() is not None
                        and result.gas.V() > 0
                    ):
                        density_gas = molecular_weight / result.gas.V()
                elif "/" in pure_phase:
                    if (
                        result.liquid0
                        and result.liquid0.V() is not None
                        and result.liquid0.V() > 0
                    ):
                        density_liquid = molecular_weight / result.liquid0.V()
                    if (
                        result.gas
                        and result.gas.V() is not None
                        and result.gas.V() > 0
                    ):
                        density_gas = molecular_weight / result.gas.V()
                elif hasattr(result, "phases"):
                    for phase in result.phases:
                        if (
                            phase.is_liquid
                            and phase.V() is not None
                            and phase.V() > 0
                        ):
                            density_liquid = molecular_weight / phase.V()
                        elif (
                            phase.is_gas
                            and phase.V() is not None
                            and phase.V() > 0
                        ):
                            density_gas = molecular_weight / phase.V()
                else:
                    return (
                        None,
                        None,
                        f"Unexpected pure phase result: {pure_phase}",
                    )
            elif hasattr(result, "phases"):
                for phase in result.phases:
                    phase_volume = phase.V()
                    if phase_volume is not None and phase_volume > 0:
                        if phase.is_liquid:
                            density_liquid = molecular_weight / phase_volume
                        elif phase.is_gas:
                            density_gas = molecular_weight / phase_volume
                    else:
                        print(
                            "Warning: Phase volume is None or zero for phase "
                            f"type (Liq={phase.is_liquid}, Gas={phase.is_gas})"
                        )
            else:
                return (
                    None,
                    None,
                    "Flash result object did not contain expected "
                    "'phases' attribute.",
                )

            if density_liquid is not None and density_gas is not None:
                print(
                    "Determined Phase: Two-Phase, "
                    f"DensLiq={density_liquid:.3f}, "
                    f"DensGas={density_gas:.3f}"
                )
                return (density_liquid, density_gas), "Two-Phase", None
            if density_liquid is not None:
                print(
                    f"Determined Phase: Liquid, DensLiq={density_liquid:.3f}"
                )
                return density_liquid, "Liquid", None
            if density_gas is not None:
                print(f"Determined Phase: Vapor, DensGas={density_gas:.3f}")
                return density_gas, "Vapor", None

            print("Warning: Could not calculate density for any phase found.")
            phase_guess = "Unknown"
            if hasattr(result, "VF"):
                if result.VF == 0:
                    phase_guess = "Liquid (calc failed)"
                elif result.VF == 1:
                    phase_guess = "Vapor (calc failed)"
                elif 0 < result.VF < 1:
                    phase_guess = "Two-Phase (calc failed)"
            return None, phase_guess, "Failed to calculate density values."
        except Exception as exc:
            import traceback

            print(
                "Error during density calculation: "
                f"{exc}\n{traceback.format_exc()}"
            )
            return None, None, str(exc)

    def calculate_bubble_point(
        self,
        pressure_pa: float,
    ) -> BubblePointCalculation:
        """Calculate the bubble-point temperature."""
        if not self.components:
            return None, "No components selected."
        total = sum(self.components.values())
        if total == 0:
            return None, "Total fraction is zero."

        zs = [value / total for value in self.components.values()]
        component_names = list(self.components.keys())
        print(
            f"Calculating bubble point for: {component_names} "
            f"at P={pressure_pa:.1f} Pa"
        )

        try:
            if self.constants is None or self.properties is None:
                self.constants, self.properties = (
                    ChemicalConstantsPackage.from_IDs(component_names)
                )

            if (
                None in self.constants.Tcs
                or None in self.constants.Pcs
                or None in self.constants.omegas
            ):
                missing_components = [
                    component_names[index]
                    for index, critical_temperature in enumerate(
                        self.constants.Tcs
                    )
                    if critical_temperature is None
                ]
                return (
                    None,
                    "Missing critical properties for: "
                    + ", ".join(missing_components),
                )

            eos_kwargs = {
                "Pcs": self.constants.Pcs,
                "Tcs": self.constants.Tcs,
                "omegas": self.constants.omegas,
            }
            heat_capacity_gases = (
                self.properties.HeatCapacityGases
                if hasattr(self.properties, "HeatCapacityGases")
                else None
            )

            if len(component_names) == 1:
                gas = CEOSGas(
                    PRMIX,
                    eos_kwargs=eos_kwargs,
                    HeatCapacityGases=heat_capacity_gases,
                )
                liquid = CEOSLiquid(
                    PRMIX,
                    eos_kwargs=eos_kwargs,
                    HeatCapacityGases=heat_capacity_gases,
                )
                flasher = FlashPureVLS(
                    self.constants,
                    self.properties,
                    gas=gas,
                    liquids=[liquid],
                    solids=[],
                )
            else:
                gas = CEOSGas(
                    PRMIX,
                    HeatCapacityGases=heat_capacity_gases,
                    eos_kwargs=eos_kwargs,
                    zs=zs,
                )
                liquid = CEOSLiquid(
                    PRMIX,
                    HeatCapacityGases=heat_capacity_gases,
                    eos_kwargs=eos_kwargs,
                    zs=zs,
                )
                flasher = FlashVL(
                    self.constants,
                    self.properties,
                    liquid=liquid,
                    gas=gas,
                )

            result = flasher.flash(P=pressure_pa, VF=0, zs=zs)
            print(f"Bubble point flash result T={result.T}")
            return result.T - 273.15, None
        except Exception as exc:
            import traceback

            print(
                "Error during bubble point calculation: "
                f"{exc}\n{traceback.format_exc()}"
            )
            return None, str(exc)

    def calculate_lhv(
        self,
        lhv_data: Mapping[str, float],
    ) -> tuple[float, list[str]]:
        """Retain the legacy convenience API while LHV remains a domain rule."""
        return calculate_mixture_lhv(self.components, lhv_data)
