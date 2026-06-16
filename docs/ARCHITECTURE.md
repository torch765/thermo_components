# Architecture

## Intent

This document defines the target architecture for Thermo Components as it moves from a single-file desktop application toward a DDD-inspired, hexagonal design.

The objective is not "DDD for its own sake." The objective is to isolate the engineering parts that change for different reasons:

- thermodynamic rules and result semantics
- application workflows
- UI concerns
- third-party library integration
- persistence and report generation

## Current State

Phases 1, 2, and 3 have extracted framework-free rules, application workflows, thermodynamics integration, LHV persistence, resource location, and Excel report export. Phase 4 has started by extracting the Qt worker bridge, calculation input collector, warning-banner controller, and results-list presenter. The following responsibilities still live in [density.py](../density.py):

- Qt controller logic and thread wiring
- most PyQt widget state and rendering
- report export button coordination and user messaging
- startup and dependency composition

The extracted domain modules currently own:

- composition identities, normalization, and basis conversion
- reference-condition values and unit conversion
- flow-unit definitions and conversion
- LHV mixture and display rules
- thermodynamic route and warning policies
- density-result interpretation

The application layer currently owns:

- typed calculation request and response DTOs
- selected, normal, and standard property calculation orchestration
- flow conversion workflow
- composition normalization and inactive-basis derivation workflows
- report result and warning projection

`CalculatePropertiesUseCase` now depends on the application-owned `ThermoPropertyGateway` port. `ThermoGateway` implements that contract in `adapters/thermo` and isolates all direct `thermo` and `chemicals` imports. `density.py` retains `MixtureCalculator` only as a compatibility alias.

Startup resolves `lhv_data.db` through the application-owned `ResourceLocator` port and loads it through `LhvRepository`. `RuntimeResourceLocator` isolates PyInstaller bundle detection, while `SqliteLhvRepository` owns all runtime and seed-script SQL access.

Report export depends on the application-owned `ReportExporter` port. `OpenPyxlReportExporter` owns workbook layout, formatting, and all direct `openpyxl` imports.

The remaining monolith still creates predictable problems:

- UI code still owns framework coordination and some presentation decisions
- export and persistence logic are coupled to screen state
- changing one area raises regression risk in unrelated areas

## Architectural Style

The target is a pragmatic hexagonal architecture with a lean domain model.

```text
Driving adapters: Qt / CLI / tests
                 |
                 v
        Application Use Cases
           |             ^
           v             |
Domain Model       Application Ports
                         ^
                         |
Driven adapters: thermo / SQLite / Excel / packaging
```

Dependency direction is always inward:

- adapters depend on application and domain
- application depends on domain
- domain depends on nothing framework-specific

## Bounded Areas

The codebase is small enough that full-blown bounded contexts would be excessive, but the domain still has clear responsibility areas:

### Composition

- component identities
- basis selection
- normalization
- validation
- mole/weight conversion

### Thermodynamic Properties

- route selection
- density semantics
- phase semantics
- bubble-point semantics
- warning policies

### Energy

- LHV lookup contract
- mixture LHV calculation
- display-unit derivation rules

### Flow Conversion

- unit definitions
- reference-basis handling
- density requirements for conversions

### Reporting

- report row definitions
- warning rows
- result projection for exports

## Target Package Layout

```text
src/thermo_components/
  domain/
    composition.py
    conditions.py
    flow_units.py
    flow_conversion.py
    lhv.py
    warnings.py
    thermo_routes.py
    results.py
  application/
    dto.py
    ports/
      persistence.py
      reporting.py
      resources.py
      thermo.py
    use_cases/
      calculate_properties.py
      convert_flow.py
      normalize_composition.py
      prepare_report.py
  adapters/
    ui/
      qt_main_window.py
      qt_worker.py
      presenters.py
    thermo/
      thermo_gateway.py
    persistence/
      sqlite_lhv_repository.py
    reporting/
      openpyxl_report_exporter.py
    packaging/
      resource_locator.py
  bootstrap/
    container.py
    main.py
```

## Core Domain Types

Prefer small typed objects over loose dicts and stringly-typed state.

Examples:

- `ComponentFraction`
- `CompositionInput`
- `NormalizedComposition`
- `OperatingConditions`
- `DensityValue`
- `ReferenceDensitySet`
- `ThermoRoute`
- `PropertyCalculationResult`
- `WarningMessage`
- `FlowConversionRequest`
- `FlowConversionResult`

Use `dataclass` and `Enum` first. Do not introduce richer tactical DDD patterns unless the invariants justify them.

## Application Ports

The application layer owns interfaces for the external capabilities required by its use cases. The domain remains independent of both ports and adapters.

### `ThermoPropertyGateway`

Responsibilities:

- calculate density for a given composition and conditions
- calculate bubble point or saturation temperature
- expose route-specific calculation capabilities

### `LhvRepository`

Responsibilities:

- return LHV values by normalized component identity

### `ReportExporter`

Responsibilities:

- export a typed report projection to a chosen output format

### `ResourceLocator`

Responsibilities:

- resolve runtime files in source and packaged modes

## Adapter Responsibilities

### Qt Adapter

- read widget inputs
- dispatch application use cases
- render presenter output
- manage threads and user feedback

The Qt layer should not calculate business results itself.

### Thermo Adapter

- translate domain/application requests into `thermo` calls
- isolate `PRMIX`, `IAPWS95`, flash objects, and library-specific error handling
- return domain-friendly result types

### Persistence Adapter

- load LHV data from SQLite
- hide SQL details from the rest of the system

### Reporting Adapter

- turn report DTOs into Excel output through `openpyxl`

### Packaging Adapter

- resolve files for source mode and PyInstaller mode

## Application Layer

Use cases should orchestrate workflows and compose domain rules with ports.

Primary use cases:

- `CalculatePropertiesUseCase`
- `ConvertFlowUseCase`
- `NormalizeCompositionUseCase`
- `ExportReportUseCase`

The application layer is the right place for:

- transaction-like orchestration
- error mapping
- response assembly
- coordinating multiple domain services and ports

It is not the right place for widget manipulation or direct third-party library calls.

## Mapping From Current Code

The current module should be decomposed as follows:

| Current responsibility | Current location | Target location |
| --- | --- | --- |
| Flow unit constants and conversion | Extracted | `domain/flow_units.py`, `domain/flow_conversion.py` |
| LHV display derivation | Extracted | `domain/lhv.py` |
| Composition and basis rules | Extracted | `domain/composition.py` |
| Reference conditions | Extracted | `domain/conditions.py` |
| Water route and warning policy | Extracted | `domain/thermo_routes.py`, `domain/warnings.py` |
| Density result interpretation | Extracted | `domain/results.py` |
| Property calculation orchestration | Extracted | `application/use_cases/calculate_properties.py` |
| Flow and composition workflows | Extracted | `application/use_cases/convert_flow.py`, `application/use_cases/normalize_composition.py` |
| Report projection | Extracted | `application/use_cases/prepare_report.py` |
| Qt thread bridge | Extracted; imported by launcher | `adapters/ui/qt_worker.py` |
| `thermo` integration | Extracted; compatibility alias remains | `application/ports/thermo.py`, `adapters/thermo/thermo_gateway.py` |
| Excel export | Extracted; button coordination remains | `application/ports/reporting.py`, `adapters/reporting/openpyxl_report_exporter.py` |
| LHV DB loading | Extracted; compatibility wrapper remains | `application/ports/persistence.py`, `adapters/persistence/sqlite_lhv_repository.py` |
| Result-list presentation | Extracted | `adapters/ui/presenters.py` |
| Calculation input collection | Extracted | `adapters/ui/input_collection.py` |
| Thermo warning banner | Extracted | `adapters/ui/warning_banner.py` |
| UI state and rendering | `MainWindow` | `adapters/ui/qt_main_window.py`, additional presenters |
| Resource lookup | Extracted; compatibility wrapper remains | `application/ports/resources.py`, `adapters/packaging/resource_locator.py` |
| Startup composition | `main` | `bootstrap/main.py` |

## Dependency Rules

- `domain` imports only the standard library unless there is a compelling reason otherwise.
- `application` may import `domain`, but not framework libraries.
- `adapters` may import `application` and `domain`.
- `bootstrap` wires concrete adapters into application services.
- Generated Qt UI files stay isolated from business code.

## Non-Goals

These are explicitly not part of the first architecture pass:

- microservices
- distributed messaging
- event sourcing
- plugin systems
- abstracting every helper into its own class

The project is still a local desktop app. The architecture should reflect that.
