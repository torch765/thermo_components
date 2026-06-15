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

Phases 1 and 2 have extracted framework-free rules and application workflows. The following responsibilities still live in [density.py](../density.py):

- `thermo` library orchestration
- Qt thread and signal setup
- PyQt widget state and rendering
- Excel export
- resource lookup and startup

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

`CalculatePropertiesUseCase` currently depends on a local calculator protocol. Phase 3 will move that contract into a formal port and place the existing `MixtureCalculator` implementation behind an adapter.

The remaining monolith still creates predictable problems:

- UI code owns business workflows
- calculation rules are hard to test without the GUI
- `thermo` implementation details leak everywhere
- export and persistence logic are coupled to screen state
- changing one area raises regression risk in unrelated areas

## Architectural Style

The target is a pragmatic hexagonal architecture with a lean domain model.

```text
UI / CLI / tests
      |
      v
Application Use Cases
      |
      v
Domain Model + Domain Services
      |
      v
Ports
      |
      v
Adapters: thermo, Qt, SQLite, Excel, packaging
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
    ports.py
  application/
    dto.py
    use_cases/
      calculate_properties.py
      convert_flow.py
      normalize_composition.py
      export_report.py
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

## Domain Ports

The domain and application layers should depend on interfaces, not libraries.

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
| Qt thread bridge | `CalculationWorker` | `adapters/ui/qt_worker.py` |
| `thermo` integration | `MixtureCalculator` | `adapters/thermo/thermo_gateway.py` |
| Excel export | `MainWindow.export_results_to_excel` | `adapters/reporting/openpyxl_report_exporter.py` |
| LHV DB loading | `load_lhv_data`, `lhv_data.py` | `adapters/persistence/sqlite_lhv_repository.py` |
| UI state and rendering | `MainWindow` | `adapters/ui/qt_main_window.py`, `adapters/ui/presenters.py` |
| Startup and resource lookup | `main`, `resource_path` | `bootstrap/main.py`, `adapters/packaging/resource_locator.py` |

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
