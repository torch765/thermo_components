# Refactor Roadmap

## Goal

Refactor Thermo Components from a working monolith into a maintainable DDD-inspired, hexagonal desktop application without breaking the current user workflow.

This roadmap is designed for incremental execution. The application must remain runnable after every phase.

## Ground Rules

- Preserve current behavior unless a change is intentional and documented.
- Prefer extraction over rewrite.
- Introduce seams before moving logic.
- Add tests around high-risk behavior before invasive moves.
- Keep PyQt and `thermo` out of new domain modules.

## Phased Plan

Current status: Phases 0, 1, 2, and 3 are complete. Phase 4 is in progress; the Qt worker bridge, result-list presenter, and calculation input collector have been extracted.

### Phase 0: Stabilize and Characterize

Status: Complete as of 2026-06-14.

Deliverables:

- Create a `tests/` layout.
- Add characterization tests for:
  - flow conversion
  - basis conversion and normalization
  - route selection and warning generation
  - LHV value derivation
  - report row assembly
- Document the target architecture and migration strategy.

Exit criteria:

- Critical rules are reproducible without launching the full GUI.
- Refactor decisions are written down and discoverable in-repo.

### Phase 1: Extract Pure Domain Modules

Status: Complete as of 2026-06-14.

Deliverables:

- Move constants, unit definitions, basis conversion, warning policies, and flow conversion into `src/thermo_components/domain/`.
- Replace loose helper dicts with typed dataclasses and enums where appropriate.
- Keep `density.py` calling the extracted functions.

Exit criteria:

- Pure business rules run without Qt, SQLite, Excel, or `thermo`.
- Existing GUI behavior is unchanged.

### Phase 2: Introduce Application Use Cases

Status: Complete as of 2026-06-14.

Deliverables:

- Add `application/use_cases/` for:
  - property calculation
  - flow conversion
  - normalization
  - report export preparation
- Convert worker result payloads from ad hoc dicts into response DTOs.

Exit criteria:

- Workflow orchestration no longer lives directly in `MainWindow`.
- The UI can call use cases instead of raw helper chains.

### Phase 3: Isolate External Dependencies Behind Ports

Status: Complete as of 2026-06-16.

Deliverables:

- Add ports for thermo properties, LHV lookup, report export, and resource location. Complete.
- Move `MixtureCalculator` logic behind a thermo adapter implementing the thermo gateway port. Complete.
- Move SQLite access behind an LHV repository adapter. Complete.
- Move Excel report generation behind an OpenPyXL report exporter adapter. Complete.

Exit criteria:

- The application layer depends only on ports and domain objects.
- `thermo` and `sqlite3` imports are confined to adapter modules.

### Phase 4: Refactor the Qt Layer

Status: In progress. Phase 4A started on 2026-06-16.

Deliverables:

- Split the current `MainWindow` into:
  - a Qt controller
  - a presenter
  - a worker/async bridge
- Remove business rules from widget event handlers.

Completed so far:

- Moved `CalculationWorker` into `adapters/ui/qt_worker.py`.
- Moved result-list formatting into `adapters/ui/presenters.py`.
- Moved calculation input collection into `adapters/ui/input_collection.py`.

Remaining:

- Introduce smaller Qt controller/presenter modules around table state, warnings, and report actions.

Exit criteria:

- Widget code mostly maps inputs and renders outputs.
- New UI behavior can be changed without touching calculation rules.

### Phase 5: Extract Reporting and Bootstrap

Deliverables:

- Move Excel export into a reporting adapter.
- Introduce a bootstrap module that wires adapters and use cases.
- Reduce `density.py` to a thin compatibility launcher or replace it with a package entry point.

Exit criteria:

- Startup and dependency wiring are explicit.
- Export logic no longer scrapes UI state directly.

### Phase 6: Consolidate and Clean Up

Deliverables:

- Remove dead helper code from the monolith.
- Update packaging paths and developer workflow docs.
- Expand test coverage around extracted adapters.

Exit criteria:

- The old monolithic structure is retired.
- The codebase reflects the documented architecture.

## Module-by-Module Move Plan

Move code out of `density.py` in this order:

1. Pure constants and conversion helpers
2. Composition normalization and basis conversion
3. Water-route selection and warning policies
4. Flow conversion rules
5. LHV derivation and display models
6. Report row builders
7. Calculation orchestration
8. `thermo` adapter
9. Qt controller/presenter split
10. Startup and packaging concerns

This order minimizes breakage because each step moves logic with progressively stronger framework coupling.

## First Three PRs

### PR 1: Project Skeleton and Characterization Tests

Scope:

- Add `src/thermo_components/` package skeleton.
- Add `tests/` package.
- Add characterization tests for existing pure behavior.
- Leave runtime imports in `density.py` for now.

Files likely introduced:

- `src/thermo_components/__init__.py`
- `tests/test_flow_conversion.py`
- `tests/test_composition_rules.py`
- `tests/test_lhv.py`
- `tests/test_warning_policies.py`

Acceptance criteria:

- Test suite covers critical helper behavior.
- No visible GUI behavior changes.

### PR 2: Extract Domain Rules

Scope:

- Move pure helper logic into domain modules.
- Add enums/dataclasses for basis, route, conditions, and result fragments.
- Update `density.py` imports to use domain modules.

Files likely introduced:

- `src/thermo_components/domain/composition.py`
- `src/thermo_components/domain/flow_units.py`
- `src/thermo_components/domain/flow_conversion.py`
- `src/thermo_components/domain/lhv.py`
- `src/thermo_components/domain/thermo_routes.py`
- `src/thermo_components/domain/warnings.py`

Acceptance criteria:

- `density.py` shrinks materially at the top-level helper section.
- All extracted modules remain free of framework imports.

### PR 3: Add Application Use Cases and DTOs

Scope:

- Introduce request/response DTOs.
- Move workflow orchestration out of `CalculationWorker`.
- Have the Qt worker call application services rather than domain helpers directly.

Files likely introduced:

- `src/thermo_components/application/dto.py`
- `src/thermo_components/application/use_cases/calculate_properties.py`
- `src/thermo_components/application/use_cases/convert_flow.py`
- `src/thermo_components/application/use_cases/normalize_composition.py`

Acceptance criteria:

- The worker becomes a thin adapter.
- UI and orchestration responsibilities are visibly separated.

## Later PR Sequence

After the first three PRs, the next sequence should be:

1. Thermo gateway port plus adapter
2. SQLite LHV repository adapter
3. Report exporter adapter
4. Qt presenter split
5. Bootstrap wiring and entry-point cleanup
6. Packaging cleanup and final monolith removal

## Risks and Controls

### Risk: Hidden behavior changes during extraction

Control:

- Characterization tests first
- Small PRs
- Keep old call sites temporarily if needed

### Risk: Over-engineering the domain

Control:

- Use dataclasses and small services first
- Introduce patterns only when they remove real coupling

### Risk: Qt refactor stalls because too much logic still sits in the window

Control:

- Extract application orchestration before tackling widget rendering

### Risk: `thermo` library details leak into the new core

Control:

- Define the gateway interface before moving the adapter
- Forbid direct `thermo` imports outside the adapter package

## Definition of Done

The refactor is complete when:

- domain code contains no Qt, `thermo`, SQLite, or Excel imports
- application workflows are expressed as explicit use cases
- adapters own all framework integration
- startup is a composition root, not a hidden side effect chain
- tests cover core rules and critical adapter behavior
- the current UI still works against the new core
