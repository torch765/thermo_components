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

Current status: Phases 0 through 6 and Web Phases 0 through 4A are complete. Web desktop-parity increments 4B through 4D are planned before report download.

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

Status: Complete as of 2026-06-16.

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
- Moved thermo warning-banner rendering into `adapters/ui/warning_banner.py`.
- Moved composition-table setup, active-basis styling, and total validation into `adapters/ui/composition_table.py`.
- Moved composition row/list add, remove, duplicate, and clear widget mutations into `adapters/ui/composition_table.py`.
- Moved normalization active-column read/write and success feedback into `adapters/ui/composition_table.py`.
- Moved report export path selection, export dialogs, and optional file opening into `adapters/ui/report_controller.py`.
- Moved concrete desktop dependency composition into `bootstrap/desktop.py`.

Deferred to Phase 5:

- Remaining calculation/progress orchestration around `QThread`, result lifecycle, and progress animation.
- Flow-tab widget setup and conversion rendering.
- Report export request assembly from current conditions and composition rows.
- Final decision on moving `MainWindow` itself into `adapters/ui`.

Exit criteria:

- Widget code mostly maps inputs and renders outputs.
- New UI behavior can be changed without touching calculation rules.

### Phase 5: Thin Launcher and MainWindow Migration

Status: Complete as of 2026-06-17.

Deliverables:

- Move `MainWindow` into `adapters/ui/qt_main_window.py` or reduce `density.py` to importing it.
- Extract remaining calculation/progress workflow coordination only after adding targeted tests.
- Extract remaining flow-tab UI coordination if it still materially clutters `MainWindow`.
- Move report request assembly behind a small UI/application boundary if the current helper wrappers become a blocker.
- Keep legacy aliases and wrappers in `density.py` until tests and packaging no longer need them.

Completed so far:

- Added characterization tests for `MainWindow` dependency injection, flow-tab rendering, report request helpers, and calculation result/error lifecycle hooks.
- Moved Flow-tab widget setup, output configuration, signal wiring, and conversion rendering into `adapters/ui/flow_tab.py`.
- Moved report condition collection, composition row collection, projection calls, and export-request assembly into `adapters/ui/report_request.py`.
- Moved calculation input orchestration, worker/thread wiring, result/error rendering, finish-state handling, and progress reset into `adapters/ui/calculation_workflow.py`.
- Moved `MainWindow` into `adapters/ui/qt_main_window.py` and reduced `density.py` to a compatibility launcher plus legacy aliases.

Next:

- Start Phase 6 consolidation and cleanup.

Exit criteria:

- `density.py` is a thin compatibility launcher plus legacy import aliases.
- `MainWindow` construction uses the bootstrap dependency container.
- Remaining UI workflows have direct regression coverage before they move.

### Phase 6: Consolidate and Clean Up

Status: Complete as of 2026-06-17.

Deliverables:

- Remove dead helper code from the monolith.
- Update packaging paths and developer workflow docs.
- Expand test coverage around extracted adapters.

Completed so far:

- Removed stale imports left in the moved Qt main-window adapter.
- Replaced absolute local `.venv` paths in `density.spec` with PyInstaller hook-based package collection.
- Clarified developer setup docs around the editable package install.
- Moved tests to real adapter/bootstrap imports and kept only a focused `density.py` launcher compatibility test.
- Removed the test-only `QMessageBox` compatibility export from `density.py`.
- Added a packaging regression test that keeps `density.spec` portable and verifies bundled resource entries remain present.
- Decided to keep `MainWindow`, `MixtureCalculator`, `load_lhv_data`, and `resource_path` as launcher compatibility aliases for the next public release.
- Confirmed the Phase 6 branch with a final GUI smoke test.

Next:

- Optional web-readiness planning, focused on UI-independent application facades.

Exit criteria:

- The old monolithic structure is retired.
- The codebase reflects the documented architecture.

## Web App Plan

Design basis:

- See `docs/WEB_APP_DESIGN.md`.
- See `docs/adr/0001-web-stack.md`.

Target:

- Build an online calculator, not a download page for the desktop executable.
- Add a new FastAPI web adapter that reuses the existing domain/application core.
- Deploy the first MVP to DigitalOcean App Platform.
- Keep the PyQt desktop app working during web development.

### Web Phase 0: Design Basis

Status: Complete as of 2026-06-17.

Deliverables:

- Document MVP scope and non-goals.
- Record web stack decision.
- Define initial route, report download, data, and deployment assumptions.

### Web Phase 1: Shared Application Facade

Status: Complete as of 2026-06-17.

Deliverables:

- Add a UI-independent calculation workflow facade.
- Keep FastAPI and PyQt details out of the facade.
- Add tests independent of Qt and FastAPI.

Completed:

- Added `CalculationSessionService` in `src/thermo_components/application/services/calculation_session.py`.
- Composed property calculation and report projection behind a UI-independent request/response.
- Added application tests that do not import Qt, FastAPI, or `density.py`.

Exit criteria:

- Web and desktop adapters can share calculation workflow assembly.
- No web route needs to call PyQt code or `density.py`.

### Web Phase 2: FastAPI Skeleton

Status: Complete as of 2026-06-17.

Deliverables:

- Add FastAPI/Uvicorn dependencies.
- Add `src/thermo_components/adapters/web/app.py`.
- Add `/health`.
- Add basic app startup tests.

Completed:

- Added a FastAPI app factory and module-level Uvicorn entry point.
- Added `GET /health` for local and deployment smoke tests.
- Added web startup tests and import-boundary coverage for FastAPI/Uvicorn.

Exit criteria:

- The web app can run locally with Uvicorn.
- The desktop app still runs.

### Web Phase 3: Calculation Form/API

Status: Complete as of 2026-06-18.

Deliverables:

- Add web schemas and request validation.
- Translate form/API input into application requests.
- Render or return calculated results and warnings.

Completed:

- Added `POST /api/calculations` with typed Pydantic request and response schemas.
- Added active-basis composition input and inactive-basis derivation through the application layer.
- Added request-scoped thermo calculation sessions to prevent cross-request state leakage.
- Added validation, error-mapping, presenter, bootstrap, architecture, and real-calculation tests.

Exit criteria:

- A user can run a calculation without PyQt.

### Web Phase 4A: Core Server-Rendered UI

Status: Complete and browser-tested as of 2026-06-18.

Deliverables:

- Add Jinja2 templates and project CSS.
- Add calculator form, result view, and validation messages.

Completed:

- Added server-rendered calculator routes at `/` and `/calculator`.
- Added dynamic composition rows, basis selection, condition inputs, and inline validation.
- Added responsive result, warning, reference-density, and report-ready views.
- Kept the HTML form and JSON API on the same shared web calculation handler.
- Added focused page, form, static-asset, and real-calculation tests.

Exit criteria:

- The MVP is usable from a browser with no JavaScript framework.

### Web Phase 4B: Dual-Basis Composition And Normalization

Deliverables:

- Show `Mol %` and `Wt %` columns simultaneously.
- Keep the selected basis editable and the derived basis read-only.
- Derive the inactive basis through `DeriveCompositionUseCase`.
- Add a Normalize action backed by `NormalizeCompositionUseCase`.
- Update both columns and the active total without running a full thermo calculation.
- Mark displayed calculation results as stale when composition changes.

Implementation notes:

- Add typed composition derive/normalize web schemas and endpoints.
- Use small debounced browser requests; do not duplicate molecular-weight formulas in JavaScript.
- On basis change, activate the displayed derived column and derive the other column.
- Preserve a no-JavaScript fallback through the existing calculator form.

Exit criteria:

- Editing either active basis updates the other basis consistently.
- Normalize scales the active basis to exactly 100% and refreshes the derived basis.
- Zero-total and incomplete-row errors are clear and non-destructive.

### Web Phase 4C: Expanded LHV Results

Deliverables:

- Display all existing volumetric LHV units:
  - `MJ/Nm³`
  - `kcal/Nm³`
  - `MMkcal/Nm³`
  - `GJ/Nm³`
  - `MMBtu/Nm³`
- Display all existing mass-basis LHV units:
  - `MJ/kg`, `MJ/t`, `GJ/kg`, `GJ/t`
  - `kcal/kg`, `kcal/t`, `MMkcal/kg`, `MMkcal/t`
  - `MMBtu/kg`, `MMBtu/t`
- Reuse `build_lhv_display_values`; do not reproduce conversion constants in templates or JavaScript.
- Render missing-LHV warnings prominently.

Exit criteria:

- Web LHV values match desktop/domain values for representative mixtures.
- Units with unavailable data render as `N/A` rather than misleading zeroes.

### Web Phase 4D: Flow Conversion Workspace

Deliverables:

- Add a dedicated Flow page or primary workspace section.
- Expose all units in `FLOW_UNIT_ORDER`.
- Add `POST /api/flow-conversions` backed by `ConvertFlowUseCase`.
- Update the result immediately when value or units change.
- Reuse the latest calculated normal and standard densities.
- Allow density-independent conversions before a property calculation.
- Show clear density-required messages for conversions that need unavailable state.

Implementation notes:

- Keep the server stateless: send densities with each conversion request.
- Store the latest calculation densities in browser session state for navigation to the Flow page.
- Provide visible normal/standard density context and optional manual overrides.
- Add swap and copy controls only after the core conversion path is tested.

Exit criteria:

- All 17 desktop flow units are available online.
- Mass, same-reference-volume, mass/volume, and cross-reference conversions match domain tests.
- Flow conversion never imports or calls PyQt code.

### Web Phase 5: Report Download

Deliverables:

- Add report download endpoint.
- Generate reports per request without persistent local storage.

Exit criteria:

- A user can download an Excel report from the web calculation.

### Web Phase 6: DigitalOcean Deployment

Deliverables:

- Add App Platform deployment notes or app spec.
- Confirm build/start command.
- Deploy from GitHub.

Exit criteria:

- The MVP is reachable online and passes a deployment smoke test.

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
