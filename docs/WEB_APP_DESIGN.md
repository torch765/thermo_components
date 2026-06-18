# Web App Design Basis

Status: Accepted for MVP planning. Web Phase 4B is implemented and awaiting a browser smoke test.

## Purpose

Build an online version of Thermo Components as a real web tool, not a download page for the desktop executable.

The web app should reuse the existing domain and application core while replacing the PyQt adapter with a new web adapter. The initial goal is a small, maintainable MVP that can be deployed to DigitalOcean App Platform and improved incrementally.

## Design Principles

- Keep business rules out of the web framework.
- Reuse `domain`, `application`, and non-UI adapters wherever practical.
- Treat the web layer as another driving adapter, parallel to the Qt adapter.
- Prefer a simple server-rendered MVP before adding frontend complexity.
- Avoid user accounts, saved projects, and databases until there is a clear need.
- Preserve desktop app behavior while the web app is introduced.

## Recommended Stack

- Backend: FastAPI
- Frontend: server-rendered HTML with Jinja2 templates
- Styling: simple project-owned CSS, no large frontend framework for MVP
- Runtime server: Uvicorn
- Deployment: DigitalOcean App Platform
- Persistence for MVP: bundled read-only `lhv_data.db`
- Report files for MVP: generated per request as temporary downloads
- Authentication for MVP: none

DigitalOcean App Platform is appropriate for the first web release because it supports Python buildpacks, deploys from Git repositories or container images, and manages infrastructure. Its filesystem should be treated as ephemeral, so generated files must be temporary and persistent user data should wait for a managed database or object storage.

Reference docs:

- DigitalOcean App Platform overview: https://docs.digitalocean.com/products/app-platform/
- Python buildpack: https://docs.digitalocean.com/products/app-platform/reference/buildpacks/python/
- App Platform filesystem limits: https://docs.digitalocean.com/products/app-platform/details/limits/
- App Platform data storage guidance: https://docs.digitalocean.com/products/app-platform/how-to/store-data/

## MVP Scope

The first online version should support:

- A calculator form for components, basis, composition, temperature, pressure, and model.
- Simultaneous `Mol %` and `Wt %` display with one editable active basis.
- Composition normalization to exactly 100%.
- Server-side validation with clear error messages.
- A results page showing the same major calculated values as the desktop app.
- The complete existing LHV unit set on volumetric and mass bases.
- Flow conversion using the latest calculated normal and standard densities where applicable.
- Excel report download for the current calculation.
- Basic health check endpoint for deployment monitoring.

## Explicit Non-Goals For MVP

- User accounts or login.
- Saved calculations.
- Shared public calculation links.
- Payment, subscriptions, or quotas.
- Admin dashboard.
- Multi-user database.
- A React/Vue single-page app.
- Replacing the desktop app.
- Building or serving a Windows executable.

## Target User Flow

1. User opens the web calculator.
2. User selects components and enters composition.
3. User selects basis, temperature, pressure, and model.
4. User submits the form.
5. Server validates the request and runs the application use case.
6. Server renders results and warnings.
7. User optionally downloads an Excel report.
8. User can adjust inputs and recalculate.

## Proposed Architecture

```text
Browser
  |
  v
FastAPI routes
  |
  v
Web schemas and presenters
  |
  v
Application facade or use cases
  |
  v
Existing application/domain core
  |
  v
Thermo, LHV, and reporting adapters
```

The web layer should live under:

```text
src/thermo_components/adapters/web/
  app.py
  routes.py
  schemas.py
  presenters.py
  templates/
  static/
```

The web app must not import from:

```text
density.py
src/thermo_components/adapters/ui/
gui.py
```

The Qt adapter and web adapter should both depend inward on application use cases and shared adapters.

## Application Facade

Before building substantial routes, add a UI-independent facade for the calculation workflow. This prevents web routes from duplicating desktop workflow assembly.

Candidate module:

```text
src/thermo_components/application/use_cases/run_calculation.py
```

or:

```text
src/thermo_components/application/services/calculation_session.py
```

Responsibilities:

- Accept a typed calculation request suitable for any UI.
- Call `CalculatePropertiesUseCase`.
- Optionally call `PrepareReportUseCase`.
- Return a structured response containing calculation data, report projection, warnings, and flow-density state.

Non-responsibilities:

- HTML rendering.
- FastAPI request/response objects.
- Qt widgets.
- File downloads.

## Initial Web Routes

Recommended route set for the first implementation:

```text
GET  /health
POST /api/calculations
GET  /
GET  /calculator
POST /calculator
POST /calculator/report
```

`GET /health` returns a simple health response.

`POST /api/calculations` accepts JSON input, runs the shared application
workflow, and returns typed calculation results. It is the first implemented
calculation route and provides the backend contract for the HTML form.

`GET /calculator` renders the form.

`POST /calculator` validates form input, runs the calculation, and renders results.

`POST /calculator/report` regenerates the calculation from submitted form data and streams an Excel file.

Planned desktop-parity API routes:

```text
POST /api/compositions/derive
POST /api/compositions/normalize
POST /api/flow-conversions
GET  /flow
```

The composition endpoints call the existing application use cases. The browser
must not contain a second copy of molecular-weight conversion or normalization
rules.

The flow endpoint accepts the input value, source/destination units, and the
latest normal/standard density values. The server remains stateless.

## Calculation API Contract

`POST /api/calculations` accepts one percentage value per component on the
selected basis:

```json
{
  "components": [
    {"name": "methane", "percentage": 100.0}
  ],
  "basis": "Mol %",
  "temperature_c": 25.0,
  "pressure_atm": 1.0,
  "model": "PRMIX",
  "include_report_projection": false
}
```

The response contains:

- the full property calculation
- selected, normal, and standard density details
- flow-density state for later conversion workflows
- warnings
- an optional report projection

The web adapter derives the inactive composition basis through
`DeriveCompositionUseCase`. It creates a fresh thermo gateway for every
calculation request because the gateway holds mutable composition state.

## Desktop-Parity Interaction Design

### Dual-Basis Composition

The calculator should render three columns:

```text
Component | Mol % | Wt %
```

- The selected basis column is editable.
- The inactive column is read-only and visually distinct.
- Component or percentage edits trigger a debounced derive request.
- Switching basis makes the currently displayed derived values editable and
  derives the other column.
- Any composition change marks existing thermo results as stale until the user
  recalculates.

This is progressive enhancement. A user without JavaScript can still submit
the existing active-basis form and receive a calculation.

### Normalization

Normalize operates only on the active basis:

1. Send active percentages to `NormalizeCompositionUseCase`.
2. Replace the active values with the normalized response.
3. Recalculate the inactive basis through `DeriveCompositionUseCase`.
4. Show a concise success or zero-total validation message.
5. Invalidate stale calculation and flow-density state.

### LHV Presentation

Use `build_lhv_display_values` as the only conversion source. Present the base
`MJ/Nm³` value prominently and the remaining volumetric and mass-basis values
in a compact engineering table or expandable section.

Do not calculate LHV unit conversions in Jinja templates or JavaScript.

### Flow Workspace

The Flow workspace should be reachable from the primary navigation and use
`ConvertFlowUseCase` through a JSON endpoint.

- Populate source and destination units from `FLOW_UNIT_ORDER`.
- Default to `kg/h` and `t/d`.
- Recalculate after value or unit changes with a short debounce.
- Store latest normal and standard densities in browser `sessionStorage`.
- Send densities with each request rather than storing user state on the
  server.
- Show the density values being used and permit manual override.
- Keep mass-to-mass and same-reference conversions available without a prior
  property calculation.
- Show the existing domain error messages when normal or standard density is
  required.

## Report Download Strategy

For MVP, reports should be generated per request and returned directly to the browser. Do not store report files permanently on App Platform local disk.

Preferred implementation:

- Build `ReportExportRequest`.
- Write to a temporary file or in-memory buffer.
- Return `FileResponse` or `StreamingResponse`.
- Delete temporary files after response where practical.

If users later need saved reports, use object storage or a database-backed record.

## Data Strategy

For MVP:

- Keep `lhv_data.db` as a bundled read-only resource.
- Load it at app startup using existing bootstrap/resource code.
- Do not allow users to modify LHV data online.

Later:

- Move LHV data and saved calculations to a managed database if needed.
- Consider DigitalOcean Managed PostgreSQL only when persistence is actually required.

## Deployment Basis

Initial deployment target:

- DigitalOcean App Platform
- GitHub repository source
- Python buildpack
- Uvicorn start command

Expected deployment command:

```text
python -m uvicorn thermo_components.adapters.web.app:app --host 0.0.0.0 --port ${PORT:-8000}
```

The exact command may need adjustment for App Platform environment variable handling.

## Risks And Unknowns

- The `thermo` dependency stack may increase build time or require App Platform build tuning.
- Report generation must avoid relying on persistent local disk.
- Form-based component entry needs careful UX to avoid making the web MVP awkward.
- Desktop and web workflows may diverge unless shared application facades are introduced early.
- Long calculations may eventually require background jobs, but the MVP should start synchronously.

## Implementation Phases

### Web Phase 0: Design Basis

- Add this design note.
- Add an ADR for the stack choice.
- Update the roadmap.

### Web Phase 1: Shared Application Facade

Status: Complete as of 2026-06-17.

- Add a UI-independent calculation workflow facade.
- Add tests independent of Qt and FastAPI.
- Facade location: `src/thermo_components/application/services/calculation_session.py`.

### Web Phase 2: FastAPI Skeleton

Status: Complete as of 2026-06-17.

- Add FastAPI dependency.
- Add `adapters/web/app.py`.
- Add `/health`.
- Add minimal app startup tests.
- Local command: `python -m uvicorn thermo_components.adapters.web.app:app --reload`.

### Web Phase 3: Calculation API/Form Handler

Status: Complete as of 2026-06-18.

- Add web schemas and validation.
- Convert form/API input into application requests.
- Return calculation results through a web presenter.
- Add request-scoped thermo composition through `bootstrap/web.py`.
- Implement `POST /api/calculations`.

### Web Phase 4A: Core Server-Rendered UI

Status: Complete and browser-tested as of 2026-06-18.

- Add calculator template and styling.
- Render results and warnings.
- Keep design functional before polishing visuals.
- Reuse the calculation API handler rather than duplicating workflow mapping.
- Serve the calculator at `/` and `/calculator`.

### Web Phase 4B: Dual-Basis Composition And Normalization

Status: Implementation complete as of 2026-06-18; browser smoke test pending.

- Render both composition bases.
- Add derive and normalize endpoints.
- Add live updates and result-staleness behavior.
- Endpoints: `POST /api/compositions/derive` and
  `POST /api/compositions/normalize`.

### Web Phase 4C: Expanded LHV Results

- Add the complete volumetric and mass-basis LHV table.
- Reuse domain display-value rules.
- Test web values against desktop/domain expectations.

### Web Phase 4D: Flow Conversion Workspace

- Add `/flow` and `POST /api/flow-conversions`.
- Reuse `ConvertFlowUseCase` and all existing flow units.
- Carry latest calculation densities in browser session state.

### Web Phase 5: Report Download

- Add report download endpoint.
- Ensure no persistent local report storage.
- Add tests around generated response.

### Web Phase 6: DigitalOcean Deployment

- Add App Platform notes or app spec.
- Confirm build and runtime command.
- Deploy from GitHub.
- Smoke test online app.

## Acceptance Criteria For MVP

- User can complete a calculation online without desktop PyQt.
- Web code does not import Qt adapter modules.
- Existing desktop app still runs.
- Test suite covers shared facade and core web routes.
- App can deploy on DigitalOcean App Platform from the repository.
- Reports can be downloaded without persistent local file assumptions.
