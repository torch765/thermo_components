# Thermo Components

Thermo Components calculates thermodynamic properties of gas and liquid mixtures. It provides a PyQt6 desktop application and an in-progress FastAPI web application that reuse the same domain and application core.

## Status

The desktop tool remains fully operational. The web adapter provides a server-rendered calculator, health endpoint, and JSON calculation API. Pure rules live in `src/thermo_components/domain`, workflow orchestration and ports live in `src/thermo_components/application`, and framework integrations live under `src/thermo_components/adapters`.

## Features

- Multi-component mixture entry in `Mol %` or `Wt %`
- Simultaneous Mol %/Wt % display with active-basis normalization
- Density at selected, normal, and standard conditions
- Phase and bubble-point calculation
- Mixture LHV reporting on volumetric and mass bases
- Flow conversion between mass and reference-volume units
- Excel report export with warnings and calculated results
- Special handling for pure-water calculations through `IAPWS-95`
- JSON calculation API with typed validation and no PyQt dependency
- Complete volumetric and mass-basis LHV unit table
- Stateless web flow-conversion workspace with all desktop units

## Quick Start

1. Create and activate a Python virtual environment.
2. Install dependencies:

```powershell
python -m pip install -r requirements.txt
```

`requirements.txt` includes the editable local package install. For an older existing virtual environment, rerun the command above if `python density.py` cannot import `thermo_components`.

3. Run the desktop app:

```powershell
python density.py
```

4. Run the web application:

```powershell
python -m uvicorn thermo_components.adapters.web.app:app --reload
```

Then open:

- `http://127.0.0.1:8000/` for the browser calculator
- `http://127.0.0.1:8000/flow` for flow conversion
- `http://127.0.0.1:8000/health` for the health check
- `http://127.0.0.1:8000/docs` to try `POST /api/calculations`

## Tests

Install the development dependencies and run the characterization suite:

```powershell
python -m pip install -r requirements-dev.txt
python -m pytest
```

The tests run Qt in offscreen mode and do not display the application window.

## Project Docs

- Architecture target: [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md)
- Refactor roadmap: [docs/REFACTOR_ROADMAP.md](docs/REFACTOR_ROADMAP.md)
- Web app design basis: [docs/WEB_APP_DESIGN.md](docs/WEB_APP_DESIGN.md)
- Web stack ADR: [docs/adr/0001-web-stack.md](docs/adr/0001-web-stack.md)
- Contribution workflow: [CONTRIBUTING.md](CONTRIBUTING.md)
- Change history: [DEVLOG.md](DEVLOG.md)

## Current Repository Layout

```text
thermo_components/
  density.py          # Compatibility launcher and legacy import aliases
  gui.ui              # Qt Designer source
  gui.py              # Generated PyQt UI module
  lhv_data.py         # SQLite seed helper for LHV data
  lhv_data.db         # Runtime LHV lookup database
  density.spec        # PyInstaller spec
  pyproject.toml       # Test runner configuration
  requirements-dev.txt
  src/
    thermo_components/
      domain/          # Extracted framework-free business rules
      application/     # Typed DTOs, services, and workflow use cases
      adapters/        # Thermo, SQLite, Excel, packaging, Qt, and web integrations
      bootstrap/       # Desktop and web dependency composition
  tests/               # Characterization tests
  docs/               # Architecture and roadmap documents
```

## Architectural Direction

The target design is a practical hexagonal architecture:

- The domain layer will own mixture rules, units, warnings, and result semantics.
- The application layer will own use cases such as property calculation, flow conversion, normalization, and report export.
- Adapters will isolate PyQt, `thermo`, SQLite, Excel, and packaging concerns.
- The current UI should remain functional throughout the migration; this is an incremental refactor, not a rewrite branch.

Desktop Phases 1 through 6 and Web Phases 0 through 4D are complete. Web Phase 5 report download is next. The HTML form and JSON API share one calculation handler and use request-scoped thermo gateways to avoid shared calculation state.

## Development Notes

- `gui.py` is generated from `gui.ui` and should not become a home for business logic.
- Runtime behavior should be preserved during refactor work unless a change is intentional and documented.
- New behavior should be covered by tests as the architecture is extracted.

## License

GPL-3.0
