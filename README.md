# Thermo Components

Thermo Components is a PyQt6 desktop application for calculating thermodynamic properties of gas and liquid mixtures. It currently provides a GUI for component selection, composition entry, density and phase calculations, bubble-point estimation, lower heating value (LHV) lookup, flow conversion, and Excel report export.

## Status

The application works today as a single-user desktop tool. Pure rules live in `src/thermo_components/domain`, workflow orchestration and ports live in `src/thermo_components/application`, and external thermodynamics, SQLite, Excel, packaging, and several Qt UI concerns live under `src/thermo_components/adapters`. The main PyQt window and startup composition are still primarily hosted by [density.py](density.py).

## Features

- Multi-component mixture entry in `Mol %` or `Wt %`
- Density at selected, normal, and standard conditions
- Phase and bubble-point calculation
- Mixture LHV reporting on volumetric and mass bases
- Flow conversion between mass and reference-volume units
- Excel report export with warnings and calculated results
- Special handling for pure-water calculations through `IAPWS-95`

## Quick Start

1. Create and activate a Python virtual environment.
2. Install dependencies:

```powershell
python -m pip install -r requirements.txt
```

3. Run the desktop app:

```powershell
python density.py
```

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
- Contribution workflow: [CONTRIBUTING.md](CONTRIBUTING.md)
- Change history: [DEVLOG.md](DEVLOG.md)

## Current Repository Layout

```text
thermo_components/
  density.py          # Main application entry point and current monolith
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
      application/     # Typed DTOs and workflow use cases
      adapters/        # Thermo, SQLite, Excel, and packaging integrations
      bootstrap/       # Reserved for dependency wiring
  tests/               # Characterization tests
  docs/               # Architecture and roadmap documents
```

## Architectural Direction

The target design is a practical hexagonal architecture:

- The domain layer will own mixture rules, units, warnings, and result semantics.
- The application layer will own use cases such as property calculation, flow conversion, normalization, and report export.
- Adapters will isolate PyQt, `thermo`, SQLite, Excel, and packaging concerns.
- The current UI should remain functional throughout the migration; this is an incremental refactor, not a rewrite branch.

Phases 1, 2, and 3 are complete, and Phase 4 is in progress. The domain package owns pure rules, application use cases coordinate workflows through formal ports, and thermodynamics, SQLite LHV persistence, Excel reporting, and resource lookup are isolated in adapters. The Qt worker bridge, result-list presenter, calculation input collector, warning-banner controller, and composition-table setup/basis/total controller have been extracted; the remaining Phase 4 work is to keep splitting `MainWindow` into clearer controller/presenter boundaries.

## Development Notes

- `gui.py` is generated from `gui.ui` and should not become a home for business logic.
- Runtime behavior should be preserved during refactor work unless a change is intentional and documented.
- New behavior should be covered by tests as the architecture is extracted.

## License

GPL-3.0
