# Thermo Components

Thermo Components is a PyQt6 desktop application for calculating thermodynamic properties of gas and liquid mixtures. It currently provides a GUI for component selection, composition entry, density and phase calculations, bubble-point estimation, lower heating value (LHV) lookup, flow conversion, and Excel report export.

## Status

The application works today as a single-user desktop tool. The UI, orchestration, and external integrations are still largely organized around [density.py](density.py), while the first architecture extraction has moved pure calculation rules into `src/thermo_components/domain`.

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
      application/     # Reserved for Phase 2 use cases
      adapters/        # Reserved for external integrations
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

Phase 1 is complete. Composition, reference conditions, flow conversion, LHV, route selection, warnings, and density-result interpretation now live in the domain package. Phase 2 will introduce application use cases and typed workflow DTOs.

## Development Notes

- `gui.py` is generated from `gui.ui` and should not become a home for business logic.
- Runtime behavior should be preserved during refactor work unless a change is intentional and documented.
- New behavior should be covered by tests as the architecture is extracted.

## License

GPL-3.0
