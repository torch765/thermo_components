## DEVLOG

### 2026-06-14
- Added project-level architecture documentation for the planned DDD-inspired, hexagonal refactor.
- Added a phased refactor roadmap with concrete PR sequencing and extraction order.
- Added a `CONTRIBUTING.md` guide to capture development rules, refactor constraints, and documentation expectations.
- Reworked `README.md` so the repo has a proper project overview, startup guide, and links to architecture documents.
- Completed refactor Phase 0 by adding the `src/thermo_components` package skeleton and pytest configuration.
- Added characterization tests for flow conversion, composition rules, route and warning policies, LHV derivation, normalization, and report projection.
- Added `requirements-dev.txt` to keep test tooling separate from runtime dependencies.
- Completed refactor Phase 1 by extracting framework-free domain modules for composition, reference conditions, flow conversion, LHV, thermo routes, warnings, and density-result interpretation.
- Updated `density.py` to delegate pure rules to the domain package while retaining its existing public names and UI behavior.
- Expanded the test suite to cover direct domain APIs and enforce the domain dependency boundary.
- Completed refactor Phase 2 by adding typed application DTOs and use cases for property calculation, flow conversion, composition handling, and report projection.
- Reduced `CalculationWorker` to a Qt execution bridge and changed it to emit typed calculation responses.
- Updated `MainWindow` handlers to delegate workflow decisions to the application layer while preserving existing rendering and export behavior.
- Expanded architecture checks to keep both domain and application layers free from framework and infrastructure imports.

### 2026-03-20
- Added an automatic pure-water model override so effectively pure water now uses `IAPWS-95` instead of the default `PRMIX` route.
- Suppressed the warning banner for pure water because it no longer uses `PRMIX`, while preserving the existing warning behavior for water-containing mixtures.
- Updated the displayed/exported model field to reflect the actual calculation route and preserved the shared selected/normal/standard density state used by the Flow tab.
- Added a persistent, non-modal warning banner above the main Results area for water-containing `PRMIX` calculations.
- Added specific thermo warnings for water-containing mixtures and water-containing two-phase behavior, and carried the same warning text into the Excel report `Warnings` section.
- Expanded the first thermo tab from one density output to three: selected-condition density, normal density, and standard density.
- Added normal reference density at `0 °C` and `1 atm` for future `Nm3`-basis work.
- Added standard reference density at `60 °F` and `1 atm` for future `Sm3` / liquid-hydrocarbon standard-volume work.
- Updated the Excel report export to include all three density rows and stored the scalar density values in shared app state for later reuse.
- Implemented the `Flow` tab with mass and reference-volume unit conversions across `kg`, `t`, `lb`, `Klb`, `Nm3`, `Sm3`, `bbl`, `SCFH`, `MSCFD`, and `MMSCFD`.
- Corrected the Flow-tab mass unit list so it uses `kg/d` instead of the previously implemented `kt/d` kilotonne-per-day unit.
- Wired the Flow tab to the densities produced by the first tab, using normal density for `Nm3` units and standard density for `Sm3`, `SCF`-family, and `bbl` units.
- Used `60 °F` and `1 atm` as the standard reference basis for standard-volume and stock-tank-style conversions, and bridged cross-basis volume conversions through mass.

### 2026-03-19
- Replaced the old clipboard-based results export with a PetraPlan-style Excel report export built with `openpyxl`.
- Renamed the UI button text to `Print Results`, regenerated `gui.py` from `gui.ui`, and removed the old clipboard copy behavior.
- Reports now save automatically beside the PyInstaller exe when frozen, or beside `density.py` when running from source.
- Expanded the Results-pane LHV block to show converted volumetric values (`kcal/Nm³`, `GJ/Nm³`, `MMBtu/Nm³`) and mass-basis values (`MJ/kg`, `MJ/t`, `GJ/kg`, `GJ/t`, `kcal/kg`, `kcal/t`, `MMBtu/kg`, `MMBtu/t`).
- Kept `MJ/Nm³` as the sole stored/base LHV unit and derived the mass-basis values from mixture MW using the ideal normal molar volume at 0 °C and 1 atm.
- Added `MMkcal/Nm³`, `MMkcal/kg`, and `MMkcal/t` display lines, derived from the existing kilocalorie conversions without changing the DB schema or GUI layout.

### 2026-02-01
- Built a new PyInstaller executable from the updated `density.py` and saved it as `dist\\density ver3.exe`.
- Updated `density.spec` so future builds output the versioned name directly (`name='density ver3'`).

### 2026-01-27
- Added `Copy Results` button support (`CopyResultButton`) to copy the full Input Composition table to the clipboard.
- Copy button is enabled only when the active composition column totals 100% (same gating as `Go!`).
- Shows a confirmation popup after copying.
- Clipboard format is tab-separated with a header row; includes the Total row.
- Updated `gui.py` to reflect the new widget added in `gui.ui`.
