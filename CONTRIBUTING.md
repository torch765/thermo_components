# Contributing

## Purpose

This repository is being converted from a working desktop-script application into a maintainable software project. The immediate goal is controlled refactoring, not rapid feature churn.

## Development Setup

1. Create and activate a Python virtual environment.
2. Install development dependencies:

```powershell
python -m pip install -r requirements-dev.txt
```

3. Run the desktop application:

```powershell
python density.py
```

4. Run the web skeleton:

```powershell
python -m uvicorn thermo_components.adapters.web.app:app --reload
```

Use `http://127.0.0.1:8000/` for the calculator and `/docs` for the JSON API.

5. Run the test suite:

```powershell
python -m pytest
```

## Working Rules

- Preserve runtime behavior unless the change is intentional and documented.
- Prefer incremental extraction over large rewrites.
- Keep framework code out of the domain layer.
- Keep generated files generated. In particular, do not hand-edit `gui.py`.
- Add or update tests before moving critical logic when practical.

## Refactor Rules

- Use branch-by-abstraction: introduce new modules behind stable interfaces, then move callers gradually.
- When extracting logic from `density.py`, keep the root file as a launcher and compatibility import surface.
- Avoid leaking `PyQt6`, `thermo`, `sqlite3`, or `openpyxl` objects across architectural boundaries.
- Keep FastAPI and Pydantic request/response types inside the web adapter.
- Do not mix UI state manipulation with calculation orchestration in the same new module.

## Proposed Target Layout

```text
src/thermo_components/
  domain/
  application/
  adapters/
  bootstrap/
tests/
```

See [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) and [docs/REFACTOR_ROADMAP.md](docs/REFACTOR_ROADMAP.md) for the full target and migration plan.

## Testing Expectations

- Pure domain logic should have unit tests.
- The domain and application import-boundary tests must remain green.
- Adapter behavior should be covered with focused integration tests where needed.
- Qt tests should stay narrow and verify controller/presenter wiring rather than widget cosmetics.
- Web tests should verify routes and adapter boundaries without calling PyQt code.
- Server-rendered page changes should test successful rendering, validation, and static assets.
- High-risk areas include basis conversion, density reference handling, warning generation, flow conversion, and report row assembly.

## Documentation Expectations

- Update `README.md` when startup, layout, or developer workflow changes.
- Update `DEVLOG.md` for meaningful project changes.
- Update architecture and roadmap docs when the refactor plan materially changes.

## Pull Request Scope

- Keep each change narrow enough to review by responsibility.
- Prefer one architectural move per change set.
- If a refactor creates temporary duplication, record the cleanup step in the roadmap or changelog.
