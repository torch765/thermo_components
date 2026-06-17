# ADR 0001: Web Stack For Online Calculator

Status: Accepted

Date: 2026-06-17

## Context

Thermo Components has been refactored into a DDD-inspired, hexagonal structure. The domain and application layers are now reusable outside the PyQt desktop adapter. The next strategic goal is to make the tool available online.

The expected traffic is low to modest. The project should remain understandable for a small team and suitable for learning. The web version should be an online calculator, not a download page for the desktop executable.

## Decision

Use:

- FastAPI for the backend web adapter.
- Server-rendered HTML with Jinja2 for the first web UI.
- Uvicorn as the runtime server.
- DigitalOcean App Platform as the first deployment target.
- Bundled read-only `lhv_data.db` for MVP data.
- Temporary report generation for downloads.

Do not introduce:

- React/Vue/Svelte for the MVP.
- User accounts for the MVP.
- A persistent web database for the MVP.
- A web dependency on PyQt, `density.py`, `gui.py`, or `adapters/ui`.

## Rationale

FastAPI fits the existing Python codebase and can expose both HTML pages and JSON APIs. Server-rendered HTML keeps the MVP simple and avoids frontend build complexity. DigitalOcean App Platform is a managed deployment target that supports Python buildpacks and Git-based deployment.

The current architecture already separates domain/application code from Qt. A new web adapter can call the same application use cases without porting PyQt code.

## Consequences

Positive:

- Fast path to a real online tool.
- Low operational burden.
- Keeps business logic in reusable Python layers.
- Leaves room for JSON APIs or a richer frontend later.

Negative:

- Server-rendered HTML will eventually be less interactive than a single-page app.
- DigitalOcean App Platform local disk is ephemeral, so saved user data and reports require later storage decisions.
- A shared application facade is needed to avoid duplicating workflow assembly in web routes.

## Alternatives Considered

### Streamlit

Rejected for the main path. It is useful for demos, but less suitable for a polished public calculator with custom report downloads and longer-term web structure.

### Flask

Viable, but FastAPI gives stronger request modeling and a clean path to JSON endpoints.

### Django

Rejected for MVP. It is powerful, but heavier than needed without user accounts, admin workflows, or persistent database models.

### React Frontend With API Backend

Deferred. It may be useful later, but it adds build tooling and frontend state complexity before the core online workflow is proven.

### DigitalOcean Droplet

Deferred. A Droplet can be cheaper and more flexible, but it requires server management. App Platform is a better first deployment target for this project.
