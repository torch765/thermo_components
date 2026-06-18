"""FastAPI application entry point."""

from pathlib import Path

from fastapi import FastAPI
from fastapi.staticfiles import StaticFiles

from thermo_components.bootstrap.web import (
    WebDependencies,
    build_web_dependencies,
)

from .pages import page_router
from .routes import router


WEB_ROOT = Path(__file__).resolve().parent


def create_app(
    dependencies: WebDependencies | None = None,
) -> FastAPI:
    """Create the web application without importing desktop UI code."""

    web_app = FastAPI(
        title="Thermo Components",
        version="1.0.0",
    )
    web_app.state.dependencies = dependencies or build_web_dependencies()
    web_app.mount(
        "/static",
        StaticFiles(directory=WEB_ROOT / "static"),
        name="static",
    )

    @web_app.get("/health", tags=["system"])
    def health_check() -> dict[str, str]:
        return {
            "service": "thermo-components",
            "status": "ok",
        }

    web_app.include_router(router)
    web_app.include_router(page_router)
    return web_app


app = create_app()
