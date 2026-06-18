"""FastAPI application entry point."""

from fastapi import FastAPI

from thermo_components.bootstrap.web import (
    WebDependencies,
    build_web_dependencies,
)

from .routes import router


def create_app(
    dependencies: WebDependencies | None = None,
) -> FastAPI:
    """Create the web application without importing desktop UI code."""

    web_app = FastAPI(
        title="Thermo Components",
        version="1.0.0",
    )
    web_app.state.dependencies = dependencies or build_web_dependencies()

    @web_app.get("/health", tags=["system"])
    def health_check() -> dict[str, str]:
        return {
            "service": "thermo-components",
            "status": "ok",
        }

    web_app.include_router(router)
    return web_app


app = create_app()
