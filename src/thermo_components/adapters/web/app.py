"""FastAPI application entry point."""

from fastapi import FastAPI


def create_app() -> FastAPI:
    """Create the web application without importing desktop UI code."""

    web_app = FastAPI(
        title="Thermo Components",
        version="1.0.0",
    )

    @web_app.get("/health", tags=["system"])
    def health_check() -> dict[str, str]:
        return {
            "service": "thermo-components",
            "status": "ok",
        }

    return web_app


app = create_app()
