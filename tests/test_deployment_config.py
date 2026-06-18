from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def test_digitalocean_app_spec_has_required_web_service_settings():
    app_spec = (PROJECT_ROOT / ".do" / "app.yaml").read_text(
        encoding="utf-8"
    )

    assert "repo: torch765/thermo_components" in app_spec
    assert "branch: main" in app_spec
    assert "environment_slug: python" in app_spec
    assert "--host 0.0.0.0 --port $PORT" in app_spec
    assert "http_path: /health" in app_spec
    assert "instance_size_slug: apps-s-1vcpu-0.5gb" in app_spec


def test_digitalocean_python_runtime_is_supported_and_pinned():
    runtime = (PROJECT_ROOT / "runtime.txt").read_text(
        encoding="utf-8"
    )

    assert runtime.strip() == "python-3.13.12"
