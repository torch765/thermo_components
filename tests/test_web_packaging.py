from pathlib import Path
import tomllib


PROJECT_ROOT = Path(__file__).resolve().parents[1]
WEB_ROOT = (
    PROJECT_ROOT
    / "src"
    / "thermo_components"
    / "adapters"
    / "web"
)


def test_web_templates_and_static_assets_are_package_data():
    configuration = tomllib.loads(
        (PROJECT_ROOT / "pyproject.toml").read_text(encoding="utf-8")
    )
    package_data = configuration["tool"]["setuptools"]["package-data"][
        "thermo_components.adapters.web"
    ]

    assert sorted(package_data) == [
        "static/*.css",
        "static/*.js",
        "templates/*.html",
    ]
    assert (WEB_ROOT / "templates" / "base.html").is_file()
    assert (WEB_ROOT / "templates" / "calculator.html").is_file()
    assert (WEB_ROOT / "templates" / "flow.html").is_file()
    assert (WEB_ROOT / "static" / "calculator.css").is_file()
    assert (WEB_ROOT / "static" / "calculator.js").is_file()
    assert (WEB_ROOT / "static" / "flow.js").is_file()
