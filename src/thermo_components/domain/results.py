"""Calculation-result interpretation helpers."""


def extract_scalar_density_value(density_result, phase, error) -> float | None:
    """Return a density only when the flash produced one scalar value."""
    if error or density_result is None:
        return None
    if phase == "Two-Phase" or isinstance(density_result, tuple):
        return None
    try:
        return float(density_result)
    except (TypeError, ValueError):
        return None


def build_density_note(phase, error) -> str:
    """Build the report note associated with a density result."""
    if error:
        return error
    if phase == "Two-Phase":
        return "Two-Phase result"
    if phase:
        return f"Phase: {phase}"
    return ""
