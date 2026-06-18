(() => {
  const workspace = document.querySelector("[data-flow-workspace]");
  if (!workspace) {
    return;
  }

  const storageKey = "thermo-components.flow-densities";
  const input = workspace.querySelector("[data-flow-input]");
  const fromUnit = workspace.querySelector("[data-flow-from]");
  const toUnit = workspace.querySelector("[data-flow-to]");
  const result = workspace.querySelector("[data-flow-result]");
  const resultUnit = workspace.querySelector("[data-flow-result-unit]");
  const error = workspace.querySelector("[data-flow-error]");
  const normalDensity = workspace.querySelector("[data-normal-density]");
  const standardDensity = workspace.querySelector("[data-standard-density]");
  const densityState = workspace.querySelector("[data-density-state]");
  const clearButton = workspace.querySelector("[data-clear-densities]");
  let conversionTimer;
  let conversionSequence = 0;

  const positiveNumberOrNull = (value) => {
    if (value == null || String(value).trim() === "") {
      return null;
    }
    const parsed = Number.parseFloat(value);
    return Number.isFinite(parsed) && parsed > 0 ? parsed : null;
  };

  const loadDensities = () => {
    try {
      const stored = JSON.parse(
        window.sessionStorage.getItem(storageKey) || "null",
      );
      if (!stored) {
        return false;
      }
      const normal = positiveNumberOrNull(stored.normal_density_kg_m3);
      const standard = positiveNumberOrNull(
        stored.standard_density_kg_m3,
      );
      normalDensity.value = normal ?? "";
      standardDensity.value = standard ?? "";
      return normal != null || standard != null;
    } catch {
      return false;
    }
  };

  const saveDensities = () => {
    try {
      window.sessionStorage.setItem(
        storageKey,
        JSON.stringify({
          normal_density_kg_m3: positiveNumberOrNull(normalDensity.value),
          standard_density_kg_m3: positiveNumberOrNull(
            standardDensity.value,
          ),
        }),
      );
    } catch {
      // Manual values remain available for the current page.
    }
  };

  const responseError = async (response) => {
    try {
      const body = await response.json();
      if (typeof body.detail === "string") {
        return body.detail;
      }
      if (Array.isArray(body.detail)) {
        return body.detail.map((item) => item.msg).join(" ");
      }
    } catch {
      return "The flow conversion request failed.";
    }
    return "The flow conversion request failed.";
  };

  const convert = async () => {
    const sequence = ++conversionSequence;
    resultUnit.textContent = toUnit.value;
    error.textContent = "";

    if (!input.value.trim()) {
      result.textContent = "—";
      return;
    }

    const response = await fetch("/api/flow-conversions", {
      method: "POST",
      headers: {"Content-Type": "application/json"},
      body: JSON.stringify({
        input_text: input.value,
        from_unit: fromUnit.value,
        to_unit: toUnit.value,
        normal_density_kg_m3: positiveNumberOrNull(normalDensity.value),
        standard_density_kg_m3: positiveNumberOrNull(
          standardDensity.value,
        ),
      }),
    });

    if (sequence !== conversionSequence) {
      return;
    }
    if (!response.ok) {
      result.textContent = "—";
      error.textContent = await responseError(response);
      return;
    }

    const body = await response.json();
    result.textContent = body.display_value || "—";
  };

  const scheduleConversion = () => {
    window.clearTimeout(conversionTimer);
    conversionTimer = window.setTimeout(() => {
      convert().catch(() => {
        result.textContent = "—";
        error.textContent = "Could not update the flow conversion.";
      });
    }, 180);
  };

  const hasStoredDensities = loadDensities();
  densityState.textContent = hasStoredDensities
    ? "Using reference densities from the latest calculator result."
    : "No calculator density is stored. Enter values when required.";

  input.addEventListener("input", scheduleConversion);
  fromUnit.addEventListener("change", scheduleConversion);
  toUnit.addEventListener("change", scheduleConversion);
  [normalDensity, standardDensity].forEach((densityInput) => {
    densityInput.addEventListener("input", () => {
      saveDensities();
      densityState.textContent = "Using manually adjusted density values.";
      scheduleConversion();
    });
  });

  clearButton.addEventListener("click", () => {
    normalDensity.value = "";
    standardDensity.value = "";
    try {
      window.sessionStorage.removeItem(storageKey);
    } catch {
      // Inputs are still cleared for the current page.
    }
    densityState.textContent = "Density values cleared.";
    scheduleConversion();
  });

  convert().catch(() => {
    result.textContent = "—";
    error.textContent = "Could not initialize the flow conversion.";
  });
})();
