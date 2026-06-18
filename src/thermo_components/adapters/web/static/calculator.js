(() => {
  const densityState = document.querySelector("[data-flow-density-state]");
  if (densityState) {
    const numberOrNull = (value) => {
      const parsed = Number.parseFloat(value);
      return Number.isFinite(parsed) && parsed > 0 ? parsed : null;
    };
    const storedDensities = {
      selected_density_kg_m3: numberOrNull(
        densityState.dataset.selectedDensity,
      ),
      normal_density_kg_m3: numberOrNull(
        densityState.dataset.normalDensity,
      ),
      standard_density_kg_m3: numberOrNull(
        densityState.dataset.standardDensity,
      ),
    };
    try {
      window.sessionStorage.setItem(
        "thermo-components.flow-densities",
        JSON.stringify(storedDensities),
      );
    } catch {
      // Conversion still works with manual density entry.
    }
  }

  const form = document.querySelector("[data-calculator-form]");
  const list = document.querySelector("[data-composition-list]");
  const template = document.querySelector("#composition-row-template");
  const addButton = document.querySelector("[data-add-row]");
  const normalizeButton = document.querySelector("[data-normalize]");
  const totalValue = document.querySelector("[data-total-value]");
  const totalIndicator = document.querySelector("[data-total-indicator]");
  const activeBasisLabel = document.querySelector("[data-active-basis-label]");
  const feedback = document.querySelector("[data-composition-feedback]");
  const resultsPanel = document.querySelector(".results-panel");

  if (!form || !list || !template) {
    return;
  }

  let deriveTimer;
  let deriveSequence = 0;

  const getBasis = () => (
    form.querySelector('input[name="basis"]:checked')?.value || "Mol %"
  );

  const getRows = () => (
    [...list.querySelectorAll("[data-composition-row]")]
  );

  const getBasisInput = (row, basis) => (
    row.querySelector(`[data-basis-input="${basis}"]`)
  );

  const setFeedback = (message = "", state = "") => {
    if (!feedback) {
      return;
    }
    feedback.textContent = message;
    feedback.classList.toggle("success", state === "success");
    feedback.classList.toggle("error", state === "error");
  };

  const markResultsStale = () => {
    resultsPanel?.classList.add("is-stale");
  };

  const setBasisState = () => {
    const activeBasis = getBasis();
    const inactiveBasis = activeBasis === "Mol %" ? "Wt %" : "Mol %";

    if (activeBasisLabel) {
      activeBasisLabel.textContent = activeBasis;
    }

    getRows().forEach((row) => {
      const activeInput = getBasisInput(row, activeBasis);
      const inactiveInput = getBasisInput(row, inactiveBasis);
      const activeCell = activeInput?.closest("[data-basis-cell]");
      const inactiveCell = inactiveInput?.closest("[data-basis-cell]");

      if (activeInput) {
        activeInput.readOnly = false;
      }
      if (inactiveInput) {
        inactiveInput.readOnly = true;
      }
      activeCell?.classList.add("is-active");
      activeCell?.classList.remove("is-derived");
      inactiveCell?.classList.add("is-derived");
      inactiveCell?.classList.remove("is-active");
    });
  };

  const updateTotal = () => {
    const activeBasis = getBasis();
    const total = getRows().reduce((sum, row) => {
      const input = getBasisInput(row, activeBasis);
      return sum + (Number.parseFloat(input?.value) || 0);
    }, 0);

    if (totalValue) {
      totalValue.textContent = `${total.toFixed(2)}%`;
    }
    if (totalIndicator) {
      totalIndicator.classList.toggle(
        "invalid",
        Math.abs(total - 100) > 0.0001,
      );
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
      return "The composition request failed.";
    }
    return "The composition request failed.";
  };

  const deriveInactiveBasis = async ({ announce = false } = {}) => {
    const activeBasis = getBasis();
    const inactiveBasis = activeBasis === "Mol %" ? "Wt %" : "Mol %";
    const selectedRows = getRows().filter(
      (row) => row.querySelector('select[name="component_name"]')?.value,
    );

    getRows()
      .filter((row) => !selectedRows.includes(row))
      .forEach((row) => {
        const inactiveInput = getBasisInput(row, inactiveBasis);
        if (inactiveInput) {
          inactiveInput.value = "";
        }
      });

    if (!selectedRows.length) {
      return;
    }

    const sequence = ++deriveSequence;
    const response = await fetch("/api/compositions/derive", {
      method: "POST",
      headers: {"Content-Type": "application/json"},
      body: JSON.stringify({
        basis: activeBasis,
        components: selectedRows.map((row) => ({
          name: row.querySelector('select[name="component_name"]').value,
          percentage: (
            Number.parseFloat(getBasisInput(row, activeBasis)?.value) || 0
          ),
        })),
      }),
    });

    if (sequence !== deriveSequence) {
      return;
    }
    if (!response.ok) {
      setFeedback(await responseError(response), "error");
      return;
    }

    const body = await response.json();
    selectedRows.forEach((row, index) => {
      const inactiveInput = getBasisInput(row, inactiveBasis);
      const value = body.percentages[index];
      if (inactiveInput) {
        inactiveInput.value = value == null ? "" : value.toFixed(4);
      }
    });
    if (announce) {
      setFeedback(`${inactiveBasis} values updated.`, "success");
    } else if (feedback?.classList.contains("error")) {
      setFeedback();
    }
  };

  const scheduleDerive = () => {
    window.clearTimeout(deriveTimer);
    deriveTimer = window.setTimeout(() => {
      deriveInactiveBasis().catch(() => {
        setFeedback("Could not update the derived basis.", "error");
      });
    }, 250);
  };

  const bindRow = (row) => {
    row.querySelector('select[name="component_name"]')
      ?.addEventListener("change", () => {
        markResultsStale();
        scheduleDerive();
      });

    row.querySelectorAll("[data-percentage-input]").forEach((input) => {
      input.addEventListener("input", () => {
        if (input.dataset.basisInput !== getBasis()) {
          return;
        }
        markResultsStale();
        updateTotal();
        scheduleDerive();
      });
    });

    row.querySelector("[data-remove-row]")?.addEventListener("click", () => {
      const rows = getRows();
      if (rows.length === 1) {
        row.querySelector("select").value = "";
        row.querySelectorAll("input").forEach((input) => {
          input.value = "";
        });
      } else {
        row.remove();
      }
      markResultsStale();
      updateTotal();
      scheduleDerive();
    });
  };

  form.querySelectorAll('input[name="basis"]').forEach((radio) => {
    radio.addEventListener("change", () => {
      setBasisState();
      updateTotal();
      markResultsStale();
      deriveInactiveBasis({announce: true}).catch(() => {
        setFeedback("Could not switch the active basis.", "error");
      });
    });
  });

  addButton?.addEventListener("click", () => {
    const fragment = template.content.cloneNode(true);
    const row = fragment.querySelector("[data-composition-row]");
    list.appendChild(fragment);
    bindRow(row);
    setBasisState();
    row.querySelector("select")?.focus();
    updateTotal();
    markResultsStale();
  });

  normalizeButton?.addEventListener("click", async () => {
    const activeBasis = getBasis();
    const selectedRows = getRows().filter(
      (row) => row.querySelector('select[name="component_name"]')?.value,
    );
    if (!selectedRows.length) {
      setFeedback("Select at least one component.", "error");
      return;
    }

    const response = await fetch("/api/compositions/normalize", {
      method: "POST",
      headers: {"Content-Type": "application/json"},
      body: JSON.stringify({
        percentages: selectedRows.map(
          (row) => Number.parseFloat(
            getBasisInput(row, activeBasis)?.value,
          ) || 0,
        ),
        precision: 4,
      }),
    });
    if (!response.ok) {
      setFeedback(await responseError(response), "error");
      return;
    }

    const body = await response.json();
    selectedRows.forEach((row, index) => {
      getBasisInput(row, activeBasis).value = body.percentages[index];
    });
    updateTotal();
    markResultsStale();
    await deriveInactiveBasis();
    setFeedback(
      `${activeBasis} normalized to 100%; derived values refreshed.`,
      "success",
    );
  });

  getRows().forEach(bindRow);
  setBasisState();
  updateTotal();
  deriveInactiveBasis().catch(() => {
    setFeedback("Could not initialize the derived basis.", "error");
  });
})();
