(() => {
  const list = document.querySelector("[data-composition-list]");
  const template = document.querySelector("#composition-row-template");
  const addButton = document.querySelector("[data-add-row]");
  const totalValue = document.querySelector("[data-total-value]");
  const totalIndicator = document.querySelector("[data-total-indicator]");

  if (!list || !template) {
    return;
  }

  const updateTotal = () => {
    const total = [...list.querySelectorAll("[data-percentage-input]")]
      .reduce((sum, input) => sum + (Number.parseFloat(input.value) || 0), 0);

    if (totalValue) {
      totalValue.textContent = `${total.toFixed(2)}%`;
    }
    if (totalIndicator) {
      totalIndicator.classList.toggle("invalid", Math.abs(total - 100) > 0.0001);
    }
  };

  const bindRow = (row) => {
    row.querySelector("[data-percentage-input]")
      ?.addEventListener("input", updateTotal);

    row.querySelector("[data-remove-row]")?.addEventListener("click", () => {
      const rows = list.querySelectorAll("[data-composition-row]");
      if (rows.length === 1) {
        row.querySelector("select").value = "";
        row.querySelector("input").value = "";
      } else {
        row.remove();
      }
      updateTotal();
    });
  };

  list.querySelectorAll("[data-composition-row]").forEach(bindRow);

  addButton?.addEventListener("click", () => {
    const fragment = template.content.cloneNode(true);
    const row = fragment.querySelector("[data-composition-row]");
    list.appendChild(fragment);
    bindRow(row);
    row.querySelector("select")?.focus();
    updateTotal();
  });

  updateTotal();
})();
