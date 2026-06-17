from pathlib import Path

from PyQt6.QtWidgets import (
    QButtonGroup,
    QMainWindow,
    QMessageBox,
    QTableWidgetItem,
)
from PyQt6.QtCore import Qt, QTimer
from PyQt6.QtGui import QColor

# --- Import the generated UI class ---
from gui import Ui_Dialog

from thermo_components.domain.composition import (
    MOLECULAR_WEIGHTS,
    PURE_WATER_WARNING_FRACTION,
    WATER_COMPONENT_ALIASES,
    active_basis_amount_rows as _active_basis_amount_rows,
    derive_inactive_percentages,
    has_water_component,
    is_effectively_pure_water,
    is_water_component,
    normalize_component_identity,
    normalize_percentages,
    percentages_to_mole_fractions,
    water_fraction_active_basis,
)
from thermo_components.domain.conditions import (
    ATM_TO_PA,
    NORMAL_P_ATM,
    NORMAL_T_C,
    R,
    STANDARD_P_ATM,
    STANDARD_T_C,
    STANDARD_T_F,
    atm_to_pa,
    celsius_to_kelvin,
)
from thermo_components.domain.lhv import (
    KCAL_PER_MJ,
    MJ_PER_MMBTU,
    NORMAL_MOLAR_VOLUME_NM3_PER_KMOL,
)
from thermo_components.domain.results import (
    build_density_note,
    extract_scalar_density_value,
)
from thermo_components.adapters.ui.calculation_workflow import (
    QtCalculationWorkflowController,
)
from thermo_components.adapters.ui.composition_table import (
    ComponentAddStatus,
    ComponentRemoveStatus,
    CompositionTableController,
)
from thermo_components.adapters.ui.flow_tab import FlowTabController
from thermo_components.adapters.ui.report_request import QtReportRequestBuilder
from thermo_components.adapters.ui.warning_banner import (
    ThermoWarningBannerController,
)
from thermo_components.application.dto import (
    DeriveCompositionRequest,
    NormalizeCompositionRequest,
)
from thermo_components.bootstrap import (
    build_desktop_dependencies,
)


# --- MainWindow Class (Modified for gui.py) ---
class MainWindow(QMainWindow):
    def __init__(self, lhv_data=None, dependencies=None):
        super().__init__()

        # --- Set up the UI from Designer ---
        self.ui = Ui_Dialog()
        self.ui.setupUi(self) # Setup the UI onto this QMainWindow
        if hasattr(self.ui, "tabWidget"):
            self.ui.tabWidget.setCurrentIndex(0)
        self.composition_table = CompositionTableController(self.ui)

        # Step 2: Progress bar animation state
        self.progress_timer = None
        self.progress_target = 100
        self.progress_animation_duration = 800  # milliseconds
        self.progress_animation_start_time = 0
        self.progress_animation_start_value = 0

        dependencies = dependencies or build_desktop_dependencies(
            lhv_data={} if lhv_data is None else lhv_data,
            source_base_dir=Path.cwd(),
        )

        # --- Store LHV data and calculator instance ---
        self.dependencies = dependencies
        self.lhv_database = dependencies.lhv_database
        self.lhv_data_loaded = bool(self.lhv_database)
        self.calculator = dependencies.calculator
        self.calculate_properties_use_case = (
            dependencies.calculate_properties_use_case
        )
        self.convert_flow_use_case = dependencies.convert_flow_use_case
        self.derive_composition_use_case = (
            dependencies.derive_composition_use_case
        )
        self.normalize_composition_use_case = (
            dependencies.normalize_composition_use_case
        )
        self.prepare_report_use_case = dependencies.prepare_report_use_case
        self.report_exporter = dependencies.report_exporter
        self.report_export_controller = (
            dependencies.report_export_controller_factory(self)
        )
        self.last_result_data = None
        self.density_actual_kg_m3 = None
        self.density_normal_kg_m3 = None
        self.density_standard_kg_m3 = None
        self.flow_tab_controller = FlowTabController(
            self.ui,
            self.convert_flow_use_case,
            normal_density_provider=lambda: self.density_normal_kg_m3,
            standard_density_provider=lambda: self.density_standard_kg_m3,
        )
        self.report_request_builder = QtReportRequestBuilder(
            self.ui,
            self.prepare_report_use_case,
            lhv_data_available_provider=lambda: self.lhv_data_loaded,
            latest_result_provider=lambda: self.last_result_data,
        )
        self.calculation_workflow_controller = QtCalculationWorkflowController(
            self.ui,
            self.calculate_properties_use_case,
            lhv_data_loaded_provider=lambda: self.lhv_data_loaded,
            fallback_model_provider=lambda: self.calculator.eos,
            invalidate_results=lambda clear_visible_results=True: (
                self.invalidate_results(clear_visible_results)
            ),
            set_result_state=lambda response, result_data: (
                self._set_calculation_result_state(response, result_data)
            ),
            set_warning_messages=lambda messages: (
                self.set_thermo_warning_messages(messages)
            ),
            update_flow_conversion=lambda: self.update_flow_conversion(),
            animate_progress_to=lambda target, duration_ms=500: (
                self.animate_progress_to(target, duration_ms)
            ),
        )

        if not self.lhv_data_loaded:
             print("Warning: MainWindow initialized with no LHV data.")
             # You could disable LHV display or show a warning label here if needed

        # --- Populate widgets ---
        self.populate_comboboxes()
        self.setup_table()
        self.setup_thermo_warning_banner()
        self.setup_flow_tab()

        # --- Connect signals ---
        self.connect_signals()

        # --- Initialize UI State ---
        self.ui.radioButton_mol_percent.setChecked(True) # Default Mol%
        self.on_input_basis_changed() # Set initial table state

        # "Export Report" would be more accurate, but keep "Print Results" as requested.
        if hasattr(self.ui, "printResultsButton"):
            self.ui.printResultsButton.setEnabled(True)

        # Set Window Title (Optional - can also be set in Designer)
        self.setWindowTitle("Thermo Calculator - v.1.0")

    def _parse_float_or_zero(self, text: str) -> float:
        """Parse float from a table cell; treat blanks/invalid as 0.0."""
        return self.composition_table.parse_float_or_zero(text)

    def _recalculate_inactive_column(self):
        """Auto-calculate the inactive composition column from the active one.

        - If Mol% is active, compute Wt% (inactive).
        - If Wt% is active, compute Mol% (inactive).
        The inactive column remains read-only and gray.
        """
        row_count = self.ui.tableWidget.rowCount()
        if row_count <= 1:
            return

        is_mol_basis = self.ui.radioButton_mol_percent.isChecked()
        active_col = 1 if is_mol_basis else 2
        inactive_col = 2 if is_mol_basis else 1
        total_row = row_count - 1

        names = []
        active_vals = []

        for row in range(total_row):
            name_item = self.ui.tableWidget.item(row, 0)
            if not name_item:
                continue
            name = name_item.text()
            names.append(name)

            active_item = self.ui.tableWidget.item(row, active_col)
            active_vals.append(self._parse_float_or_zero(active_item.text() if active_item else ""))

        basis = "Mol %" if is_mol_basis else "Wt %"
        derived = self.derive_composition_use_case.execute(
            DeriveCompositionRequest(
                component_names=tuple(names),
                active_percentages=tuple(active_vals),
                basis=basis,
            )
        ).percentages

        # Write derived values into inactive column (avoid recursion)
        self.ui.tableWidget.blockSignals(True)
        try:
            name_idx = 0
            for row in range(total_row):
                name_item = self.ui.tableWidget.item(row, 0)
                if not name_item:
                    continue

                inactive_item = self.ui.tableWidget.item(row, inactive_col)
                if inactive_item is None:
                    inactive_item = QTableWidgetItem("")
                    self.ui.tableWidget.setItem(row, inactive_col, inactive_item)

                val = derived[name_idx]
                inactive_item.setText("" if val is None else f"{val:.4f}")

                # Keep inactive styling/flags consistent
                inactive_item.setBackground(QColor("lightGray"))
                inactive_item.setFlags(Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable)

                name_idx += 1

            # Keep inactive total cell blank
            total_inactive_item = self.ui.tableWidget.item(total_row, inactive_col)
            if total_inactive_item is not None:
                total_inactive_item.setText("")
                total_inactive_item.setBackground(QColor("lightGray"))
                total_inactive_item.setFlags(Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable)
        finally:
            self.ui.tableWidget.blockSignals(False)

    def invalidate_results(self, clear_visible_results: bool = True):
        """Clear the latest successful result so exports cannot use stale data."""
        self.last_result_data = None
        self.density_actual_kg_m3 = None
        self.density_normal_kg_m3 = None
        self.density_standard_kg_m3 = None
        self.set_thermo_warning_messages([])
        if clear_visible_results:
            self.ui.results_list.clear()
        self.update_flow_conversion()

    def _set_calculation_result_state(self, response, result_data):
        """Store latest calculation state used by export and flow conversion."""
        self.last_result_data = response
        self.density_actual_kg_m3 = result_data.get("density_actual_kg_m3")
        self.density_normal_kg_m3 = result_data.get("density_normal_kg_m3")
        self.density_standard_kg_m3 = result_data.get(
            "density_standard_kg_m3"
        )

    def setup_thermo_warning_banner(self):
        """Create a persistent, non-modal warning banner above the main results area."""
        self.thermo_warning_banner = ThermoWarningBannerController(self.ui)
        self.thermo_warning_label = self.thermo_warning_banner.label

    def set_thermo_warning_messages(self, warning_messages: list[str]):
        """Show or hide the persistent thermo warning banner."""
        if not hasattr(self, "thermo_warning_banner"):
            return
        self.thermo_warning_banner.set_messages(warning_messages)

    def get_export_base_dir(self):
        """Return the directory where generated reports should be saved."""
        return str(self.report_export_controller.get_export_base_dir())

    def build_report_filename(self, timestamp=None):
        """Build the timestamped Excel report filename."""
        return self.report_export_controller.build_report_filename(timestamp)

    def get_current_conditions(self):
        """Compatibility wrapper for report condition collection."""
        return self.report_request_builder.get_current_conditions()

    def _coerce_report_value(self, text: str):
        """Compatibility wrapper for report cell value coercion."""
        return self.report_request_builder.coerce_report_value(text)

    def get_input_composition_rows(self):
        """Compatibility wrapper for report composition row collection."""
        return self.report_request_builder.get_input_composition_rows()

    def build_report_projection(self, result_data):
        """Compatibility wrapper for typed report projection."""
        return self.report_request_builder.build_report_projection(result_data)

    def build_report_warning_rows(self, result_data):
        """Compatibility wrapper for structured report warning rows."""
        return self.report_request_builder.build_report_warning_rows(result_data)

    def build_results_rows(self, result_data):
        """Compatibility wrapper for structured report result rows."""
        return self.report_request_builder.build_results_rows(result_data)

    def export_results_to_excel(self):
        """Export the latest calculation results to a formatted Excel workbook."""
        return self.report_export_controller.export_latest_report(
            has_results=bool(self.last_result_data),
            build_request=self._build_report_export_request,
        )

    def _build_report_export_request(self, report_path, timestamp):
        """Compatibility wrapper for report export request construction."""
        return self.report_request_builder.build_export_request(
            report_path,
            timestamp,
        )

    def setup_flow_tab(self):
        """Initialize the Flow tab controls from the UI adapter."""
        self.flow_tab_controller.setup()

    def _configure_copyable_flow_output(self, line_edit):
        """Compatibility wrapper for the extracted Flow-tab controller."""
        self.flow_tab_controller.configure_copyable_output(line_edit)

    def update_flow_conversion(self):
        """Compatibility wrapper for Flow-tab conversion rendering."""
        self.flow_tab_controller.update_conversion()

    def populate_comboboxes(self):
        """Populate the ComboBox widgets with options."""
        # Component Selection ComboBox
        self.ui.comboBox_select_components.addItems(list(MOLECULAR_WEIGHTS.keys()))

        # Temperature ComboBox
        cryo_temps = [f"{t} °C" for t in range(-250, 1, 10)]  # -250..0 step 10
        existing_temps = ["0 °C", "5 °C", "10 °C", "15 °C", "25 °C", "50 °C", "100 °C", "150 °C", "200 °C"]
        temps = []
        for t in cryo_temps + existing_temps:
            if t not in temps:
                temps.append(t)
        self.ui.comboBox_select_temperature.addItems(temps)
        self.ui.comboBox_select_temperature.setCurrentText("0 °C") # Default

        # Pressure ComboBox
        self.ui.comboBox_select_pressure.addItems(["1 atm", "5 atm", "10 atm", "20 atm", "50 atm", "100 atm"])
        self.ui.comboBox_select_pressure.setCurrentText("1 atm") # Default

        # EOS ComboBox
        self.ui.comboBox_select_EOS.addItems(["PRMIX"]) # Only PRMIX for now
        # Add more if the thermo gateway supports them, e.g. ["PRMIX", "SRKMIX"].


    def setup_table(self):
        """Set up the composition table headers and initial row."""
        self.composition_table.setup_table()


    def connect_signals(self):
        """Connect UI element signals to methods."""
        self.ui.comboBox_select_components.currentIndexChanged.connect(self.add_component)
        self.ui.comboBox_select_EOS.currentIndexChanged.connect(self.update_eos)

        # Connect Radio Buttons
        self.ui.radioButton_mol_percent.toggled.connect(self.on_input_basis_changed)
        # Wt% toggled is implicitly handled by Mol% toggling if in a button group,
        # but connecting both is safer if group isn't used or for explicit clarity.
        self.ui.radioButton_wt_percent.toggled.connect(self.on_input_basis_changed)

        # Connect Buttons
        self.ui.remove_component_button.clicked.connect(self.remove_component)
        self.ui.clear_all_button.clicked.connect(self.clear_all)
        self.ui.delete_all_button.clicked.connect(self.clear_all) # Connect the second clear button
        self.ui.go_button.clicked.connect(self.calculate_and_display)
        self.ui.pushButton_normalize.clicked.connect(self.normalize_composition)
        if hasattr(self.ui, "printResultsButton"):
            self.ui.printResultsButton.clicked.connect(self.export_results_to_excel)
        self.flow_tab_controller.connect_signals()

        # Connect Table Change Signal
        self.ui.tableWidget.itemChanged.connect(self.on_table_item_changed)

        # Optional: Group Radio Buttons programmatically (good practice)
        self.basis_group = QButtonGroup(self) # Create group associated with MainWindow
        self.basis_group.addButton(self.ui.radioButton_mol_percent)
        self.basis_group.addButton(self.ui.radioButton_wt_percent)

        # Step 3: Connect temperature and pressure combo boxes to clear results_list
        self.ui.comboBox_select_temperature.currentIndexChanged.connect(self.on_temperature_or_pressure_changed)
        self.ui.comboBox_select_pressure.currentIndexChanged.connect(self.on_temperature_or_pressure_changed)


    def update_eos(self):
        """Update the EOS in the calculator when the ComboBox changes."""
        # Step 5: Reset progress bar to 0 on EOS change
        if hasattr(self.ui, 'progressBar'):
            self.ui.progressBar.setValue(0)
        self.invalidate_results()
        selected_eos = self.ui.comboBox_select_EOS.currentText()
        self.calculator.set_eos(selected_eos)

    def add_component(self):
        """Add the selected component to the list and table."""
        # Step 5: Reset progress bar to 0 on add component
        if hasattr(self.ui, 'progressBar'):
            self.ui.progressBar.setValue(0)
        result = self.composition_table.add_selected_component()
        if result.status is ComponentAddStatus.EMPTY_SELECTION:
            return
        if result.status is ComponentAddStatus.DUPLICATE:
            QMessageBox.warning(
                self,
                "Duplicate",
                f"Component '{result.component_name}' is already selected.",
            )
            self.composition_table.reset_component_selection()
            return

        # Re-apply basis logic to set correct editability/colors for the new row
        self.on_input_basis_changed() # This will handle flags and background

        # Recalculate totals
        self.update_table_total()

        # Reset combo index to dummy/placeholder
        self.composition_table.reset_component_selection()

        # Clear results after successful add because prior calculations are now stale.
        self.invalidate_results()

    def remove_component(self):
        """Remove the selected component from the list and table."""
        # Step 5: Reset progress bar to 0 on remove component
        if hasattr(self.ui, 'progressBar'):
            self.ui.progressBar.setValue(0)
        result = self.composition_table.remove_selected_components()
        if result.status is ComponentRemoveStatus.NO_COMPONENTS:
            QMessageBox.information(self, "Remove", "No components to remove.")
            return

        self.update_table_total()

        # Clear results after successful remove because prior calculations are now stale.
        self.invalidate_results()

    def _remove_component_from_table(self, comp_name):
        """Helper to remove the row matching comp_name from the table."""
        self.composition_table.remove_component_from_table(comp_name)

    def clear_all(self):
        """Clear all inputs, selections, and results."""
        # Step 5: Reset progress bar to 0 on clear all
        if hasattr(self.ui, 'progressBar'):
            self.ui.progressBar.setValue(0)
        # Clear results at the start because all prior calculations are now stale.
        self.invalidate_results()
        self.composition_table.clear_component_rows_and_selection()

        # Clear results list
        # self.ui.results_list.clear() # Already cleared above

        # Reset calculator
        self.calculator.set_components({})

        # Reset UI states
        self.ui.go_button.setEnabled(True) # Re-enable Go button
        # Reset T/P/EOS combos to defaults maybe?
        self.ui.comboBox_select_temperature.setCurrentText("0 °C")
        self.ui.comboBox_select_pressure.setCurrentText("1 atm")
        self.ui.comboBox_select_EOS.setCurrentIndex(0)

        # Reset basis to Mol %
        self.ui.radioButton_mol_percent.setChecked(True)
        # on_input_basis_changed will be triggered by setChecked, which calls update_table_total

        # Explicitly update total fields after clearing
        self.update_table_total()

    def on_input_basis_changed(self):
        """Update table column appearance and editability based on radio buttons."""
        self.composition_table.apply_basis_state()
        # Recalculate derived (inactive) column values and then totals for the active column
        self._recalculate_inactive_column()
        self.invalidate_results()
        self.update_table_total() # Update totals after changing basis

    def on_table_item_changed(self, item):
        """Recalculate total when a percentage value changes."""
        # Step 5: Reset progress bar to 0 on composition change
        if hasattr(self.ui, 'progressBar'):
            self.ui.progressBar.setValue(0)
        # Only update if the changed item is in an active percentage column and not the total row
        active_col = 1 if self.ui.radioButton_mol_percent.isChecked() else 2
        if item.column() == active_col and item.row() != self.ui.tableWidget.rowCount() - 1:
             self.invalidate_results()
             # Basic validation: try converting to float, set red background if invalid?
             try:
                 float(item.text())
                 item.setBackground(QColor("white")) # Reset background on valid input
             except ValueError:
                 if item.text().strip(): # Don't color empty cells red
                      item.setBackground(QColor("pink")) # Indicate invalid number format
             # Live update of the derived inactive column
             self._recalculate_inactive_column()
             self.update_table_total()

    def update_table_total(self):
        """Sum the percentages in the active column and update the Total cell."""
        return self.composition_table.update_table_total()

    def calculate_and_display(self):
        """Compatibility wrapper for calculation workflow orchestration."""
        return self.calculation_workflow_controller.calculate_and_display()

    def on_calculation_result(self, result_data):
        """Compatibility wrapper for calculation-result rendering."""
        return self.calculation_workflow_controller.on_calculation_result(
            result_data
        )

    def on_calculation_error(self, error_msg):
        """Compatibility wrapper for calculation-error rendering."""
        return self.calculation_workflow_controller.on_calculation_error(
            error_msg
        )

    def on_calculation_finished(self):
        """Compatibility wrapper for calculation-finished UI state."""
        return self.calculation_workflow_controller.on_calculation_finished()

    def normalize_composition(self):
        """Normalize the active composition column (mol% or wt%) so the sum is 100%."""
        values = self.composition_table.read_active_percentages_for_normalization()
        try:
            response = self.normalize_composition_use_case.execute(
                NormalizeCompositionRequest(percentages=values)
            )
        except ValueError:
            QMessageBox.warning(self, "Normalize", "Cannot normalize: total is zero.")
            return

        self.composition_table.write_normalized_active_percentages(
            response.percentages
        )
        self.invalidate_results()
        self.update_table_total()
        self.composition_table.show_normalization_success()

    def on_temperature_or_pressure_changed(self):
        """Clear the results_list when temperature or pressure is changed."""
        # Step 5: Reset progress bar to 0 on T/P change
        if hasattr(self.ui, 'progressBar'):
            self.ui.progressBar.setValue(0)
        self.invalidate_results()

    def clear_results_list(self):
        """Clear the results_list widget."""
        self.ui.results_list.clear()

    def animate_progress_to(self, target_value, duration_ms=500):
        """Animate the progress bar to the target value over the given duration (ms)."""
        if not hasattr(self.ui, 'progressBar'):
            return
        # Stop any existing timer
        if self.progress_timer is not None:
            self.progress_timer.stop()
            self.progress_timer.deleteLater()
            self.progress_timer = None
        # Set up animation state
        self.progress_target = target_value
        self.progress_animation_duration = duration_ms
        self.progress_animation_start_value = self.ui.progressBar.value()
        self.progress_animation_start_time = QTimer().remainingTime()  # Not used, will use QElapsedTimer below

        # Use QElapsedTimer to track elapsed time
        from PyQt6.QtCore import QElapsedTimer
        self._progress_animation_timer = QElapsedTimer()
        self._progress_animation_timer.start()

        def update_progress():
            elapsed = self._progress_animation_timer.elapsed()
            if elapsed >= self.progress_animation_duration:
                self.ui.progressBar.setValue(self.progress_target)
                if self.progress_timer is not None:
                    self.progress_timer.stop()
                    self.progress_timer.deleteLater()
                    self.progress_timer = None
                return
            # Linear interpolation
            start = self.progress_animation_start_value
            end = self.progress_target
            value = start + (end - start) * (elapsed / self.progress_animation_duration)
            self.ui.progressBar.setValue(int(value))

        self.progress_timer = QTimer(self)
        self.progress_timer.timeout.connect(update_progress)
        self.progress_timer.start(20)
