import sys
import os # To check if DB file exists
from datetime import datetime
from pathlib import Path

from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QComboBox, QListWidget,
    QPushButton, QLabel, QGridLayout, QVBoxLayout, QTableWidget,
    QTableWidgetItem, QMessageBox, QHeaderView, QHBoxLayout, QButtonGroup, QRadioButton, QStyleFactory
)
from PyQt6.QtCore import Qt, QTimer, QObject, pyqtSignal, QThread
from PyQt6.QtGui import QPalette, QColor

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
from thermo_components.domain.flow_conversion import (
    convert_flow,
    format_flow_value,
    parse_flow_input,
)
from thermo_components.domain.flow_units import (
    FLOW_UNIT_DEFINITIONS,
    FLOW_UNIT_ORDER,
    FT3_PER_M3,
    HOURS_PER_DAY,
    KG_PER_LB,
    KG_PER_TONNE,
    LB_PER_KLB,
    M3_PER_BBL,
    M3_PER_FT3,
)
from thermo_components.domain.lhv import (
    KCAL_PER_MJ,
    MJ_PER_MMBTU,
    NORMAL_MOLAR_VOLUME_NM3_PER_KMOL,
    build_lhv_display_values,
    format_lhv_display_value,
)
from thermo_components.domain.results import (
    build_density_note,
    extract_scalar_density_value,
)
from thermo_components.adapters.thermo import ThermoGateway
from thermo_components.adapters.packaging import RuntimeResourceLocator
from thermo_components.adapters.persistence import SqliteLhvRepository
from thermo_components.adapters.reporting import OpenPyxlReportExporter
from thermo_components.application.dto import (
    DeriveCompositionRequest,
    FlowConversionRequest,
    NormalizeCompositionRequest,
    ReportCompositionRow,
    ReportConditionRow,
    ReportExportRequest,
    PropertyCalculationRequest,
    ReportPreparationRequest,
    coerce_property_response,
)
from thermo_components.application.use_cases import (
    CalculatePropertiesUseCase,
    ConvertFlowUseCase,
    DeriveCompositionUseCase,
    NormalizeCompositionUseCase,
    PrepareReportUseCase,
)

MixtureCalculator = ThermoGateway

def load_lhv_data(db_path='lhv_data.db'):
    """Load LHV data through the SQLite persistence adapter."""
    return SqliteLhvRepository(db_path).load_all()

# Helper to find resource files in both dev and PyInstaller bundle
def resource_path(relative_path):
    """Resolve a resource in source mode or a PyInstaller bundle."""
    return str(RuntimeResourceLocator().resolve(relative_path))


# Step 1: Worker class for threaded calculation
class CalculationWorker(QObject):
    result = pyqtSignal(object)
    error = pyqtSignal(str)
    finished = pyqtSignal()

    def __init__(self, use_case, request):
        super().__init__()
        self.use_case = use_case
        self.request = request

    def run(self):
        try:
            self.result.emit(self.use_case.execute(self.request))
        except ValueError as exc:
            self.error.emit(f"Error: {exc}")
        except Exception as e:
            import traceback
            self.error.emit(f"Worker Exception: {e}\n{traceback.format_exc()}")
        finally:
            self.finished.emit()

# --- MainWindow Class (Modified for gui.py) ---
class MainWindow(QMainWindow):
    def __init__(self, lhv_data):
        super().__init__()

        # --- Set up the UI from Designer ---
        self.ui = Ui_Dialog()
        self.ui.setupUi(self) # Setup the UI onto this QMainWindow
        if hasattr(self.ui, "tabWidget"):
            self.ui.tabWidget.setCurrentIndex(0)

        # Step 2: Progress bar animation state
        self.progress_timer = None
        self.progress_target = 100
        self.progress_animation_duration = 800  # milliseconds
        self.progress_animation_start_time = 0
        self.progress_animation_start_value = 0

        # --- Store LHV data and calculator instance ---
        self.lhv_database = lhv_data
        self.lhv_data_loaded = bool(self.lhv_database)
        self.calculator = MixtureCalculator()
        self.calculate_properties_use_case = CalculatePropertiesUseCase(
            self.calculator,
            self.lhv_database,
        )
        self.convert_flow_use_case = ConvertFlowUseCase()
        self.derive_composition_use_case = DeriveCompositionUseCase()
        self.normalize_composition_use_case = NormalizeCompositionUseCase()
        self.prepare_report_use_case = PrepareReportUseCase()
        self.report_exporter = OpenPyxlReportExporter()
        self.last_result_data = None
        self.density_actual_kg_m3 = None
        self.density_normal_kg_m3 = None
        self.density_standard_kg_m3 = None

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
        try:
            return float(text) if text and str(text).strip() else 0.0
        except (TypeError, ValueError):
            return 0.0

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

    def setup_thermo_warning_banner(self):
        """Create a persistent, non-modal warning banner above the main results area."""
        self.thermo_warning_label = QLabel(self.ui.tab)
        self.thermo_warning_label.setObjectName("thermoWarningLabel")
        self.thermo_warning_label.setWordWrap(True)
        self.thermo_warning_label.setTextFormat(Qt.TextFormat.PlainText)
        self.thermo_warning_label.setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop)
        self.thermo_warning_label.setMargin(8)
        self.thermo_warning_label.setStyleSheet(
            "QLabel#thermoWarningLabel {"
            "background-color: rgb(255, 244, 204);"
            "color: rgb(122, 63, 0);"
            "border: 1px solid rgb(230, 184, 0);"
            "border-radius: 4px;"
            "}"
        )
        self.thermo_warning_label.hide()

        self._results_label_base_geometry = self.ui.results_label.geometry()
        self._results_list_base_geometry = self.ui.results_list.geometry()
        self._progress_bar_base_geometry = self.ui.progressBar.geometry() if hasattr(self.ui, "progressBar") else None
        self._warning_banner_y = self.ui.groupBox.geometry().bottom() + 8

    def set_thermo_warning_messages(self, warning_messages: list[str]):
        """Show or hide the persistent thermo warning banner."""
        if not hasattr(self, "thermo_warning_label"):
            return

        cleaned_messages = [str(message).strip() for message in warning_messages if str(message).strip()]
        if not cleaned_messages:
            self.thermo_warning_label.hide()
            self.ui.results_label.setGeometry(self._results_label_base_geometry)
            self.ui.results_list.setGeometry(self._results_list_base_geometry)
            if self._progress_bar_base_geometry is not None:
                self.ui.progressBar.setGeometry(self._progress_bar_base_geometry)
            return

        warning_text = "\n".join(cleaned_messages)
        banner_width = self._results_list_base_geometry.width()
        text_rect = self.thermo_warning_label.fontMetrics().boundingRect(
            0,
            0,
            banner_width - 16,
            1000,
            int(Qt.TextFlag.TextWordWrap),
            warning_text,
        )
        banner_height = max(44, text_rect.height() + 16)
        self.thermo_warning_label.setText(warning_text)
        self.thermo_warning_label.setGeometry(
            self._results_list_base_geometry.x(),
            self._warning_banner_y,
            banner_width,
            banner_height,
        )
        self.thermo_warning_label.show()
        self.thermo_warning_label.raise_()

        results_label_y = self._warning_banner_y + banner_height + 8
        self.ui.results_label.setGeometry(
            self._results_label_base_geometry.x(),
            results_label_y,
            self._results_label_base_geometry.width(),
            self._results_label_base_geometry.height(),
        )

        list_y_offset = self._results_list_base_geometry.y() - self._results_label_base_geometry.y()
        results_list_y = results_label_y + list_y_offset
        base_results_bottom = self._results_list_base_geometry.y() + self._results_list_base_geometry.height()
        results_list_height = max(100, base_results_bottom - results_list_y)
        self.ui.results_list.setGeometry(
            self._results_list_base_geometry.x(),
            results_list_y,
            self._results_list_base_geometry.width(),
            results_list_height,
        )

        if self._progress_bar_base_geometry is not None:
            self.ui.progressBar.setGeometry(
                self._progress_bar_base_geometry.x(),
                base_results_bottom + 8,
                self._progress_bar_base_geometry.width(),
                self._progress_bar_base_geometry.height(),
            )

    def get_export_base_dir(self):
        """Return the directory where generated reports should be saved."""
        if getattr(sys, "frozen", False):
            return os.path.dirname(sys.executable)
        return os.path.dirname(os.path.abspath(__file__))

    def build_report_filename(self, timestamp: datetime | None = None):
        """Build the timestamped Excel report filename."""
        timestamp = timestamp or datetime.now()
        return f"Thermo_Report_{timestamp.strftime('%Y%m%d_%H%M%S')}.xlsx"

    def get_current_conditions(self):
        """Return the currently active calculation conditions for reporting."""
        if self.last_result_data:
            response = coerce_property_response(self.last_result_data)
            basis = response.basis
            model_display = response.model_display
        else:
            basis = (
                "Mol %"
                if self.ui.radioButton_mol_percent.isChecked()
                else "Wt %"
            )
            model_display = self.ui.comboBox_select_EOS.currentText()
        return [
            ("Basis", basis),
            ("Temperature", self.ui.comboBox_select_temperature.currentText()),
            ("Pressure", self.ui.comboBox_select_pressure.currentText()),
            ("Model / EOS", model_display),
        ]

    def _coerce_report_value(self, text: str):
        """Convert report cell text to float when possible, preserving blanks/text."""
        cleaned = str(text).strip() if text is not None else ""
        if not cleaned:
            return ""
        try:
            return float(cleaned)
        except ValueError:
            return cleaned

    def get_input_composition_rows(self):
        """Return the input composition table rows, including the Total row."""
        rows = []
        for row_index in range(self.ui.tableWidget.rowCount()):
            component_item = self.ui.tableWidget.item(row_index, 0)
            mol_item = self.ui.tableWidget.item(row_index, 1)
            wt_item = self.ui.tableWidget.item(row_index, 2)
            rows.append({
                "Component": component_item.text() if component_item else "",
                "Mol %": self._coerce_report_value(mol_item.text() if mol_item else ""),
                "Wt %": self._coerce_report_value(wt_item.text() if wt_item else ""),
            })
        return rows

    def build_report_projection(self, result_data):
        """Build the typed report projection used by display and export paths."""
        response = coerce_property_response(result_data)
        return self.prepare_report_use_case.execute(
            ReportPreparationRequest(
                calculation=response,
                lhv_data_available=self.lhv_data_loaded,
            )
        )

    def build_report_warning_rows(self, result_data):
        """Build structured warning rows for the export report."""
        projection = self.build_report_projection(result_data)
        return [row.to_dict() for row in projection.warning_rows]

    def build_results_rows(self, result_data):
        """Build structured result rows for report export without scraping the UI."""
        projection = self.build_report_projection(result_data)
        return [row.to_dict() for row in projection.result_rows]

    def export_results_to_excel(self):
        """Export the latest calculation results to a formatted Excel workbook."""
        if not self.last_result_data:
            QMessageBox.information(self, "Print Results", "No calculation results available. Please click Go first.")
            return None

        timestamp = datetime.now()
        report_path = Path(
            self.get_export_base_dir(),
            self.build_report_filename(timestamp),
        )
        export_request = ReportExportRequest(
            report_path=report_path,
            exported_at=timestamp,
            conditions=tuple(
                ReportConditionRow(setting, value)
                for setting, value in self.get_current_conditions()
            ),
            input_rows=tuple(
                ReportCompositionRow(
                    row["Component"],
                    row["Mol %"],
                    row["Wt %"],
                )
                for row in self.get_input_composition_rows()
            ),
            projection=self.build_report_projection(self.last_result_data),
        )
        try:
            exported_path = self.report_exporter.export(export_request)
        except ImportError:
            QMessageBox.warning(
                self,
                "Print Results",
                "Excel export requires the 'openpyxl' package.\n"
                "Install it in this environment and try again.\n\n"
                "Command:\npython -m pip install openpyxl",
            )
            return None
        except Exception as exc:
            QMessageBox.critical(self, "Print Results", f"Failed to save Excel report:\n{exc}")
            return None

        QMessageBox.information(self, "Print Results", f"Excel report saved to:\n{exported_path}")

        if hasattr(os, "startfile"):
            try:
                os.startfile(exported_path)
            except OSError:
                pass

        return str(exported_path)

    def setup_flow_tab(self):
        """Initialize the Flow tab controls from the current UI definition."""
        if not hasattr(self.ui, "comboBox_select_units"):
            return

        self.ui.comboBox_select_units.clear()
        self.ui.comboBox_select_units.addItems(FLOW_UNIT_ORDER)
        self.ui.comboBox_select_desired_units.clear()
        self.ui.comboBox_select_desired_units.addItems(FLOW_UNIT_ORDER)
        self.ui.comboBox_select_units.setCurrentText("kg/h")
        self.ui.comboBox_select_desired_units.setCurrentText("t/d")
        self._configure_copyable_flow_output(self.ui.lineEdit_result)
        self._configure_copyable_flow_output(self.ui.lineEdit_in_desired_units)
        self.update_flow_conversion()

    def _configure_copyable_flow_output(self, line_edit):
        """Keep flow outputs read-only while preserving selection and copy behavior."""
        line_edit.setReadOnly(True)
        line_edit.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        line_edit.setContextMenuPolicy(Qt.ContextMenuPolicy.DefaultContextMenu)
        line_edit.setCursor(Qt.CursorShape.IBeamCursor)
        line_edit.setDragEnabled(True)

    def update_flow_conversion(self):
        """Recalculate the Flow-tab conversion using the latest density state."""
        if not hasattr(self.ui, "lineEdit_enter_flow"):
            return

        from_unit = self.ui.comboBox_select_units.currentText().strip()
        to_unit = self.ui.comboBox_select_desired_units.currentText().strip()
        self.ui.lineEdit_in_desired_units.setText(to_unit)

        if not from_unit or not to_unit:
            self.ui.lineEdit_result.clear()
            return

        input_text = self.ui.lineEdit_enter_flow.text()
        try:
            response = self.convert_flow_use_case.execute(
                FlowConversionRequest(
                    input_text=input_text,
                    from_unit=from_unit,
                    to_unit=to_unit,
                    normal_density_kg_m3=self.density_normal_kg_m3,
                    standard_density_kg_m3=self.density_standard_kg_m3,
                )
            )
        except ValueError as exc:
            self.ui.lineEdit_result.setText(str(exc))
            return

        self.ui.lineEdit_result.setText(response.display_value)

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
        # Add more if MixtureCalculator supports them, e.g., ["PRMIX", "SRKMIX"]


    def setup_table(self):
        """Set up the composition table headers and initial row."""
        self.ui.tableWidget.setColumnCount(3)
        self.ui.tableWidget.setHorizontalHeaderLabels(["Component", "Mol %", "Wt %"])
        # Set reasonable column widths (adjust as needed)
        self.ui.tableWidget.setColumnWidth(0, 150)
        self.ui.tableWidget.setColumnWidth(1, 80)
        self.ui.tableWidget.setColumnWidth(2, 80)
        header = self.ui.tableWidget.horizontalHeader()
        header.setSectionResizeMode(0, QHeaderView.ResizeMode.Interactive) # Component name resizable
        header.setSectionResizeMode(1, QHeaderView.ResizeMode.Stretch) # Stretch % columns
        header.setSectionResizeMode(2, QHeaderView.ResizeMode.Stretch)

        # Add the non-editable "Total" row
        self.ui.tableWidget.setRowCount(1)
        total_item = QTableWidgetItem("Total")
        total_item.setFlags(Qt.ItemFlag.ItemIsEnabled) # Make non-editable
        self.ui.tableWidget.setItem(0, 0, total_item)
        # Initialize Total value cells
        zero_item1 = QTableWidgetItem("0.00")
        zero_item1.setFlags(Qt.ItemFlag.ItemIsEnabled)
        self.ui.tableWidget.setItem(0, 1, zero_item1)
        zero_item2 = QTableWidgetItem("0.00")
        zero_item2.setFlags(Qt.ItemFlag.ItemIsEnabled)
        self.ui.tableWidget.setItem(0, 2, zero_item2)


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
        if hasattr(self.ui, "lineEdit_enter_flow"):
            self.ui.lineEdit_enter_flow.textChanged.connect(self.update_flow_conversion)
            self.ui.comboBox_select_units.currentIndexChanged.connect(self.update_flow_conversion)
            self.ui.comboBox_select_desired_units.currentIndexChanged.connect(self.update_flow_conversion)

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
        component = self.ui.comboBox_select_components.currentText()

        # Don't add the dummy entry or duplicates
        if not component.strip():
            return # Ignore the dummy empty selection

        # Check for duplicates in the list_widget (selected_components_list)
        current_items = [self.ui.selected_components_list.item(i).text() for i in range(self.ui.selected_components_list.count())]
        if component in current_items:
            # Highlight combo red and disable Go (or show message box)
            palette = self.ui.comboBox_select_components.palette()
            palette.setColor(QPalette.ColorRole.Base, QColor('red'))
            self.ui.comboBox_select_components.setPalette(palette)
            # self.ui.go_button.setEnabled(False) # Disabling Go might be too strict
            QMessageBox.warning(self, "Duplicate", f"Component '{component}' is already selected.")
            # Reset combo index to dummy after showing message
            self.ui.comboBox_select_components.setCurrentIndex(0)
            return

        # Reset combo box color if previously red
        palette = self.ui.comboBox_select_components.palette()
        palette.setColor(QPalette.ColorRole.Base, QColor('white'))
        self.ui.comboBox_select_components.setPalette(palette)

        # Add to the list_widget
        self.ui.selected_components_list.addItem(component)

        # Add to the tableWidget
        total_row_index = self.ui.tableWidget.rowCount() - 1
        self.ui.tableWidget.insertRow(total_row_index)
        # Put the component name in column 0 (non-editable)
        name_item = QTableWidgetItem(component)
        name_item.setFlags(name_item.flags() & ~Qt.ItemFlag.ItemIsEditable) # Make name non-editable
        self.ui.tableWidget.setItem(total_row_index, 0, name_item)

        # Add empty items for Mol% and Wt% columns
        mol_item = QTableWidgetItem("") # Start empty
        wt_item = QTableWidgetItem("") # Start empty
        self.ui.tableWidget.setItem(total_row_index, 1, mol_item)
        self.ui.tableWidget.setItem(total_row_index, 2, wt_item)

        # Re-apply basis logic to set correct editability/colors for the new row
        self.on_input_basis_changed() # This will handle flags and background

        # Recalculate totals
        self.update_table_total()

        # Reset combo index to dummy/placeholder
        self.ui.comboBox_select_components.setCurrentIndex(0)

        # Clear results after successful add because prior calculations are now stale.
        self.invalidate_results()

    def remove_component(self):
        """Remove the selected component from the list and table."""
        # Step 5: Reset progress bar to 0 on remove component
        if hasattr(self.ui, 'progressBar'):
            self.ui.progressBar.setValue(0)
        selected_list_items = self.ui.selected_components_list.selectedItems()

        if not selected_list_items:
            # If nothing selected, maybe remove the last one added? Or show message.
            if self.ui.selected_components_list.count() > 0:
                 current_row = self.ui.selected_components_list.currentRow()
                 if current_row < 0: # If nothing is actively selected, use last item
                      current_row = self.ui.selected_components_list.count() - 1
                 item_to_remove = self.ui.selected_components_list.takeItem(current_row) # Remove from list
                 if item_to_remove:
                      self._remove_component_from_table(item_to_remove.text()) # Remove from table
            else:
                 QMessageBox.information(self, "Remove", "No components to remove.")
                 return
        else:
            # Remove all selected items
            for item in selected_list_items:
                row = self.ui.selected_components_list.row(item)
                comp_name = item.text()
                self.ui.selected_components_list.takeItem(row) # Remove from list
                self._remove_component_from_table(comp_name) # Remove from table by name

        self.update_table_total()
        # Reset combo color if it was red
        if self.ui.selected_components_list.count() == 0:
            palette = self.ui.comboBox_select_components.palette()
            palette.setColor(QPalette.ColorRole.Base, QColor('white'))
            self.ui.comboBox_select_components.setPalette(palette)

        # Clear results after successful remove because prior calculations are now stale.
        self.invalidate_results()

    def _remove_component_from_table(self, comp_name):
        """Helper to remove the row matching comp_name from the table."""
        # Iterate backwards to avoid index issues when removing rows
        for row_index in range(self.ui.tableWidget.rowCount() - 2, -1, -1): # Stop before Total row
            table_item = self.ui.tableWidget.item(row_index, 0)
            if table_item and table_item.text() == comp_name:
                self.ui.tableWidget.removeRow(row_index)
                break # Assume component names are unique

    def clear_all(self):
        """Clear all inputs, selections, and results."""
        # Step 5: Reset progress bar to 0 on clear all
        if hasattr(self.ui, 'progressBar'):
            self.ui.progressBar.setValue(0)
        # Clear results at the start because all prior calculations are now stale.
        self.invalidate_results()
        # Clear list
        self.ui.selected_components_list.clear()

        # Clear table (remove all rows except the header and 'Total' row)
        # Iterate backwards from second-to-last row up to first row
        while self.ui.tableWidget.rowCount() > 1:
             self.ui.tableWidget.removeRow(0) # Remove the first component row repeatedly

        # Clear results list
        # self.ui.results_list.clear() # Already cleared above

        # Reset calculator
        self.calculator.set_components({})

        # Reset UI states
        self.ui.go_button.setEnabled(True) # Re-enable Go button
        palette = self.ui.comboBox_select_components.palette()
        palette.setColor(QPalette.ColorRole.Base, QColor('white'))
        self.ui.comboBox_select_components.setPalette(palette)
        self.ui.comboBox_select_components.setCurrentIndex(0) # Reset selection
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
        is_mol_basis = self.ui.radioButton_mol_percent.isChecked()
        active_col = 1 if is_mol_basis else 2
        inactive_col = 2 if is_mol_basis else 1

        # Block signals temporarily to prevent itemChanged triggers during modification
        self.ui.tableWidget.blockSignals(True)

        row_count = self.ui.tableWidget.rowCount()
        total_row = row_count - 1 # Index of the 'Total' row

        for row in range(row_count):
            active_item = self.ui.tableWidget.item(row, active_col)
            inactive_item = self.ui.tableWidget.item(row, inactive_col)

            # Ensure items exist (might not for newly added row before full setup)
            if active_item is None:
                active_item = QTableWidgetItem("")
                self.ui.tableWidget.setItem(row, active_col, active_item)
            if inactive_item is None:
                inactive_item = QTableWidgetItem("")
                self.ui.tableWidget.setItem(row, inactive_col, inactive_item)

            # Set background colors
            active_item.setBackground(QColor("white"))
            inactive_item.setBackground(QColor("lightGray"))

            # Set flags (editability)
            if row != total_row: # Component rows
                active_item.setFlags(Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable | Qt.ItemFlag.ItemIsEditable)
                inactive_item.setFlags(Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable) # Not editable
            else: # Total row
                active_item.setFlags(Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable) # Not editable
                inactive_item.setFlags(Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable) # Not editable

        # Unblock signals
        self.ui.tableWidget.blockSignals(False)
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
        active_col = 1 if self.ui.radioButton_mol_percent.isChecked() else 2
        row_count = self.ui.tableWidget.rowCount()
        total_row_index = row_count - 1
        total_percent = 0.0
        valid_input = True

        # Iterate only component rows (up to total_row_index)
        for r in range(total_row_index):
            percent_item = self.ui.tableWidget.item(r, active_col)
            if percent_item:
                try:
                    value = float(percent_item.text()) if percent_item.text().strip() else 0.0
                    total_percent += value
                except ValueError:
                     # If text isn't a valid float, treat as error for sum validation
                     if percent_item.text().strip(): # Only count as error if not empty
                          valid_input = False
                     # Don't add to total, let the 100% check fail

        # Update the total cell in the active column
        total_item_active = self.ui.tableWidget.item(total_row_index, active_col)
        if not total_item_active:
            total_item_active = QTableWidgetItem()
            self.ui.tableWidget.setItem(total_row_index, active_col, total_item_active)
        total_item_active.setText(f"{total_percent:.4f}")
        total_item_active.setFlags(Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable) # Ensure non-editable

        # Clear the total cell in the inactive column
        inactive_col = 2 if active_col == 1 else 1
        total_item_inactive = self.ui.tableWidget.item(total_row_index, inactive_col)
        if not total_item_inactive:
             total_item_inactive = QTableWidgetItem()
             self.ui.tableWidget.setItem(total_row_index, inactive_col, total_item_inactive)
        total_item_inactive.setText("") # Clear inactive total
        total_item_inactive.setFlags(Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable)

        # Check if total is 100% and input is valid
        is_total_ok = abs(total_percent - 100.0) < 1e-4 # Use tolerance for float comparison
        if is_total_ok and valid_input:
            total_item_active.setForeground(QColor('black'))
            self.ui.go_button.setEnabled(True)
        else:
            total_item_active.setForeground(QColor('red'))
            self.ui.go_button.setEnabled(False)

    def calculate_and_display(self):
        """Gather inputs, run calculations, and display results."""
        # Step 3: Set progress bar to 0 at the start of calculation
        if hasattr(self.ui, 'progressBar'):
            self.ui.progressBar.setValue(0)
        # Start animating progress bar to 80% over 1 second as soon as Go! is clicked
        self.animate_progress_to(80, 1000)
        self.invalidate_results()

        # --- Gather inputs for worker ---
        comp_names = []
        mol_percents = []
        wt_percents = []
        row_count = self.ui.tableWidget.rowCount() - 1 # Exclude total row
        valid_input = True
        for r in range(row_count):
            comp_item = self.ui.tableWidget.item(r, 0)
            mol_item = self.ui.tableWidget.item(r, 1)
            wt_item = self.ui.tableWidget.item(r, 2)
            if not comp_item: continue
            comp_names.append(comp_item.text())
            try:
                mol_val = float(mol_item.text()) if mol_item and mol_item.text().strip() else 0.0
                mol_percents.append(mol_val)
            except ValueError:
                mol_percents.append(0.0)
                if self.ui.radioButton_mol_percent.isChecked(): valid_input = False
            try:
                wt_val = float(wt_item.text()) if wt_item and wt_item.text().strip() else 0.0
                wt_percents.append(wt_val)
            except ValueError:
                wt_percents.append(0.0)
                if self.ui.radioButton_wt_percent.isChecked(): valid_input = False

        if not valid_input:
            self.ui.results_list.addItem("Error: Invalid numeric input in composition table.")
            self.animate_progress_to(100, 500)
            return
        if not comp_names:
            self.ui.results_list.addItem("Error: No components selected.")
            self.animate_progress_to(100, 500)
            return

        basis = "Mol %" if self.ui.radioButton_mol_percent.isChecked() else "Wt %"
        try:
            temp_text = self.ui.comboBox_select_temperature.currentText()
            T_c = float(temp_text.split("°")[0])
            T_k = T_c + 273.15
        except ValueError:
            self.ui.results_list.addItem("Error: Invalid Temperature selection.")
            self.animate_progress_to(100, 500)
            return
        try:
            pressure_text = self.ui.comboBox_select_pressure.currentText()
            pressure_atm = float(pressure_text.split(" ")[0])
            pressure_pa = pressure_atm * 101325.0
            if pressure_pa <= 0: raise ValueError("Pressure must be positive")
        except ValueError:
            self.ui.results_list.addItem("Error: Invalid Pressure selection.")
            self.animate_progress_to(100, 500)
            return

        # --- Threaded calculation setup ---
        self.ui.go_button.setEnabled(False)
        if hasattr(self.ui, "printResultsButton"):
            self.ui.printResultsButton.setEnabled(False)
        calculation_request = PropertyCalculationRequest.from_sequences(
            component_names=comp_names,
            mole_percents=mol_percents,
            weight_percents=wt_percents,
            basis=basis,
            temperature_k=T_k,
            pressure_pa=pressure_pa,
            pressure_atm=pressure_atm,
        )
        self.worker_thread = QThread()
        self.worker = CalculationWorker(
            use_case=self.calculate_properties_use_case,
            request=calculation_request,
        )
        self.worker.moveToThread(self.worker_thread)
        self.worker.result.connect(self.on_calculation_result)
        self.worker.error.connect(self.on_calculation_error)
        self.worker.finished.connect(self.on_calculation_finished)
        self.worker_thread.started.connect(self.worker.run)
        self.worker.finished.connect(self.worker_thread.quit)
        self.worker.finished.connect(self.worker.deleteLater)
        self.worker_thread.finished.connect(self.worker_thread.deleteLater)
        self.worker_thread.start()
        return

    def on_calculation_result(self, result_data):
        # Update the UI with calculation results (runs in main thread)
        response = coerce_property_response(result_data)
        self.last_result_data = response
        result_data = response.to_legacy_dict()
        self.density_actual_kg_m3 = result_data.get("density_actual_kg_m3")
        self.density_normal_kg_m3 = result_data.get("density_normal_kg_m3")
        self.density_standard_kg_m3 = result_data.get("density_standard_kg_m3")
        self.set_thermo_warning_messages(result_data.get("warnings") or [])
        self.update_flow_conversion()
        self.ui.results_list.clear()
        self.ui.results_list.addItem(f"Calculating for {result_data['basis']} basis...")
        self.ui.results_list.addItem(f"Model / EOS: {result_data.get('model_display', self.calculator.eos)}")
        self.ui.results_list.addItem(f"Avg. Molecular Wt: {result_data['mw']:.2f} g/mol")
        # Density/Phase
        density_result = result_data['density_result']
        phase = result_data['phase']
        error = result_data['density_error']
        if error:
            self.ui.results_list.addItem("Phase @ selected conditions: N.A.")
            self.ui.results_list.addItem("Density @ selected conditions: N.A.")
            self.ui.results_list.addItem(f"Density/Phase Error: {error}")
        elif phase == "Two-Phase" and density_result is not None:
            if isinstance(density_result, tuple) and len(density_result) == 2:
                density_liq, density_gas = density_result
                self.ui.results_list.addItem(f"Phase @ selected conditions: {phase}")
                if density_liq is not None:
                    self.ui.results_list.addItem(f"  Density @ selected conditions (Liq): {density_liq:.3f} kg/m³")
                else:
                    self.ui.results_list.addItem("  Density @ selected conditions (Liq): N.A.")
                if density_gas is not None:
                    self.ui.results_list.addItem(f"  Density @ selected conditions (Vap): {density_gas:.3f} kg/m³")
                else:
                    self.ui.results_list.addItem("  Density @ selected conditions (Vap): N.A.")
            else:
                self.ui.results_list.addItem(f"Phase @ selected conditions: {phase} (Result format unexpected)")
        elif density_result is not None:
            self.ui.results_list.addItem(f"Phase @ selected conditions: {phase or 'Could not determine.'}")
            self.ui.results_list.addItem(f"Density @ selected conditions: {density_result:.3f} kg/m³")
        elif phase:
            self.ui.results_list.addItem(f"Phase @ selected conditions: {phase}")
            self.ui.results_list.addItem("Density @ selected conditions: N.A.")
        else:
            self.ui.results_list.addItem("Phase @ selected conditions: Could not determine.")
            self.ui.results_list.addItem("Density @ selected conditions: N.A.")

        def add_reference_density_line(label, density_value, phase_value, error_value):
            if density_value is not None:
                self.ui.results_list.addItem(f"{label}: {density_value:.3f} kg/m³")
                return
            detail = ""
            if error_value:
                detail = f" ({error_value})"
            elif phase_value == "Two-Phase":
                detail = " (Two-Phase)"
            self.ui.results_list.addItem(f"{label}: N.A.{detail}")

        add_reference_density_line(
            "Density @ normal conditions",
            result_data.get("density_normal_kg_m3"),
            result_data.get("density_normal_phase"),
            result_data.get("density_normal_error"),
        )
        add_reference_density_line(
            "Density @ standard conditions",
            result_data.get("density_standard_kg_m3"),
            result_data.get("density_standard_phase"),
            result_data.get("density_standard_error"),
        )
        # Bubble Point
        bubble_point = result_data['bubble_point']
        bp_error = result_data['bp_error']
        pressure_atm = result_data['pressure_atm']
        if bp_error:
            self.ui.results_list.addItem(f"Bubble Point Error: {bp_error}")
        elif bubble_point is not None:
            self.ui.results_list.addItem(f"Bubble Point @ {pressure_atm} atm: {bubble_point:.2f} °C")
        else:
            self.ui.results_list.addItem(f"Bubble Point: Calculation failed.")
        # LHV
        if self.lhv_data_loaded:
            mixture_lhv = result_data['mixture_lhv']
            missing_lhv = result_data['missing_lhv']
            if missing_lhv:
                self.ui.results_list.addItem(f"LHV Warning: No data for: {', '.join(missing_lhv)}")
            self.ui.results_list.addItem("-" * 40)
            lhv_display_values = build_lhv_display_values(mixture_lhv, result_data['mw'])
            volumetric_lines = [
                ("Mixture LHV = ", "MJ/Nm³"),
                ("= ", "kcal/Nm³"),
                ("= ", "MMkcal/Nm³"),
                ("= ", "GJ/Nm³"),
                ("= ", "MMBtu/Nm³"),
            ]
            for prefix, unit in volumetric_lines:
                self.ui.results_list.addItem(
                    f"{prefix}{format_lhv_display_value(lhv_display_values['volumetric'][unit])} {unit}"
                )

            self.ui.results_list.addItem("")
            self.ui.results_list.addItem("Mass basis:")
            for unit in ["MJ/kg", "MJ/t", "GJ/kg", "GJ/t", "kcal/kg", "kcal/t", "MMkcal/kg", "MMkcal/t", "MMBtu/kg", "MMBtu/t"]:
                self.ui.results_list.addItem(
                    f"= {format_lhv_display_value(lhv_display_values['mass_basis'][unit])} {unit}"
                )
        else:
            self.ui.results_list.addItem("-" * 40)
            self.ui.results_list.addItem("Mixture LHV = N/A (DB not loaded)")
        self.animate_progress_to(100, 500)

    def on_calculation_error(self, error_msg):
        self.invalidate_results(clear_visible_results=False)
        self.ui.results_list.clear()
        self.ui.results_list.addItem(error_msg)
        self.animate_progress_to(100, 500)

    def on_calculation_finished(self):
        self.ui.go_button.setEnabled(True)
        if hasattr(self.ui, "printResultsButton"):
            self.ui.printResultsButton.setEnabled(True)

    def normalize_composition(self):
        """Normalize the active composition column (mol% or wt%) so the sum is 100%."""
        # Determine active column
        active_col = 1 if self.ui.radioButton_mol_percent.isChecked() else 2
        row_count = self.ui.tableWidget.rowCount()
        total_row = row_count - 1

        # Gather values from the active column (excluding total row)
        values = []
        for row in range(total_row):
            item = self.ui.tableWidget.item(row, active_col)
            try:
                val = float(item.text()) if item and item.text().strip() else 0.0
            except ValueError:
                val = 0.0
            values.append(val)

        try:
            response = self.normalize_composition_use_case.execute(
                NormalizeCompositionRequest(percentages=tuple(values))
            )
        except ValueError:
            QMessageBox.warning(self, "Normalize", "Cannot normalize: total is zero.")
            return
        norm_values = response.percentages

        # Block signals to avoid recursion during update
        self.ui.tableWidget.blockSignals(True)
        # Set the normalized values in the table
        for row, norm_val in enumerate(norm_values):
            item = self.ui.tableWidget.item(row, active_col)
            if item:
                item.setText(f"{norm_val:.4f}")
        self.ui.tableWidget.blockSignals(False)

        self.invalidate_results()
        self.update_table_total()
        # Optionally, provide user feedback in the results list
        if hasattr(self.ui, 'results_list'):
            self.ui.results_list.addItem("Composition normalized to 100%.")

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


# --- Main Execution Block (Modified) ---
def main():
    app = QApplication(sys.argv)
    # It's good practice to set application name and version if distributing
    # app.setApplicationName("ThermoCalculator")
    # app.setApplicationVersion("1.1")

    app.setStyle(QStyleFactory.create("Windows"))
    app.setPalette(app.style().standardPalette())

    #app.setStyle(QStyleFactory.create("WindowsVista"))

    # Resolve and load reference data before creating the window.
    resource_locator = RuntimeResourceLocator()
    lhv_repository = SqliteLhvRepository(
        resource_locator.resolve("lhv_data.db")
    )
    lhv_data_for_app = lhv_repository.load_all()

    # Create and show the main window, passing LHV data
    window = MainWindow(lhv_data=lhv_data_for_app)
    window.show()

    sys.exit(app.exec())

if __name__ == "__main__":
    main()
