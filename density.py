import sys
import sqlite3
import os # To check if DB file exists

from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QComboBox, QListWidget,
    QPushButton, QLabel, QGridLayout, QVBoxLayout, QTableWidget,
    QTableWidgetItem, QMessageBox, QHeaderView, QHBoxLayout, QButtonGroup, QRadioButton, QStyleFactory
)
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QPalette, QColor

# --- Import the generated UI class ---
from gui import Ui_Dialog

# THERMO FLASH INTERFACE IMPORTS
# NOTE: Make sure PRMIX is imported if it's the only option for now
# If you add SRKMIX later, you'll need: from thermo.eos_mix import PRMIX, SRKMIX
from thermo import ChemicalConstantsPackage, CEOSGas, CEOSLiquid, FlashVL, PRMIX, FlashPureVLS
# from thermo.eos_mix import PRMIX # Explicit import if needed


# --- Constants and Data Loading (Keep as is) ---
MOLECULAR_WEIGHTS = {
    "        ": 0.0,  # Dummy entry
    "hydrogen": 2.016, "carbon dioxide": 44.01, "carbon monoxide": 28.01,
    "hydrogen sulfide": 34.08, "methane": 16.04, "ethane": 30.07,
    "ethylene": 28.05, "propane": 44.10, "propylene": 42.08,
    "isobutane": 58.12, "n-butane": 58.12, "isobutylene": 56.11,
    "n-butylene": 56.11, "1-butene": 56.11, "2-butene": 56.11,
    "butadiene": 54.09, "isopentane": 72.15, "n-pentane": 72.15,
    "hexane": 86.18, "heptane": 100.20, "benzene": 78.11,
    "toluene": 92.14, "o-xylene": 106.17, "m-xylene": 106.17,
    "p-xylene": 106.17, "ethylbenzene": 106.17, "cumene": 120.19,
    "octane": 114.23, "nonane": 128.26, "oxygen": 32.00,
    "nitrogen": 28.01,
}

R = 0.082057  # L·atm/(mol·K)

def load_lhv_data(db_path='lhv_data.db'):
    """Loads LHV data from the SQLite database and returns it."""
    loaded_data = {}
    if not os.path.exists(db_path):
        print(f"Warning: LHV Database file not found at {db_path}. LHV calculations will not be available.")
        return loaded_data

    conn = None
    try:
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        cursor.execute("SELECT component_name, lhv_mj_nm3 FROM lhv_values")
        rows = cursor.fetchall()
        for row in rows:
            loaded_data[row[0]] = row[1]
        print(f"Successfully loaded {len(loaded_data)} LHV entries from {db_path}.")
    except sqlite3.Error as e:
        print(f"Failed to load LHV data from database: {e}")
    finally:
        if conn:
            conn.close()
    return loaded_data

# --- MixtureCalculator Class (do not change this code) ---
class MixtureCalculator:
    """Handles the thermodynamic calculations, decoupled from the UI."""
    def __init__(self, components=None, eos='PRMIX'):
        self.components = components if components is not None else {}
        self.eos = eos
        self.constants = None
        self.properties = None

    def set_components(self, components):
        self.components = components
        self.constants = None # Reset cached constants when components change
        self.properties = None # Reset cached properties

    def set_eos(self, eos):
        # Add logic here if you support more EOS types later
        self.eos = eos
        print(f"EOS set to: {self.eos}") # Debug print

    def calculate_molecular_weight(self):
        if not self.components: return 0.0
        mw_sum = 0.0
        frac_sum = sum(self.components.values())
        if frac_sum == 0: return 0.0

        normalized_components = {comp: frac / frac_sum for comp, frac in self.components.items()}

        for comp_name, mole_frac in normalized_components.items():
            mw = MOLECULAR_WEIGHTS.get(comp_name, 0.0)
            mw_sum += mw * mole_frac
        # The average MW doesn't depend on the sum of fractions if normalized
        return mw_sum # / frac_sum if frac_sum > 0 else 0.0 -> removed as normalized


    def calculate_density(self, temperature_k, pressure_pa):
        """Calculates density, accounting for phase (liquid/vapor/two-phase).
           (Reverted to original phase iteration logic)"""
        if not self.components:
            # Return format matches original expectation: density, phase_str, error_str
            return None, None, "No components selected."

        total = sum(self.components.values())
        if total == 0:
            return None, None, "Total fraction is zero."

        zs = [x / total for x in self.components.values()]
        comp_list = list(self.components.keys())
        print(f"Calculating density (original logic) for: {comp_list} at T={temperature_k:.2f} K, P={pressure_pa:.1f} Pa with EOS={self.eos}")

        try:
            if self.constants is None or self.properties is None:
                print("Fetching chemical constants...")
                self.constants, self.properties = ChemicalConstantsPackage.from_IDs(comp_list)
                print("Constants fetched.")

            # Basic check for constants needed by EOS
            if None in self.constants.Tcs or None in self.constants.Pcs or None in self.constants.omegas:
                 missing_data_comps = [comp_list[i] for i, tc in enumerate(self.constants.Tcs) if tc is None]
                 return None, None, f"Missing critical properties for: {', '.join(missing_data_comps)}"

            eos_kwargs = {
                'Pcs': self.constants.Pcs,
                'Tcs': self.constants.Tcs,
                'omegas': self.constants.omegas,
            }

            eos_class = PRMIX # Assuming PRMIX only for now
            hcap_gases = self.properties.HeatCapacityGases if hasattr(self.properties, 'HeatCapacityGases') else None

            if len(comp_list) == 1:
                 # Original code used FlashPureVLS here, keep that part
                 gas = CEOSGas(eos_class, eos_kwargs=eos_kwargs, HeatCapacityGases=hcap_gases)
                 liquid = CEOSLiquid(eos_class, eos_kwargs=eos_kwargs, HeatCapacityGases=hcap_gases)
                 flasher = FlashPureVLS(self.constants, self.properties, gas=gas, liquids=[liquid], solids=[])
            else:
                 gas = CEOSGas(eos_class, HeatCapacityGases=hcap_gases, eos_kwargs=eos_kwargs, zs=zs)
                 liquid = CEOSLiquid(eos_class, HeatCapacityGases=hcap_gases, eos_kwargs=eos_kwargs, zs=zs)
                 flasher = FlashVL(self.constants, self.properties, liquid=liquid, gas=gas)

            print(f"Performing flash T={temperature_k}, P={pressure_pa}, zs={zs}")
            result = flasher.flash(T=temperature_k, P=pressure_pa, zs=zs)
            # Debug: Print result attributes relevant to original logic
            print(f"Flash result: Phases present = {len(result.phases)}, PurePhase='{getattr(result, 'phase', 'N/A')}'")
            if hasattr(result, 'phases'):
                 for i, p in enumerate(result.phases):
                      print(f"  Phase {i}: Liquid={p.is_liquid}, Gas={p.is_gas}, V={p.V()}")

            # Calculate Overall Molecular Weight (convert g/mol to kg/mol) - used by original logic
            MWs = self.constants.MWs if self.constants is not None else []
            if not MWs or None in MWs:
                return None, None, "Could not retrieve molecular weight for all components."
             # Use overall mole fractions (zs) for average MW calculation
            MW = sum([zi * Mwi for zi, Mwi in zip(zs, MWs)]) / 1000.0
            if MW <= 0:
                 return None, None, "Calculated average molecular weight is zero or negative."


            density_liq = None
            density_gas = None

            # --- Original Logic based on iterating result.phases (for mixtures) ---
            # --- and result.phase (for pure) ---
            if len(comp_list) == 1 and isinstance(flasher, FlashPureVLS):
                 # Handle pure component flash results based on 'phase' attribute
                 pure_phase_attr = getattr(result, 'phase', '').lower()
                 if pure_phase_attr in ['l', 'liquid']:
                      if result.liquid0 and result.liquid0.V() is not None and result.liquid0.V() > 0:
                           density_liq = MW / result.liquid0.V()
                 elif pure_phase_attr in ['g', 'v', 'gas', 'vapor']:
                      if result.gas and result.gas.V() is not None and result.gas.V() > 0:
                           density_gas = MW / result.gas.V()
                 elif '/' in pure_phase_attr: # Two phase like 'l/g'
                      if result.liquid0 and result.liquid0.V() is not None and result.liquid0.V() > 0:
                           density_liq = MW / result.liquid0.V()
                      if result.gas and result.gas.V() is not None and result.gas.V() > 0:
                           density_gas = MW / result.gas.V()
                 else:
                      # If phase attribute is unexpected, try checking phases list
                      if hasattr(result, 'phases'):
                           for phase in result.phases:
                                if phase.is_liquid and phase.V() is not None and phase.V() > 0:
                                    density_liq = MW / phase.V()
                                elif phase.is_gas and phase.V() is not None and phase.V() > 0:
                                    density_gas = MW / phase.V()
                      else:
                          return None, None, f"Unexpected pure phase result: {pure_phase_attr}"


            elif hasattr(result, 'phases'): # Mixture case: Iterate phases list
                 for phase in result.phases:
                      phase_vol = phase.V()
                      if phase_vol is not None and phase_vol > 0:
                           # Original logic used overall MW for phase density estimate
                           if phase.is_liquid:
                                density_liq = MW / phase_vol
                           elif phase.is_gas:
                                density_gas = MW / phase_vol
                      else:
                           print(f"Warning: Phase volume is None or zero for phase type (Liq={phase.is_liquid}, Gas={phase.is_gas})")

            else:
                 return None, None, "Flash result object did not contain expected 'phases' attribute."

            # --- Determine final return based on calculated densities ---
            if density_liq is not None and density_gas is not None:
                print(f"Determined Phase: Two-Phase, DensLiq={density_liq:.3f}, DensGas={density_gas:.3f}")
                return (density_liq, density_gas), "Two-Phase", None
            elif density_liq is not None:
                print(f"Determined Phase: Liquid, DensLiq={density_liq:.3f}")
                return density_liq, "Liquid", None
            elif density_gas is not None:
                print(f"Determined Phase: Vapor, DensGas={density_gas:.3f}")
                return density_gas, "Vapor", None
            else:
                # If no densities could be calculated
                print(f"Warning: Could not calculate density for any phase found.")
                # Try to determine phase based on VF if available, otherwise return Unknown
                phase_guess = "Unknown"
                if hasattr(result, 'VF'):
                     if result.VF == 0: phase_guess = "Liquid (calc failed)"
                     elif result.VF == 1: phase_guess = "Vapor (calc failed)"
                     elif 0 < result.VF < 1: phase_guess = "Two-Phase (calc failed)"
                return None, phase_guess, "Failed to calculate density values."

        except Exception as e:
            import traceback
            print(f"Error during density calculation: {e}\n{traceback.format_exc()}")
            # Return format matches original expectation
            return None, None, str(e)
    
    def calculate_bubble_point(self, pressure_pa):
        """Calculates the bubble point temperature."""
        if not self.components: return None, "No components selected."
        total = sum(self.components.values())
        if total == 0: return None, "Total fraction is zero."

        zs = [x / total for x in self.components.values()]
        comp_list = list(self.components.keys())
        print(f"Calculating bubble point for: {comp_list} at P={pressure_pa:.1f} Pa")

        try:
            if self.constants is None or self.properties is None:
                 self.constants, self.properties = ChemicalConstantsPackage.from_IDs(comp_list)

            if None in self.constants.Tcs or None in self.constants.Pcs or None in self.constants.omegas:
                 missing_data_comps = [comp_list[i] for i, tc in enumerate(self.constants.Tcs) if tc is None]
                 return None, f"Missing critical properties for: {', '.join(missing_data_comps)}"

            eos_kwargs = {'Pcs': self.constants.Pcs, 'Tcs': self.constants.Tcs, 'omegas': self.constants.omegas,}
            hcap_gases = self.properties.HeatCapacityGases if hasattr(self.properties, 'HeatCapacityGases') else None
            eos_class = PRMIX # Assuming PRMIX for now

            if len(comp_list) == 1:
                gas = CEOSGas(eos_class, eos_kwargs=eos_kwargs, HeatCapacityGases=hcap_gases)
                liquid = CEOSLiquid(eos_class, eos_kwargs=eos_kwargs, HeatCapacityGases=hcap_gases)
                flasher = FlashPureVLS(self.constants, self.properties, gas=gas, liquids=[liquid], solids=[])
            else:
                gas = CEOSGas(eos_class, HeatCapacityGases=hcap_gases, eos_kwargs=eos_kwargs, zs=zs)
                liquid = CEOSLiquid(eos_class, HeatCapacityGases=hcap_gases, eos_kwargs=eos_kwargs, zs=zs)
                flasher = FlashVL(self.constants, self.properties, liquid=liquid, gas=gas)

            result = flasher.flash(P=pressure_pa, VF=0, zs=zs) # VF=0 for bubble point
            print(f"Bubble point flash result T={result.T}")
            return result.T - 273.15, None

        except Exception as e:
             import traceback
             print(f"Error during bubble point calculation: {e}\n{traceback.format_exc()}")
             return None, str(e)

    def calculate_lhv(self, lhv_data):
        """Calculates the LHV of the mixture based on mole fractions (MJ/Nm³)."""
        if not self.components: return 0.0, []
        mixture_lhv = 0.0
        missing_components = []
        total_fraction = sum(self.components.values())
        if total_fraction == 0: return 0.0, []

        normalized_components = {comp: frac / total_fraction for comp, frac in self.components.items()}
        for comp_name, mole_frac in normalized_components.items():
            component_lhv = lhv_data.get(comp_name)
            if component_lhv is not None:
                mixture_lhv += mole_frac * component_lhv
            else:
                if comp_name.strip(): missing_components.append(comp_name)
        return mixture_lhv, missing_components


# --- MainWindow Class (Modified for gui.py) ---
class MainWindow(QMainWindow):
    def __init__(self, lhv_data):
        super().__init__()

        # --- Set up the UI from Designer ---
        self.ui = Ui_Dialog()
        self.ui.setupUi(self) # Setup the UI onto this QMainWindow

        # --- Store LHV data and calculator instance ---
        self.lhv_database = lhv_data
        self.lhv_data_loaded = bool(self.lhv_database)
        self.calculator = MixtureCalculator()

        if not self.lhv_data_loaded:
             print("Warning: MainWindow initialized with no LHV data.")
             # You could disable LHV display or show a warning label here if needed

        # --- Populate widgets ---
        self.populate_comboboxes()
        self.setup_table()

        # --- Connect signals ---
        self.connect_signals()

        # --- Initialize UI State ---
        self.ui.radioButton_mol_percent.setChecked(True) # Default Mol%
        self.on_input_basis_changed() # Set initial table state

        # Set Window Title (Optional - can also be set in Designer)
        self.setWindowTitle("Thermo Calculator")

    def populate_comboboxes(self):
        """Populate the ComboBox widgets with options."""
        # Component Selection ComboBox
        self.ui.comboBox_select_components.addItems(list(MOLECULAR_WEIGHTS.keys()))

        # Temperature ComboBox
        self.ui.comboBox_select_temperature.addItems(["0 °C", "5 °C", "10 °C", "15 °C", "25 °C", "50 °C", "100 °C", "150 °C", "200 °C"])
        self.ui.comboBox_select_temperature.setCurrentText("15 °C") # Default

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
        selected_eos = self.ui.comboBox_select_EOS.currentText()
        self.calculator.set_eos(selected_eos)

    def add_component(self):
        """Add the selected component to the list and table."""
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

        # Clear results_list after successful add
        self.clear_results_list()

    def remove_component(self):
        """Remove the selected component from the list and table."""
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

        # Clear results_list after successful remove
        self.clear_results_list()

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
        # Clear results_list at the start
        self.clear_results_list()
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
        self.ui.comboBox_select_temperature.setCurrentText("15 °C")
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
                inactive_item.setText("") # Clear inactive column value
                inactive_item.setFlags(Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable) # Not editable
            else: # Total row
                active_item.setFlags(Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable) # Not editable
                inactive_item.setFlags(Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable) # Not editable

        # Unblock signals
        self.ui.tableWidget.blockSignals(False)
        self.update_table_total() # Update totals after changing basis

    def on_table_item_changed(self, item):
        """Recalculate total when a percentage value changes."""
        # Only update if the changed item is in an active percentage column and not the total row
        active_col = 1 if self.ui.radioButton_mol_percent.isChecked() else 2
        if item.column() == active_col and item.row() != self.ui.tableWidget.rowCount() - 1:
             # Basic validation: try converting to float, set red background if invalid?
             try:
                 float(item.text())
                 item.setBackground(QColor("white")) # Reset background on valid input
             except ValueError:
                 if item.text().strip(): # Don't color empty cells red
                      item.setBackground(QColor("pink")) # Indicate invalid number format
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
        self.ui.results_list.clear() # Use the new name from gui.py

        # --- Get components and compositions from tableWidget ---
        comp_names = []
        mol_percents = [] # Store as percentages first
        wt_percents = []
        row_count = self.ui.tableWidget.rowCount() - 1 # Exclude total row
        valid_input = True

        for r in range(row_count):
            comp_item = self.ui.tableWidget.item(r, 0)
            mol_item = self.ui.tableWidget.item(r, 1)
            wt_item = self.ui.tableWidget.item(r, 2)

            if not comp_item: continue # Should not happen if added correctly

            comp_names.append(comp_item.text())
            try:
                # Get Mol % (even if column is inactive, read stored value or default to 0)
                mol_val = float(mol_item.text()) if mol_item and mol_item.text().strip() else 0.0
                mol_percents.append(mol_val)
            except ValueError:
                 mol_percents.append(0.0)
                 if self.ui.radioButton_mol_percent.isChecked(): valid_input = False

            try:
                 # Get Wt %
                 wt_val = float(wt_item.text()) if wt_item and wt_item.text().strip() else 0.0
                 wt_percents.append(wt_val)
            except ValueError:
                 wt_percents.append(0.0)
                 if self.ui.radioButton_wt_percent.isChecked(): valid_input = False

        if not valid_input:
             self.ui.results_list.addItem("Error: Invalid numeric input in composition table.")
             return
        if not comp_names:
             self.ui.results_list.addItem("Error: No components selected.")
             return

        # --- Convert to mole fractions based on input basis ---
        mole_fracs_dict = {} # Final {name: mole_fraction}
        if self.ui.radioButton_mol_percent.isChecked():
            basis = "Mol %"
            total_mol_percent = sum(mol_percents)
            if abs(total_mol_percent - 100.0) > 1e-4 : # Check sum is 100%
                self.ui.results_list.addItem("Error: Total Mol % does not sum to 100.")
                return
            if total_mol_percent == 0:
                 self.ui.results_list.addItem("Error: Total Mol % is zero.")
                 return
            # Convert valid % to fractions
            for name, percent in zip(comp_names, mol_percents):
                 mole_fracs_dict[name] = percent / 100.0

        else: # Wt % basis
            basis = "Wt %"
            total_wt_percent = sum(wt_percents)
            if abs(total_wt_percent - 100.0) > 1e-4:
                self.ui.results_list.addItem("Error: Total Wt % does not sum to 100.")
                return
            if total_wt_percent == 0:
                 self.ui.results_list.addItem("Error: Total Wt % is zero.")
                 return

            # Convert Wt % to mole fractions
            moles = {}
            total_moles = 0.0
            conversion_possible = True
            for name, percent in zip(comp_names, wt_percents):
                mw = MOLECULAR_WEIGHTS.get(name)
                if mw is None or mw <= 0:
                    self.ui.results_list.addItem(f"Error: Missing or invalid MW for {name}.")
                    conversion_possible = False
                    break
                # Calculate moles based on 100g total mass assumption
                moles[name] = (percent / 100.0) * 100.0 / mw # (wt_frac * total_mass) / mw
                total_moles += moles[name]

            if not conversion_possible: return
            if total_moles == 0:
                self.ui.results_list.addItem("Error: Total moles is zero after Wt% conversion.")
                return

            for name in comp_names:
                mole_fracs_dict[name] = moles[name] / total_moles

        # --- Set components in calculator ---
        self.calculator.set_components(mole_fracs_dict)

        # --- Get Temperature and Pressure ---
        try:
            temp_text = self.ui.comboBox_select_temperature.currentText()
            T_c = float(temp_text.split("°")[0])
            T_k = T_c + 273.15
        except ValueError:
            self.ui.results_list.addItem("Error: Invalid Temperature selection.")
            return # Or use default T_k = 273.15 + 15.0

        try:
            pressure_text = self.ui.comboBox_select_pressure.currentText()
            pressure_atm = float(pressure_text.split(" ")[0])
            pressure_pa = pressure_atm * 101325.0
            if pressure_pa <=0: raise ValueError("Pressure must be positive")
        except ValueError:
            self.ui.results_list.addItem("Error: Invalid Pressure selection.")
            return # Or use default pressure_pa = 101325.0

        # --- Perform Calculations ---
        self.ui.results_list.addItem(f"Calculating for {basis} basis...")

        mw = self.calculator.calculate_molecular_weight()
        self.ui.results_list.addItem(f"Avg. Molecular Wt: {mw:.2f} g/mol")



        # --- Density Calculation Call ---
        # Now expects return: density_result, phase_str, error_str
        density_result, phase, error = self.calculator.calculate_density(T_k, pressure_pa)

        # --- Adjusted Result Display ---
        if error:
            self.ui.results_list.addItem(f"Density/Phase Error: {error}")
        # Check phase string AND if density_result is not None (important!)
        elif phase == "Two-Phase" and density_result is not None:
             # density_result should be a tuple (liq, gas)
             if isinstance(density_result, tuple) and len(density_result) == 2:
                  density_liq, density_gas = density_result
                  self.ui.results_list.addItem(f"Phase @(T,P): {phase}") # Display phase
                  # Only display densities if they are valid numbers
                  if density_liq is not None:
                      self.ui.results_list.addItem(f"  Density (Liq): {density_liq:.3f} kg/m³")
                  else:
                      self.ui.results_list.addItem(f"  Density (Liq): Calculation Failed")
                  if density_gas is not None:
                      self.ui.results_list.addItem(f"  Density (Vap): {density_gas:.3f} kg/m³")
                  else:
                       self.ui.results_list.addItem(f"  Density (Vap): Calculation Failed")
             else:
                  # Should not happen if calculate_density returns correctly
                  self.ui.results_list.addItem(f"Phase: {phase} (Result format unexpected)")

        elif phase == "Liquid" and density_result is not None:
             self.ui.results_list.addItem(f"Phase @(T,P): {phase}")
             self.ui.results_list.addItem(f"  Density: {density_result:.3f} kg/m³")
        elif phase == "Vapor" and density_result is not None:
             self.ui.results_list.addItem(f"Phase @(T,P): {phase}")
             self.ui.results_list.addItem(f"  Density: {density_result:.3f} kg/m³")
        elif phase: # Phase determined but density calculation failed
             self.ui.results_list.addItem(f"Phase @(T,P): {phase}")
             self.ui.results_list.addItem(f"  Density: Calculation Failed")
        else: # No phase determined and likely an error occurred before returning
             self.ui.results_list.addItem(f"Density/Phase: Could not determine.")

        # Bubble Point Calculation
        bubble_point, bp_error = self.calculator.calculate_bubble_point(pressure_pa)
        if bp_error:
            self.ui.results_list.addItem(f"Bubble Point Error: {bp_error}")
        elif bubble_point is not None:
            self.ui.results_list.addItem(f"Bubble Point @ {pressure_atm} atm: {bubble_point:.2f} °C")
        else:
            self.ui.results_list.addItem(f"Bubble Point: Calculation failed.")


        # LHV Calculation
        if self.lhv_data_loaded:
            mixture_lhv, missing_lhv = self.calculator.calculate_lhv(self.lhv_database)
            if missing_lhv:
                self.ui.results_list.addItem(f"LHV Warning: No data for: {', '.join(missing_lhv)}")
            self.ui.results_list.addItem("-" * 40) # Separator
            self.ui.results_list.addItem(f"Mixture LHV = {mixture_lhv:.2f} MJ/Nm³")
        else:
            self.ui.results_list.addItem("-" * 40)
            self.ui.results_list.addItem("Mixture LHV = N/A (DB not loaded)")


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

        total = sum(values)
        if total == 0:
            QMessageBox.warning(self, "Normalize", "Cannot normalize: total is zero.")
            return

        # Block signals to avoid recursion during update
        self.ui.tableWidget.blockSignals(True)
        norm_values = []
        # Normalize and round all but the last value
        for row, val in enumerate(values):
            if row < total_row - 1:
                norm_val = round((val / total) * 100, 4)
                norm_values.append(norm_val)
            else:
                # Last value: set to 100 - sum of others, rounded to 4 decimals
                last_val = round(100.0 - sum(norm_values), 4)
                norm_values.append(last_val)
        # Set the normalized values in the table
        for row, norm_val in enumerate(norm_values):
            item = self.ui.tableWidget.item(row, active_col)
            if item:
                item.setText(f"{norm_val:.4f}")
        self.ui.tableWidget.blockSignals(False)

        self.update_table_total()
        # Optionally, provide user feedback in the results list
        if hasattr(self.ui, 'results_list'):
            self.ui.results_list.addItem("Composition normalized to 100%.")

    def on_temperature_or_pressure_changed(self):
        """Clear the results_list when temperature or pressure is changed."""
        self.ui.results_list.clear()

    def clear_results_list(self):
        """Clear the results_list widget."""
        self.ui.results_list.clear()


# --- Main Execution Block (Modified) ---
def main():
    app = QApplication(sys.argv)
    # It's good practice to set application name and version if distributing
    # app.setApplicationName("ThermoCalculator")
    # app.setApplicationVersion("1.1")

    app.setStyle(QStyleFactory.create("Windows"))
    app.setPalette(app.style().standardPalette())

    #app.setStyle(QStyleFactory.create("WindowsVista"))

    # Load the LHV data before creating the window
    lhv_data_for_app = load_lhv_data()

    # Create and show the main window, passing LHV data
    window = MainWindow(lhv_data=lhv_data_for_app)
    window.show()

    sys.exit(app.exec())

if __name__ == "__main__":
    main()