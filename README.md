# Thermo Components

   A PyQt6-based desktop application for calculating thermodynamic properties of gas and liquid mixtures.  
   The program allows users to select chemical components, specify their composition, and compute properties such as density, phase, bubble point temperature, and lower heating value (LHV) using the `thermo` Python library.

   ## Features
   - User-friendly GUI
   - Supports multiple components 
   - Calculates density at selected conditions, normal conditions, and standard conditions
   - Calculates phase, bubble point, and LHV

   ## Requirements
   - Python 3.x
   - PyQt6
   - thermo
   - openpyxl (required for Excel report export)

   ## Usage
   1. Run `python density.py`
   2. Select components and enter compositions
   3. Set temperature and pressure
   4. Enter composition in mol% or wt%
   5. Click Go to see results for:
      - density at the selected temperature and pressure
      - density at normal conditions for `Nm3` basis support
      - density at standard conditions for `Sm3` / standard liquid basis support
      - phase, bubble point, and LHV

   ## Reference Density Bases
   - Normal conditions = `0 °C`, `1 atm`
   - Standard conditions = `60 °F`, `1 atm` (`15.5555556 °C`, `1 atm`)

   ## Flow Tab
   The `Flow` tab converts between mass-flow units and reference-volume-flow units using the densities produced by the first tab.

   Supported units:
   - Mass: `kg/h`, `kg/d`, `t/h`, `t/d`, `lb/h`, `lb/d`, `Klb/h`, `Klb/d`
   - Reference volume: `Nm3/h`, `Nm3/d`, `Sm3/h`, `Sm3/d`, `bbl/h`, `bbl/d`, `SCFH`, `MSCFD`, `MMSCFD`

   Basis logic:
   - `Nm3` uses the normal reference density at `0 °C` and `1 atm`
   - `Sm3`, `SCFH`, `MSCFD`, `MMSCFD`, and `bbl` use the standard reference density at `60 °F` and `1 atm`
   - Cross-basis conversions such as `Nm3` to `Sm3` convert through mass and may require both reference densities

   Density usage:
   - Mass-to-mass and same-basis reference-volume conversions do not require density
   - Mass-to-`Nm3` and `Nm3`-to-mass conversions use normal density
   - Mass-to-standard-basis volume and standard-basis-volume-to-mass conversions use standard density

   ## License
   GPL-3.0
