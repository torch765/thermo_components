# LHV Data Basis

## Stored Unit And Convention

`lhv_data.db` stores lower heating values in `MJ/Nm³`, where the normal
volume basis is `0 °C` and `1 atm`.

For the components added in June 2026, values were derived from NIST Chemistry
WebBook standard gas-phase formation enthalpies at 298.15 K. Complete
combustion products were taken as:

- carbon dioxide gas
- water vapor, so the result is LHV rather than HHV
- sulfur dioxide gas for sulfur-containing fuels
- nitrogen gas for ammonia

The molar combustion result was divided by the project normal molar volume of
`22.414 Nm³/kmol` and rounded to one decimal place for consistency with the
existing database.

Product formation enthalpies used:

- [`CO₂(g)`](https://webbook.nist.gov/cgi/cbook.cgi?ID=C124389&Mask=1):
  `-393.51 kJ/mol`
- [`H₂O(g)`](https://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185&Mask=1):
  `-241.826 kJ/mol`
- [`SO₂(g)`](https://webbook.nist.gov/cgi/cbook.cgi?ID=C7446095&Mask=1):
  `-296.81 kJ/mol`

## Added Component Values

| Component | NIST gas formation enthalpy (kJ/mol) | LHV (MJ/Nm³) | NIST source |
| --- | ---: | ---: | --- |
| ammonia | -45.94 | 14.1 | [NIST](https://webbook.nist.gov/cgi/cbook.cgi?ID=C7664417&Mask=1) |
| carbonyl sulfide | -138.41 | 24.6 | [NIST](https://webbook.nist.gov/cgi/cbook.cgi?ID=C463581&Mask=1) |
| sulfur dioxide | -296.81 | 0.0 | [NIST](https://webbook.nist.gov/cgi/cbook.cgi?ID=C7446095&Mask=1) |
| carbon disulfide | 116.94 | 49.3 | [NIST](https://webbook.nist.gov/cgi/cbook.cgi?ID=C75150&Mask=1) |
| methyl mercaptan | -22.8 | 51.4 | [NIST](https://webbook.nist.gov/cgi/cbook.cgi?ID=C74931&Mask=1) |
| methanol | -201.49 | 30.1 | [NIST](https://webbook.nist.gov/cgi/cbook.cgi?ID=C67561&Mask=1) |
| MTBE | -283.2 | 139.9 | [NIST](https://webbook.nist.gov/cgi/cbook.cgi?ID=C1634044&Mask=1) |
| dimethyl ether | -184.1 | 59.3 | [NIST](https://webbook.nist.gov/cgi/cbook.cgi?ID=C115106&Mask=1) |
| decane | -249.7 | 283.1 | [NIST](https://webbook.nist.gov/cgi/cbook.cgi?ID=C124185&Mask=1) |
| dodecane | -290.9 | 338.0 | [NIST](https://webbook.nist.gov/cgi/cbook.cgi?ID=C112403&Mask=1) |
| cyclohexane | -123.1 | 164.6 | [NIST](https://webbook.nist.gov/cgi/cbook.cgi?ID=C110827&Mask=1) |
| methylcyclohexane | -154.8 | 191.5 | [NIST](https://webbook.nist.gov/cgi/cbook.cgi?ID=C108872&Mask=1) |
| cis-2-butene | -7.7 | 113.0 | [NIST](https://webbook.nist.gov/cgi/cbook.cgi?ID=C590181&Mask=1) |
| trans-2-butene | -10.8 | 112.9 | [NIST](https://webbook.nist.gov/cgi/cbook.cgi?ID=C624646&Mask=1) |

Sulfur dioxide is explicitly stored as zero because it is already the selected
fully oxidized sulfur product under this convention.
