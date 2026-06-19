from pathlib import Path

from thermo_components.adapters.persistence import SqliteLhvRepository


# LHV values in MJ/Nm³ at 0 °C and 1 atm. New values added in June 2026 are
# derived from NIST gas-phase standard formation enthalpies with water vapor
# as the combustion product. See docs/LHV_DATA.md.
LHV_DATA_RAW = {
    "1-butene": 111.4,
    "ammonia": 14.1,
    "argon": 0.0,
    "benzene": 165.0,
    "butadiene": 111.1,
    "carbon dioxide": 0.0,
    "carbon disulfide": 49.3,
    "carbon monoxide": 12.6,
    "carbonyl sulfide": 24.6,
    "cis-2-butene": 113.0,
    "cumene": 220.0,
    "cyclohexane": 164.6,
    "decane": 283.1,
    "dimethyl ether": 59.3,
    "dodecane": 338.0,
    "ethane": 65.0,
    "ethylbenzene": 204.0,
    "ethylene": 59.0,
    "heptane": 210.0,
    "hexane": 179.0,
    "hydrogen": 10.8,
    "hydrogen sulfide": 26.3,
    "isobutane": 119.6,
    "isobutylene": 111.4,
    "isopentane": 148.0,
    "m-xylene": 204.0,
    "methane": 35.8,
    "methanol": 30.1,
    "methyl mercaptan": 51.4,
    "methylcyclohexane": 191.5,
    "MTBE": 139.9,
    "n-butane": 119.6,
    "n-pentane": 148.0,
    "nitrogen": 0.0,
    "nonane": 263.0,
    "o-xylene": 204.0,
    "octane": 232.0,
    "oxygen": 0.0,
    "p-xylene": 204.0,
    "propane": 93.0,
    "propylene": 85.8,
    "sulfur dioxide": 0.0,
    "toluene": 185.0,
    "trans-2-butene": 112.9,
    "water": 0.0,
}

DB_NAME = "lhv_data.db"


def create_database():
    database_path = Path(DB_NAME)
    database_path.unlink(missing_ok=True)
    repository = SqliteLhvRepository(database_path)
    if repository.upsert_all(LHV_DATA_RAW):
        print(f"Database '{DB_NAME}' created and populated successfully.")
        for component_name, lhv in repository.load_all().items():
            print(f"Component: {component_name}, LHV: {lhv}")


if __name__ == "__main__":
    create_database()
