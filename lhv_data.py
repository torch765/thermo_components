from thermo_components.adapters.persistence import SqliteLhvRepository

# Data provided by the user (ensure keys match MOLECULAR_WEIGHTS keys)
LHV_DATA_RAW = {
    "hydrogen": 10.8,
    "carbon dioxide": 0.0,
    "carbon monoxide": 12.6,
    "hydrogen sulfide": 26.3, # Corrected value based on common sources might be ~23.3, but using user value
    "methane": 35.8,
    "ethane": 65.0,
    "ethylene": 59.0,
    "propane": 93.0,
    "propylene": 85.8,
    "isobutane": 119.6,
    "n-butane": 119.6,
    "isobutylene": 111.4,
    "n-butylene": 111.4, # Note: Often n-butylene refers to 1-butene or 2-butene
    "1-butene": 111.4,
    "2-butene": 111.4, # Assuming cis/trans avg is same as 1-butene for this data
    "butadiene": 111.1, # Often higher (~117), using user value
    "isopentane": 148.0,
    "n-pentane": 148.0,
    "hexane": 179.0, # Represents n-hexane typically
    "heptane": 210.0, # Represents n-heptane typically
    "benzene": 165.0, # Often lower (~140), using user value
    "toluene": 185.0, # Often lower (~167), using user value
    "o-xylene": 204.0, # Often lower (~192), using user value
    "m-xylene": 204.0, # Often lower (~192), using user value
    "p-xylene": 204.0, # Often lower (~192), using user value
    "ethylbenzene": 204.0, # Often lower (~195), using user value
    "cumene": 220.0, # Often lower (~216), using user value
    "octane": 232.0, # Represents n-octane typically
    "nonane": 263.0, # Represents n-nonane typically
    "oxygen": 0.0,
    "nitrogen": 0.0,
    # Inerts / non-fuels (set to 0.0 to avoid “missing LHV” warnings)
    "argon": 0.0,
    "water": 0.0,
    # Ensure any other components in MOLECULAR_WEIGHTS are added if needed
    # Add the dummy entry if you want it stored, though it's not needed for calculation
    # "        ": 0.0
}

DB_NAME = 'lhv_data.db'

def create_database():
    repository = SqliteLhvRepository(DB_NAME)
    if repository.upsert_all(LHV_DATA_RAW):
        print(f"Database '{DB_NAME}' created and populated successfully.")
        for component_name, lhv in repository.load_all().items():
            print(f"Component: {component_name}, LHV: {lhv}")

if __name__ == "__main__":
    create_database()
