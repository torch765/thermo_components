import sqlite3

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
    # Ensure any other components in MOLECULAR_WEIGHTS are added if needed
    # Add the dummy entry if you want it stored, though it's not needed for calculation
    # "        ": 0.0
}

DB_NAME = 'lhv_data.db'
TABLE_NAME = 'lhv_values'

def create_database():
    conn = None
    try:
        conn = sqlite3.connect(DB_NAME)
        cursor = conn.cursor()

        # Create table - Component name is primary key
        cursor.execute(f'''
        CREATE TABLE IF NOT EXISTS {TABLE_NAME} (
            component_name TEXT PRIMARY KEY,
            lhv_mj_nm3 REAL NOT NULL
        )
        ''')

        # Insert data - Use INSERT OR IGNORE to avoid errors if run multiple times
        # or INSERT OR REPLACE to update values if they change
        for name, lhv in LHV_DATA_RAW.items():
             # Ensure name matches the keys in your main script's MOLECULAR_WEIGHTS exactly
            component_key = name.strip().lower() # Normalize name just in case
            if component_key: # Avoid inserting blank keys
                cursor.execute(f'''
                INSERT OR REPLACE INTO {TABLE_NAME} (component_name, lhv_mj_nm3)
                VALUES (?, ?)
                ''', (component_key, lhv))

        conn.commit()
        print(f"Database '{DB_NAME}' created and populated successfully.")

    except sqlite3.Error as e:
        print(f"Database error: {e}")
    except Exception as e:
        print(f"An error occurred: {e}")
    finally:
        if conn:
            conn.close()

            def print_database():
                conn = None
                try:
                    conn = sqlite3.connect(DB_NAME)
                    cursor = conn.cursor()
                    cursor.execute(f"SELECT * FROM {TABLE_NAME}")
                    rows = cursor.fetchall()
                    for row in rows:
                        print(f"Component: {row[0]}, LHV: {row[1]}")
                except sqlite3.Error as e:
                    print(f"Database error: {e}")
                except Exception as e:
                    print(f"An error occurred: {e}")
                finally:
                    if conn:
                        conn.close()

            print_database()

if __name__ == "__main__":
    create_database()