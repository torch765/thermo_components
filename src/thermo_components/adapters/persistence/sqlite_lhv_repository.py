"""SQLite implementation of the LHV repository port."""

from collections.abc import Mapping
from pathlib import Path
import sqlite3


class SqliteLhvRepository:
    """Store and retrieve component LHV values from SQLite."""

    TABLE_NAME = "lhv_values"

    def __init__(self, database_path: str | Path):
        self.database_path = Path(database_path)

    def load_all(self) -> dict[str, float]:
        if not self.database_path.exists():
            print(
                "Warning: LHV Database file not found at "
                f"{self.database_path}. LHV calculations will not be available."
            )
            return {}

        try:
            with sqlite3.connect(self.database_path) as connection:
                rows = connection.execute(
                    "SELECT component_name, lhv_mj_nm3 FROM lhv_values"
                ).fetchall()
        except sqlite3.Error as exc:
            print(f"Failed to load LHV data from database: {exc}")
            return {}

        values = {str(name): float(lhv) for name, lhv in rows}
        print(
            f"Successfully loaded {len(values)} LHV entries from "
            f"{self.database_path}."
        )
        return values

    def upsert_all(self, values: Mapping[str, float]) -> bool:
        """Create the schema and upsert a normalized set of LHV values."""
        normalized_values = [
            (str(name).strip().lower(), float(lhv))
            for name, lhv in values.items()
            if str(name).strip()
        ]

        try:
            with sqlite3.connect(self.database_path) as connection:
                connection.execute(
                    """
                    CREATE TABLE IF NOT EXISTS lhv_values (
                        component_name TEXT PRIMARY KEY,
                        lhv_mj_nm3 REAL NOT NULL
                    )
                    """
                )
                connection.executemany(
                    """
                    INSERT OR REPLACE INTO lhv_values (
                        component_name,
                        lhv_mj_nm3
                    )
                    VALUES (?, ?)
                    """,
                    normalized_values,
                )
        except sqlite3.Error as exc:
            print(f"Database error: {exc}")
            return False

        return True
