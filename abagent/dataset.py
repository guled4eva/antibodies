"""Reader for the SAbDab2 ML dataset split CSVs.

The split files (ab_split.csv, abag_split.csv, and their _sd variants) carry
ground-truth IMGT annotations for every structure -- sequences, CDR loops and
numbering -- which makes them usable for validation without parsing any .cif.

Dataset: Capel et al. 2026, CC-BY-4.0, DOI 10.5281/zenodo.23343734
Fetch it with scripts/get_data.py.
"""

import csv
import sys
from dataclasses import dataclass
from pathlib import Path

# Split files shipped with the dataset, most generally useful first.
SPLIT_FILES = ("ab_split.csv", "abag_split.csv", "ab_split_sd.csv", "abag_split_sd.csv")

# Columns this module relies on (documented in the dataset README).
REQUIRED = ("INSTANCE", "Hseq", "Lseq", "CDRH1", "CDRH2", "CDRH3",
            "CDRL1", "CDRL2", "CDRL3")


class DatasetError(RuntimeError):
    """Raised when the dataset is missing or not in the expected shape."""


@dataclass
class Record:
    """One antibody structure's sequences and ground-truth IMGT CDRs."""

    instance: str
    pdb_id: str
    heavy: str
    light: str
    cdrh1: str
    cdrh2: str
    cdrh3: str
    cdrl1: str
    cdrl2: str
    cdrl3: str
    split: str  # "train" / "test" as recorded in the file

    @property
    def is_paired(self) -> bool:
        return bool(self.heavy) and bool(self.light)


def find_split_csv(data_dir="data", name=None):
    """Locate a split CSV under data_dir (searched recursively)."""
    root = Path(data_dir)
    if not root.exists():
        raise DatasetError(
            f"No '{root}' directory. Download the dataset first:\n"
            f"    python scripts/get_data.py"
        )
    wanted = (name,) if name else SPLIT_FILES
    for filename in wanted:
        hits = sorted(root.rglob(filename))
        if hits:
            return hits[0]
    raise DatasetError(
        f"Found '{root}' but no split CSV inside it (looked for: "
        f"{', '.join(wanted)}). Is the archive extracted?"
    )


def iter_records(csv_path, limit=None, split=None, paired_only=False):
    """Yield Record objects from a split CSV.

    split: "train"/"test" to filter, or None for all.
    paired_only: skip single-domain entries that lack a chain.
    """
    csv.field_size_limit(min(sys.maxsize, 2**31 - 1))  # numbering lists are long
    path = Path(csv_path)
    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        missing = [c for c in REQUIRED if c not in (reader.fieldnames or [])]
        if missing:
            raise DatasetError(
                f"{path.name} is missing expected columns: {', '.join(missing)}.\n"
                f"Columns present: {', '.join((reader.fieldnames or [])[:12])}..."
            )
        split_col = "ab_split" if "ab_split" in reader.fieldnames else "ab_ag_split"
        count = 0
        for row in reader:
            rec = Record(
                instance=row.get("INSTANCE", ""),
                pdb_id=row.get("PDB_ID", ""),
                heavy=(row.get("Hseq") or "").strip(),
                light=(row.get("Lseq") or "").strip(),
                cdrh1=(row.get("CDRH1") or "").strip(),
                cdrh2=(row.get("CDRH2") or "").strip(),
                cdrh3=(row.get("CDRH3") or "").strip(),
                cdrl1=(row.get("CDRL1") or "").strip(),
                cdrl2=(row.get("CDRL2") or "").strip(),
                cdrl3=(row.get("CDRL3") or "").strip(),
                split=(row.get(split_col) or "").strip(),
            )
            if split and rec.split != split:
                continue
            if paired_only and not rec.is_paired:
                continue
            yield rec
            count += 1
            if limit and count >= limit:
                return


def describe(csv_path):
    """Quick summary of a split file: row count and column names."""
    path = Path(csv_path)
    csv.field_size_limit(min(sys.maxsize, 2**31 - 1))
    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        cols = reader.fieldnames or []
        rows = sum(1 for _ in reader)
    return {"file": str(path), "rows": rows, "n_columns": len(cols), "columns": cols}
