#!/usr/bin/env python3
"""Build the sequence-identity index the website uses to recognise an antibody.

Answers "have I seen this antibody before?" — an exact hit means an
experimental structure exists in SAbDab2 and can be looked up by PDB ID; a near
hit tells you how close the sequence is to anything solved; no hit suggests a
novel sequence.

The dataset carries PDB accessions, not therapeutic names, so the index reports
accessions and lets the user look the name up.

Usage:
    python scripts/build_identity_index.py --csv splits_final/ab_split.csv \\
        --out identity.json
"""

import argparse
import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", required=True)
    ap.add_argument("--out", default="identity.json")
    ap.add_argument("--max-pdbs", type=int, default=4,
                    help="accessions kept per unique sequence")
    args = ap.parse_args(argv)

    csv.field_size_limit(min(sys.maxsize, 2**31 - 1))
    heavy, light, pdb_meta = {}, {}, {}
    with open(args.csv, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            acc = row["PDB_ID"].replace("pdb_0000", "").upper()
            if row.get("Hseq"):
                heavy.setdefault(row["Hseq"], set()).add(acc)
            if row.get("Lseq"):
                light.setdefault(row["Lseq"], set()).add(acc)
            if acc not in pdb_meta:
                pdb_meta[acc] = {
                    "m": row.get("method") or "",
                    "r": row.get("resolution") or "",
                    "t": row.get("type") or "",
                    "ag": row.get("agtypes") or "",
                }

    def pack(table):
        return [[seq, sorted(accs)[:args.max_pdbs]] for seq, accs in table.items()]

    payload = {
        "_meta": {
            "source": "SAbDab2 ML dataset, DOI 10.5281/zenodo.23343734 (CC-BY-4.0)",
            "heavy_sequences": len(heavy),
            "light_sequences": len(light),
            "structures": len(pdb_meta),
            "note": "PDB accessions only; the dataset carries no therapeutic names.",
        },
        "H": pack(heavy),
        "L": pack(light),
        "pdb": pdb_meta,
    }
    out = Path(args.out)
    with open(out, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, separators=(",", ":"))
    print(f"wrote {out}  ({out.stat().st_size / 1024 / 1024:.2f} MB)")
    print(f"  heavy {len(heavy)}  light {len(light)}  structures {len(pdb_meta)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
