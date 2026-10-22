#!/usr/bin/env python3
"""Build a per-IMGT-position solvent-exposure table from SAbDab2 structures.

Why: a sequence liability only matters if the residue is on the surface, but a
user pasting a sequence has no structure. Averaging exposure per IMGT position
across many experimental structures gives a lookup table that answers "is this
position usually buried?" instantly, with no structure prediction.

Method notes that matter:
  * One representative per `ab_cluster`, not every structure. 16,415 structures
    come from only 2,656 clusters, so averaging over all of them would
    over-weight antibodies that happen to have been crystallised many times.
  * Train split only by default, so the table can be validated on test later.
  * Antigen chains are EXCLUDED from the SASA calculation. Developability
    (oxidation, deamidation, aggregation) happens in the vial, where no antigen
    is bound, so the paratope should count as exposed.
  * The .cif files are already IMGT-numbered (verified: auth_seq_id + insertion
    code match the CSV's VH_numbering_list exactly), so no remapping is needed.

Usage:
    python scripts/build_exposure_table.py --cif-dir splits_final \\
        --csv splits_final/ab_split.csv --out abagent/data/exposure.json
"""

import argparse
import csv
import json
import sys
import time
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from abagent.sasa import parse_cif_atoms, residue_sasa  # noqa: E402


def find_cif(cif_dir, name):
    direct = Path(cif_dir) / name
    if direct.exists():
        return direct
    hits = list(Path(cif_dir).rglob(name))
    return hits[0] if hits else None


def representatives(csv_path, split, cluster_col="ab_cluster"):
    """One row per antibody cluster, so common antibodies don't dominate."""
    csv.field_size_limit(min(sys.maxsize, 2**31 - 1))
    seen = set()
    out = []
    with open(csv_path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            if split and row.get("ab_split") != split:
                continue
            key = row.get(cluster_col) or row.get("INSTANCE")
            if key in seen:
                continue
            seen.add(key)
            out.append(row)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", required=True)
    ap.add_argument("--cif-dir", required=True)
    ap.add_argument("--out", default="abagent/data/exposure.json")
    ap.add_argument("--split", default="train", choices=["train", "test", "all"])
    ap.add_argument("--limit", type=int, default=0, help="0 = all representatives")
    ap.add_argument("--n-points", type=int, default=100)
    args = ap.parse_args(argv)

    split = None if args.split == "all" else args.split
    rows = representatives(args.csv, split)
    if args.limit:
        rows = rows[:args.limit]
    print(f"{len(rows)} cluster representatives ({args.split} split)")

    acc = {"H": defaultdict(list), "L": defaultdict(list)}
    done = missing = failed = 0
    start = time.time()

    for i, row in enumerate(rows, 1):
        path = find_cif(args.cif_dir, row["cif_file"])
        if path is None:
            missing += 1
            continue
        chain_of = {}
        if row.get("Hchain") and row["Hchain"] != "+":
            chain_of[row["Hchain"]] = "H"
        if row.get("Lchain") and row["Lchain"] != "+":
            chain_of[row["Lchain"]] = "L"
        if not chain_of:
            continue
        try:
            atoms = [a for a in parse_cif_atoms(path) if a.chain in chain_of]
            if not atoms:
                failed += 1
                continue
            for res in residue_sasa(atoms, n_points=args.n_points):
                ctype = chain_of[res["chain"]]
                pos = f"{res['resseq']}{res['icode']}".strip()
                if res["rsa"] == res["rsa"]:  # not NaN
                    acc[ctype][pos].append(res["rsa"])
            done += 1
        except Exception as exc:
            failed += 1
            if failed <= 3:
                print(f"  ! {path.name}: {exc}")
        if i % 50 == 0 or i == len(rows):
            rate = i / max(time.time() - start, 1e-9)
            eta = (len(rows) - i) / max(rate, 1e-9)
            print(f"  {i}/{len(rows)}  ok={done} missing={missing} failed={failed} "
                  f"  {rate:.1f}/s  eta {eta / 60:.1f} min", flush=True)

    table = {}
    for ctype, positions in acc.items():
        table[ctype] = {}
        for pos, values in positions.items():
            values.sort()
            n = len(values)
            table[ctype][pos] = {
                "median_rsa": round(values[n // 2], 4),
                "mean_rsa": round(sum(values) / n, 4),
                "n": n,
            }

    meta = {
        "source": "SAbDab2 ML dataset, DOI 10.5281/zenodo.23343734 (CC-BY-4.0)",
        "split": args.split,
        "structures_used": done,
        "selection": "one representative per ab_cluster",
        "antigen_excluded": True,
        "n_points": args.n_points,
        "note": "rsa = SASA / Tien et al. 2013 max-ASA; antibody chains only",
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", encoding="utf-8") as fh:
        json.dump({"_meta": meta, "positions": table}, fh, indent=1, sort_keys=True)

    print(f"\nwrote {out}  ({out.stat().st_size / 1024:.0f} kB)")
    print(f"  structures used {done}, missing {missing}, failed {failed}")
    print(f"  H positions {len(table['H'])}, L positions {len(table['L'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
