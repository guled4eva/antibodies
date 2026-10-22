#!/usr/bin/env python3
"""Check abagent.sasa against Biopython's reference Shrake-Rupley.

Biopython is NOT a runtime dependency -- it is used here only to verify the
bundled implementation. Run this if you change anything in abagent/sasa.py.

Usage:
    pip install biopython
    python scripts/validate_sasa.py path/to/structure.cif [more.cif ...]
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np  # noqa: E402

from abagent.sasa import parse_cif_atoms, residue_sasa  # noqa: E402


N_POINTS = 500  # compare where both implementations have converged


def reference(path):
    from Bio.PDB import MMCIFParser
    from Bio.PDB.SASA import ShrakeRupley

    structure = MMCIFParser(QUIET=True).get_structure("s", path)
    model = next(structure.get_models())
    ShrakeRupley(probe_radius=1.4, n_points=N_POINTS).compute(model, level="R")
    out = {}
    for chain in model:
        for res in chain:
            if res.id[0] != " ":  # skip waters / hetero
                continue
            out[(chain.id, str(res.id[1]), res.id[2].strip())] = res.sasa
    return out


def main(argv):
    paths = argv or []
    if not paths:
        print(__doc__)
        return 1
    all_rel = []
    for path in paths:
        mine = {(r["chain"], r["resseq"], r.get("icode", "")): r["sasa"]
                for r in residue_sasa(parse_cif_atoms(path), n_points=N_POINTS)}
        try:
            theirs = reference(path)
        except ImportError:
            print("biopython not installed; run: pip install biopython")
            return 1
        shared = sorted(set(mine) & set(theirs))
        if not shared:
            print(f"{Path(path).name}: no residues in common -- parser mismatch!")
            return 1
        a = np.array([mine[k] for k in shared])
        b = np.array([theirs[k] for k in shared])
        denom = np.where(b > 1.0, b, np.nan)
        rel = np.abs(a - b) / denom
        rel = rel[~np.isnan(rel)]
        all_rel.append(rel)
        corr = np.corrcoef(a, b)[0, 1]
        print(f"{Path(path).name:34s} n={len(shared):4d}  "
              f"r={corr:.5f}  median|rel err|={np.median(rel) * 100:5.2f}%  "
              f"max={np.max(rel) * 100:5.1f}%  "
              f"mine/theirs total={a.sum() / b.sum():.4f}")
    pooled = np.concatenate(all_rel)
    ok = np.median(pooled) < 0.02
    print(f"\npooled median relative error: {np.median(pooled) * 100:.2f}%  "
          f"-> {'OK' if ok else 'TOO HIGH'}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
