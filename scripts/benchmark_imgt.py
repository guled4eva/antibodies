#!/usr/bin/env python3
"""Score abagent.imgt against SAbDab2's ground-truth IMGT numbering.

Reports exact position accuracy overall and per IMGT region, plus the
per-sequence distribution — the overall figure hides that the errors are
bimodal rather than evenly spread.

Usage:
    python scripts/benchmark_imgt.py --csv splits_final/ab_split.csv --limit 200
"""

import argparse
import csv
import json
import random
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from abagent.imgt import number  # noqa: E402


def region(label):
    digits = "".join(c for c in label if c.isdigit())
    p = int(digits or 0)
    for limit, name in ((26, "FR1"), (38, "CDR1"), (55, "FR2"), (65, "CDR2"),
                        (104, "FR3"), (117, "CDR3")):
        if p <= limit:
            return name
    return "FR4"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", required=True)
    ap.add_argument("--split", default="test")
    ap.add_argument("--limit", type=int, default=200)
    ap.add_argument("--seed", type=int, default=11)
    args = ap.parse_args(argv)

    csv.field_size_limit(min(sys.maxsize, 2**31 - 1))
    rows = []
    with open(args.csv, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            if row.get("ab_split") == args.split and row.get("VH_numbering_list"):
                rows.append(row)
    if not rows:
        print("no rows with numbering found", file=sys.stderr)
        return 1
    random.seed(args.seed)
    sample = random.sample(rows, min(args.limit, len(rows)))

    ok, tot = Counter(), Counter()
    per_seq = []
    for row in sample:
        truth = json.loads(row["VH_numbering_list"])
        seq = "".join(aa for aa, _ in truth)
        mine = number(seq)
        want = [f"{n[0]}{(n[1] or ' ').strip()}" for _, n in truth]
        got = [f"{p}{i}" for _, p, i in mine]
        hits = 0
        for i in range(min(len(want), len(got))):
            reg = region(want[i]); tot[reg] += 1
            if want[i] == got[i]:
                ok[reg] += 1; hits += 1
        per_seq.append(hits / max(1, len(want)))

    print(f"{args.split} split, n={len(sample)}\n")
    print(f"{'region':6s} {'correct':>9s} {'total':>8s}  accuracy")
    for reg in ("FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4"):
        if tot[reg]:
            print(f"{reg:6s} {ok[reg]:9d} {tot[reg]:8d}  {100*ok[reg]/tot[reg]:5.1f}%")
    total, correct = sum(tot.values()), sum(ok.values())
    print(f"\noverall {100*correct/total:.1f}%")
    per_seq.sort()
    pick = lambda q: per_seq[min(len(per_seq) - 1, int(len(per_seq) * q))]  # noqa: E731
    print(f"per-sequence: p10 {100*pick(.10):.1f}%  median {100*pick(.50):.1f}%  "
          f"p90 {100*pick(.90):.1f}%")
    print(f"  >95% correct: {sum(1 for a in per_seq if a > .95)}/{len(per_seq)}")
    print(f"  <50% correct: {sum(1 for a in per_seq if a < .50)}/{len(per_seq)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
