#!/usr/bin/env python3
"""Measure the CDR finder against SAbDab2's ground-truth IMGT annotations.

Answers a question a reviewer will ask: how often does the CDR extraction
actually get the loop right? Compares abagent's heavy-chain CDR3 call against
the dataset's IMGT CDRH3 column.

Usage:
    python scripts/benchmark_cdrs.py                 # first 500 paired antibodies
    python scripts/benchmark_cdrs.py --limit 5000
    python scripts/benchmark_cdrs.py --split test    # held-out split only
    python scripts/benchmark_cdrs.py --show 20       # print example mismatches
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from abagent.dataset import DatasetError, find_split_csv, iter_records  # noqa: E402
from abagent.numbering import cdrs  # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-dir", default="data")
    ap.add_argument("--csv", default=None, help="specific split CSV to use")
    ap.add_argument("--limit", type=int, default=500, help="antibodies to test")
    ap.add_argument("--split", default=None, choices=["train", "test"])
    ap.add_argument("--show", type=int, default=5, help="example mismatches to print")
    args = ap.parse_args(argv)

    try:
        path = args.csv or find_split_csv(args.data_dir)
    except DatasetError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    print(f"Using {path}\n")

    exact = contained = 0
    total = 0
    empty = 0
    mismatches = []

    try:
        records = iter_records(path, limit=args.limit, split=args.split,
                               paired_only=True)
        for rec in records:
            if not rec.heavy or not rec.cdrh3:
                continue
            total += 1
            predicted = cdrs(rec.heavy, "H").cdr3
            truth = rec.cdrh3
            if not predicted:
                empty += 1
            elif predicted == truth:
                exact += 1
            elif truth in predicted or predicted in truth:
                contained += 1
                if len(mismatches) < args.show:
                    mismatches.append((rec.instance, truth, predicted, "overlap"))
            elif len(mismatches) < args.show:
                mismatches.append((rec.instance, truth, predicted, "wrong"))
    except DatasetError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print("\n(interrupted)")

    if not total:
        print("No usable records found.", file=sys.stderr)
        return 1

    wrong = total - exact - contained - empty
    pct = lambda n: f"{100 * n / total:5.1f}%"  # noqa: E731
    print(f"Heavy-chain CDR3 vs IMGT ground truth   (n = {total})")
    print(f"  exact match      {exact:6d}   {pct(exact)}")
    print(f"  partial overlap  {contained:6d}   {pct(contained)}")
    print(f"  wrong            {wrong:6d}   {pct(wrong)}")
    print(f"  found nothing    {empty:6d}   {pct(empty)}")

    if mismatches:
        print(f"\nExamples (truth vs predicted):")
        for inst, truth, pred, kind in mismatches:
            print(f"  [{kind:7s}] {inst}")
            print(f"      truth: {truth}")
            print(f"      ours : {pred}")

    engine = cdrs("EVQLVESGGG", "H").engine
    print(f"\nEngine in use: {engine}")
    if engine == "heuristic":
        print("Install the real numbering engine to improve this:")
        print("    pip install anarcii")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
