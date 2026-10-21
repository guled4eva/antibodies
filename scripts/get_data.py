#!/usr/bin/env python3
"""Download the SAbDab2 ML dataset from Zenodo into ./data.

The dataset is CC-BY-4.0 (Capel et al., 2026; DOI 10.5281/zenodo.23343734).
It is ~958 MB, so it is NOT committed to git — this script fetches it on demand.

Usage:
    python scripts/get_data.py            # download + verify + extract
    python scripts/get_data.py --no-extract
"""

from __future__ import annotations

import argparse
import hashlib
import sys
import tarfile
import urllib.request
from pathlib import Path

RECORD = "23343734"
FILENAME = "splits.tar.gz"
URL = f"https://zenodo.org/records/{RECORD}/files/{FILENAME}?download=1"
EXPECTED_MD5 = "75662fe08cb542177971fa3ded2cd800"
DATA_DIR = Path(__file__).resolve().parent.parent / "data"


def download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"Downloading {url}\n  -> {dest}")

    def _progress(block: int, block_size: int, total: int) -> None:
        if total > 0:
            pct = min(100, block * block_size * 100 // total)
            print(f"\r  {pct:3d}%", end="", flush=True)

    urllib.request.urlretrieve(url, dest, _progress)  # noqa: S310 (trusted DOI host)
    print()


def verify_md5(path: Path, expected: str) -> bool:
    h = hashlib.md5()  # noqa: S324 (Zenodo publishes md5; integrity check only)
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    got = h.hexdigest()
    ok = got == expected
    print(f"  md5 {'OK' if ok else 'MISMATCH'}: {got}")
    return ok


def extract(archive: Path, out_dir: Path) -> None:
    print(f"Extracting -> {out_dir}")
    with tarfile.open(archive, "r:gz") as tar:
        tar.extractall(out_dir, filter="data")  # filter guards against path traversal
    print("  done")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--no-extract", action="store_true", help="download only")
    ap.add_argument("--skip-md5", action="store_true", help="skip integrity check")
    args = ap.parse_args(argv)

    archive = DATA_DIR / FILENAME
    if not archive.exists():
        try:
            download(URL, archive)
        except Exception as exc:
            print(f"\nDownload failed: {exc}", file=sys.stderr)
            print("If you are behind a restricted network, download the file "
                  f"manually from https://zenodo.org/records/{RECORD} and place it at "
                  f"{archive}", file=sys.stderr)
            return 1
    else:
        print(f"Already present: {archive}")

    if not args.skip_md5 and not verify_md5(archive, EXPECTED_MD5):
        print("Integrity check failed; delete the file and re-run.", file=sys.stderr)
        return 1

    if not args.no_extract:
        extract(archive, DATA_DIR)

    print("\nReady. Data in", DATA_DIR)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
