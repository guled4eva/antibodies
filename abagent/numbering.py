"""Antibody numbering and CDR identification.

Primary engine is `anarcii` (modern ANARCI re-implementation) when it is
installed. Because that pulls in PyTorch, the module also ships a
dependency-free fallback that anchors CDRs on the conserved framework
residues. The fallback is approximate and is only meant to keep the
skeleton runnable until the real engine is installed on hackathon day.
"""

from __future__ import annotations

import re
from dataclasses import dataclass


# Tuned on the SAbDab2 train split only; see _cdrs_heuristic for the numbers.
# A CDR3 longer than this means the start anchor is wrong, not that the loop is
# genuinely long (the sweep was flat from 20-50, so this is not a sharp fit).
_LONG_SPAN = 25
# Conserved framework-3 motif immediately preceding Cys104, e.g. ...YYC / YFC.
_CYS104_MOTIF = r"[YF]C"


@dataclass
class CDRSet:
    chain: str
    scheme: str
    cdr1: str
    cdr2: str
    cdr3: str
    engine: str  # "anarcii" or "heuristic"


def _anarcii_available() -> bool:
    try:
        import anarcii  # noqa: F401
        return True
    except Exception:
        return False


def _cdrs_with_anarcii(seq: str, chain: str, scheme: str = "imgt") -> CDRSet:
    from anarcii import Anarcii

    model = Anarcii(seq_type="antibody", mode="accuracy", scheme=scheme, verbose=False)
    numbered = model.number([seq])[0]
    # IMGT CDR position ranges (inclusive)
    ranges = {"cdr1": (27, 38), "cdr2": (56, 65), "cdr3": (105, 117)}
    loops = {k: [] for k in ranges}
    for (pos, _ins), aa in numbered["numbering"]:
        if aa == "-":
            continue
        for name, (lo, hi) in ranges.items():
            if lo <= pos <= hi:
                loops[name].append(aa)
    return CDRSet(
        chain=chain, scheme=scheme,
        cdr1="".join(loops["cdr1"]),
        cdr2="".join(loops["cdr2"]),
        cdr3="".join(loops["cdr3"]),
        engine="anarcii",
    )


def _cdrs_heuristic(seq: str, chain: str) -> CDRSet:
    """Approximate CDRs by anchoring on conserved framework residues.

    Heavy: CDR3 sits between the second conserved Cys and the WGxG motif.
    This is deliberately rough; use the anarcii path for anything real.
    """
    # Accuracy on the SAbDab2 HELD-OUT test split (n=2479, exact CDRH3 match),
    # measured with scripts/benchmark_cdrs.py:
    #     98.99%   this rule
    #     96.41%   without the re-anchor step below
    # Parameters (_LONG_SPAN, _CYS104_MOTIF) were chosen on the TRAIN split
    # only; train scored 99.00% vs 98.99% test, so they are not overfitted.
    #
    # Also tried and rejected, both scoring worse on train:
    #   - anchoring on the last J motif and walking back to the nearest
    #     cysteine (94.5%): gives up entirely when no WGxG/FGxG is present.
    #   - always starting at the last cysteine before the end (94.3%): CDR3
    #     loops often contain their own cysteine, which truncates the loop.
    # Re-run the benchmark before changing any of this.
    cys = [m.start() for m in re.finditer("C", seq)]
    cdr3 = ""
    if len(cys) >= 2:
        # IMGT CDR3 starts at position 105, immediately after the conserved
        # Cys104 — not after the CAR/CAS motif. The second cysteine is usually
        # Cys104, which is why this is the primary anchor.
        start = cys[1] + 1
        m = re.search(r"WG.G", seq[start:])
        if m is None:
            m = re.search(r"FG.G", seq[start:])  # light-chain J motif
        end = start + m.start() if m else min(start + 15, len(seq))

        # If that yields an implausibly long loop, the second cysteine was not
        # Cys104 (the chain carries an extra cysteine earlier). Re-anchor on
        # the conserved [YF]-Cys of framework 3, which marks Cys104 directly.
        if end - start > _LONG_SPAN:
            anchors = [m.start() + 1 for m in re.finditer(_CYS104_MOTIF, seq[:end])]
            if anchors:
                start = anchors[-1] + 1
        cdr3 = seq[start:end]
    # CDR1/CDR2 are only crudely placed; left blank rather than guessed wrongly
    return CDRSet(chain=chain, scheme="heuristic", cdr1="", cdr2="",
                  cdr3=cdr3, engine="heuristic")


def cdrs(seq: str, chain: str, scheme: str = "imgt") -> CDRSet:
    """Return the CDR loops for a chain, using anarcii if available."""
    if _anarcii_available():
        try:
            return _cdrs_with_anarcii(seq, chain, scheme)
        except Exception:
            pass
    return _cdrs_heuristic(seq, chain)
