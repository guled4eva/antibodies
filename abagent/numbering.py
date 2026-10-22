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
    # Anchor on the J-region motif first (WGxG heavy / FGxG light), taking the
    # LAST occurrence since the J region is C-terminal, then walk back to the
    # nearest preceding cysteine — that is the conserved Cys104. Indexing by
    # "the second cysteine" fails on the many real antibodies that carry extra
    # cysteines, which showed up as huge framework-spanning calls when scored
    # against SAbDab2 (scripts/benchmark_cdrs.py).
    cdr3 = ""
    motifs = list(re.finditer(r"WG.G|FG.G", seq))
    if motifs:
        end = motifs[-1].start()
        preceding_cys = [i for i, aa in enumerate(seq[:end]) if aa == "C"]
        if preceding_cys:
            start = preceding_cys[-1] + 1  # IMGT CDR3 begins right after Cys104
            candidate = seq[start:end]
            # Real CDR3 loops run ~2-35 residues; anything else means the
            # anchors were wrong, and a wrong loop is worse than none.
            if 2 <= len(candidate) <= 35:
                cdr3 = candidate
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
