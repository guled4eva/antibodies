"""Developability analysis of antibody variable regions.

Pure-Python, no external dependencies. Flags sequence liabilities and
computes basic biophysical descriptors that a drug-discovery team uses to
judge whether a candidate is worth carrying forward.

References:
  - Deamidation (NG/NS/NT) and isomerisation (DG/DS/DT): Lu et al. 2019,
    mAbs 11(1).
  - Kyte & Doolittle hydropathy: J. Mol. Biol. 157 (1982) 105-132.
  - pKa set: EMBOSS.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

# --- liability motifs: name -> (compiled regex, why it matters) ---------------
_LIABILITY_PATTERNS: dict[str, tuple[str, str]] = {
    "N-glycosylation": (r"N[^P][ST]", "N-linked glycosylation sequon; adds heterogeneity"),
    "Deamidation (high)": (r"N[GS]", "Fast Asn deamidation; charge variant + activity loss"),
    "Deamidation (low)": (r"N[THA]", "Slower Asn deamidation site"),
    "Isomerisation": (r"D[GSTDH]", "Asp isomerisation; can sit in CDR and cut potency"),
    "Fragmentation": (r"DP", "Acid-labile Asp-Pro peptide bond"),
    "Oxidation (Met)": (r"M", "Met oxidation, esp. if solvent-exposed"),
    "Oxidation (Trp)": (r"W", "Trp oxidation under stress/light"),
    "Free cysteine": (r"C", "Unpaired Cys -> mispairing / aggregation"),
}

# Kyte-Doolittle hydropathy
_KD = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

# pKa values (EMBOSS) for isoelectric point
_PKA_POS = {"K": 10.8, "R": 12.5, "H": 6.5}
_PKA_NEG = {"D": 3.9, "E": 4.1, "C": 8.5, "Y": 10.1}
_PKA_NTERM = 8.6
_PKA_CTERM = 3.6


@dataclass
class Liability:
    kind: str
    chain: str
    position: int  # 1-based index into the raw sequence
    motif: str
    note: str


@dataclass
class ChainReport:
    chain: str
    length: int
    pI: float
    net_charge_pH7: float
    gravy: float  # grand average of hydropathy
    liabilities: list[Liability] = field(default_factory=list)


def _net_charge(seq: str, ph: float) -> float:
    charge = 0.0
    charge += 1.0 / (1.0 + 10 ** (ph - _PKA_NTERM))
    charge -= 1.0 / (1.0 + 10 ** (_PKA_CTERM - ph))
    for aa, pka in _PKA_POS.items():
        charge += seq.count(aa) / (1.0 + 10 ** (ph - pka))
    for aa, pka in _PKA_NEG.items():
        charge -= seq.count(aa) / (1.0 + 10 ** (pka - ph))
    return charge


def _isoelectric_point(seq: str) -> float:
    lo, hi = 0.0, 14.0
    for _ in range(100):
        mid = (lo + hi) / 2
        if _net_charge(seq, mid) > 0:
            lo = mid
        else:
            hi = mid
    return round((lo + hi) / 2, 2)


def _gravy(seq: str) -> float:
    return round(sum(_KD[a] for a in seq) / len(seq), 3)


def find_liabilities(seq: str, chain: str) -> list[Liability]:
    found: list[Liability] = []
    for kind, (pattern, note) in _LIABILITY_PATTERNS.items():
        for m in re.finditer(f"(?={pattern})", seq):  # overlapping
            found.append(
                Liability(kind=kind, chain=chain, position=m.start() + 1,
                          motif=seq[m.start(): m.start() + len(m.group(0) or "X") + 2][:3],
                          note=note)
            )
    found.sort(key=lambda x: x.position)
    return found


def analyse_chain(seq: str, chain: str) -> ChainReport:
    return ChainReport(
        chain=chain,
        length=len(seq),
        pI=_isoelectric_point(seq),
        net_charge_pH7=round(_net_charge(seq, 7.0), 2),
        gravy=_gravy(seq),
        liabilities=find_liabilities(seq, chain),
    )
