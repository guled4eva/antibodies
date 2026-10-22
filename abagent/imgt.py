"""Assign IMGT positions to a variable-domain sequence, without ANARCI.

Why this exists: the structural exposure table is keyed by IMGT position, but a
user pastes a bare sequence, so something has to say which IMGT position residue
83 of that string is. `anarcii` does this properly with germline HMMs; this is a
dependency-free approximation that also runs in a browser, and its accuracy is
measured rather than assumed (scripts/benchmark_imgt.py).

Method: numbering by transfer. Antibody frameworks are highly conserved, and the
SAbDab2 dataset ships thousands of sequences already numbered by IMGT. Pick the
most similar reference, align to it, and carry its numbering across. Residues
that align to a gap in the reference take insertion codes.

MEASURED, AND NOT GOOD ENOUGH TO GATE ON (scripts/benchmark_imgt.py, held-out
test split, n=200, residue positions exactly correct):

    hand-coded IMGT gap rules      14.9%
    numbering by transfer (this)   74.1%   framework 80.0%, CDR3 39.9%

The error is bimodal: ~20% of sequences come out almost entirely wrong (0-1%
correct, the alignment anchors badly), and the rest land at 82-100%, median
91.6%. Mislabelling a fifth of inputs is not a basis for telling someone a
liability is buried, so this is NOT wired into the liability report.

Use `pip install anarcii` for real numbering; it uses germline HMMs and is the
right tool for this. This module stays as a dependency-free fallback and as a
record of what the cheap approach actually scores.
"""

import json
from pathlib import Path

REFERENCE_FILE = Path(__file__).parent / "data" / "imgt_references.json"

# Needleman-Wunsch parameters. Frameworks are near-identical between antibodies,
# so a simple identity matrix with an affine-ish gap cost is enough.
MATCH, MISMATCH, GAP_OPEN, GAP_EXTEND = 2, -1, -6, -1

_REFS = None


def _load_references():
    global _REFS
    if _REFS is None:
        try:
            with open(REFERENCE_FILE, encoding="utf-8") as fh:
                _REFS = json.load(fh)["references"]
        except (OSError, ValueError, KeyError):
            _REFS = []
    return _REFS


def _kmers(seq, k=3):
    return {seq[i:i + k] for i in range(len(seq) - k + 1)}


def _pick_reference(seq, refs):
    """Cheap k-mer similarity to avoid aligning against every reference."""
    q = _kmers(seq)
    best, best_score = None, -1.0
    for ref in refs:
        rk = ref.get("_k")
        if rk is None:
            rk = ref["_k"] = _kmers(ref["seq"])
        score = len(q & rk) / max(1, len(q | rk))
        if score > best_score:
            best, best_score = ref, score
    return best


def _align(a, b):
    """Global alignment; returns pairs of (index_in_a | None, index_in_b | None)."""
    n, m = len(a), len(b)
    score = [[0] * (m + 1) for _ in range(n + 1)]
    ptr = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        score[i][0] = GAP_OPEN + (i - 1) * GAP_EXTEND
        ptr[i][0] = 1
    for j in range(1, m + 1):
        score[0][j] = GAP_OPEN + (j - 1) * GAP_EXTEND
        ptr[0][j] = 2
    for i in range(1, n + 1):
        ai = a[i - 1]
        row, prev = score[i], score[i - 1]
        prow = ptr[i]
        for j in range(1, m + 1):
            diag = prev[j - 1] + (MATCH if ai == b[j - 1] else MISMATCH)
            up = prev[j] + (GAP_EXTEND if ptr[i - 1][j] == 1 else GAP_OPEN)
            left = row[j - 1] + (GAP_EXTEND if prow[j - 1] == 2 else GAP_OPEN)
            best = diag
            which = 0
            if up > best:
                best, which = up, 1
            if left > best:
                best, which = left, 2
            row[j] = best
            prow[j] = which
    pairs = []
    i, j = n, m
    while i > 0 or j > 0:
        if i > 0 and j > 0 and ptr[i][j] == 0:
            pairs.append((i - 1, j - 1)); i -= 1; j -= 1
        elif i > 0 and (j == 0 or ptr[i][j] == 1):
            pairs.append((i - 1, None)); i -= 1
        else:
            pairs.append((None, j - 1)); j -= 1
    pairs.reverse()
    return pairs


def number(seq):
    """Return [(residue, imgt_position, insertion_code), ...] or [] if unavailable."""
    seq = (seq or "").upper()
    refs = _load_references()
    if not seq or not refs:
        return []
    ref = _pick_reference(seq, refs)
    if ref is None:
        return []
    pairs = _align(seq, ref["seq"])
    out = []
    last_pos, ins_n = None, 0
    for qi, ri in pairs:
        if qi is None:
            continue                      # reference residue the query lacks
        if ri is None:
            # query residue with no reference partner -> insertion after last_pos
            if last_pos is None:
                continue
            out.append((seq[qi], last_pos, chr(ord("A") + ins_n)))
            ins_n += 1
            continue
        pos, ins = ref["num"][ri]
        out.append((seq[qi], pos, ins))
        last_pos, ins_n = pos, 0
    return out


def position_map(seq):
    """{0-based sequence index: "positionInsertion"} for the given sequence."""
    seq = (seq or "").upper()
    refs = _load_references()
    if not seq or not refs:
        return {}
    ref = _pick_reference(seq, refs)
    if ref is None:
        return {}
    out = {}
    last_pos, ins_n = None, 0
    for qi, ri in _align(seq, ref["seq"]):
        if qi is None:
            continue
        if ri is None:
            if last_pos is not None:
                out[qi] = f"{last_pos}{chr(ord('A') + ins_n)}"
                ins_n += 1
            continue
        pos, ins = ref["num"][ri]
        out[qi] = f"{pos}{ins}"
        last_pos, ins_n = pos, 0
    return out
