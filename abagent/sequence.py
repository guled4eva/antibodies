"""Parsing and validation of antibody sequences."""

from __future__ import annotations

from dataclasses import dataclass

AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")


class SequenceError(ValueError):
    """Raised when an input is not a usable protein sequence."""


def clean(seq: str) -> str:
    """Strip whitespace, digits and FASTA headers; uppercase; validate letters."""
    lines = [ln for ln in seq.strip().splitlines() if not ln.startswith(">")]
    s = "".join(ch for ln in lines for ch in ln if not ch.isspace() and not ch.isdigit())
    s = s.upper().rstrip("*")
    if not s:
        raise SequenceError("Empty sequence.")
    bad = sorted(set(s) - AMINO_ACIDS)
    if bad:
        raise SequenceError(f"Non-standard residues in sequence: {''.join(bad)}")
    return s


@dataclass
class Antibody:
    """An antibody variable region: heavy chain and/or light chain."""

    heavy: str | None = None
    light: str | None = None
    name: str = "query"

    def __post_init__(self) -> None:
        if self.heavy:
            self.heavy = clean(self.heavy)
        if self.light:
            self.light = clean(self.light)
        if not self.heavy and not self.light:
            raise SequenceError("Give at least a heavy or a light chain.")

    def chains(self) -> dict[str, str]:
        out = {}
        if self.heavy:
            out["H"] = self.heavy
        if self.light:
            out["L"] = self.light
        return out


def parse_fasta(text: str) -> dict[str, str]:
    """Parse FASTA text into {header: sequence}."""
    records: dict[str, list[str]] = {}
    current = None
    for ln in text.splitlines():
        ln = ln.strip()
        if not ln:
            continue
        if ln.startswith(">"):
            current = ln[1:].strip() or f"seq{len(records) + 1}"
            records[current] = []
        else:
            if current is None:
                current = "seq1"
                records[current] = []
            records[current].append(ln)
    return {k: clean("".join(v)) for k, v in records.items()}
