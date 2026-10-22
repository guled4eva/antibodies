"""Tests for the bundled SASA implementation (no network, no big files)."""

import math
import tempfile
from pathlib import Path

from abagent.sasa import PROBE_RADIUS, parse_cif_atoms, residue_sasa, shrake_rupley

MINI_CIF = """data_TEST
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_asym_id
_atom_site.pdbx_PDB_model_num
ATOM 1 N N . GLN A 1 1 ? 0.000 0.000 0.000 1 10.0 ? 112 H 1
ATOM 2 C CA . GLN A 1 1 ? 1.500 0.000 0.000 1 10.0 ? 112 H 1
ATOM 3 N N . TYR A 1 2 A 20.000 0.000 0.000 1 10.0 ? 112 H 1
ATOM 4 C CA . TYR A 1 2 A 21.500 0.000 0.000 1 10.0 ? 112 H 1
ATOM 5 O O . HOH A 1 3 ? 40.000 0.000 0.000 1 10.0 ? 300 H 1
#
"""


def _write(text):
    path = Path(tempfile.mkdtemp()) / "mini.cif"
    path.write_text(text)
    return path


def test_parser_reads_atoms_and_skips_water():
    atoms = parse_cif_atoms(_write(MINI_CIF))
    assert len(atoms) == 4  # water excluded
    assert {a.resname for a in atoms} == {"GLN", "TYR"}


def test_parser_separates_insertion_codes():
    """112 and 112A are DIFFERENT residues; merging them inflates CDR3 exposure."""
    residues = residue_sasa(parse_cif_atoms(_write(MINI_CIF)), n_points=50)
    assert len(residues) == 2, "insertion code must not be collapsed"
    keys = {(r["resseq"], r["icode"]) for r in residues}
    assert keys == {("112", ""), ("112", "A")}


def test_isolated_atom_matches_analytic_area():
    r = 1.70
    area = shrake_rupley([(0.0, 0.0, 0.0)], [r], n_points=2000)[0]
    expected = 4 * math.pi * (r + PROBE_RADIUS) ** 2
    assert abs(area - expected) / expected < 0.01


def test_enclosed_atom_is_buried():
    """An atom ringed by close neighbours should have near-zero exposure."""
    coords = [(0.0, 0.0, 0.0)]
    for dx in (-2.0, 2.0):
        for dy in (-2.0, 2.0):
            for dz in (-2.0, 2.0):
                coords.append((dx, dy, dz))
    areas = shrake_rupley(coords, [1.7] * len(coords), n_points=500)
    assert areas[0] < 0.05 * areas[1]


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print("ok", name)
