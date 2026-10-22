"""Solvent-accessible surface area from mmCIF structures.

Dependency-free apart from numpy: a small mmCIF `_atom_site` reader plus a
vectorised Shrake-Rupley implementation. Avoids pulling in a full structural
biology stack, which matters because the only thing this package needs from a
structure is per-residue exposure.

Validated against Biopython's Bio.PDB.SASA.ShrakeRupley (see
scripts/validate_sasa.py).

Why exposure matters here: a sequence liability only becomes a real risk if the
residue is on the surface. A buried methionine does not oxidise.
"""

import math

import numpy as np

# van der Waals radii (Angstrom), Bondi 1964.
VDW_RADII = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "H": 1.20, "P": 1.80}
DEFAULT_RADIUS = 1.70
PROBE_RADIUS = 1.4  # water

# Maximum accessible surface area per residue, Tien et al. 2013 (theoretical),
# used to turn absolute SASA into relative solvent accessibility (RSA).
MAX_ASA = {
    "ALA": 129.0, "ARG": 274.0, "ASN": 195.0, "ASP": 193.0, "CYS": 167.0,
    "GLN": 225.0, "GLU": 223.0, "GLY": 104.0, "HIS": 224.0, "ILE": 197.0,
    "LEU": 201.0, "LYS": 236.0, "MET": 224.0, "PHE": 240.0, "PRO": 159.0,
    "SER": 155.0, "THR": 172.0, "TRP": 285.0, "TYR": 263.0, "VAL": 174.0,
}

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}


class Atom:
    __slots__ = ("chain", "resseq", "icode", "resname", "element", "xyz")

    def __init__(self, chain, resseq, resname, element, xyz, icode=""):
        self.chain = chain
        self.resseq = resseq
        # Insertion code. Antibodies rely on these heavily -- a CDR3 routinely
        # contains 112, 112A, 112B... as DISTINCT residues. Ignoring it merges
        # them into one, which silently inflates CDR3 exposure.
        self.icode = icode
        self.resname = resname
        self.element = element
        self.xyz = xyz

    @property
    def residue_key(self):
        return (self.chain, self.resseq, self.icode)


def parse_cif_atoms(path, include_hetatm=False):
    """Read the _atom_site loop of an mmCIF file into a list of Atom.

    Column order varies between files, so the header names are read rather
    than assumed. Waters and hydrogens are skipped; altloc keeps only the
    first conformer.
    """
    atoms = []
    cols = {}
    in_loop = False
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("_atom_site."):
                cols[line.strip().split(".", 1)[1]] = len(cols)
                in_loop = True
                continue
            if not in_loop:
                continue
            if not line or line.startswith(("#", "loop_", "_")):
                if atoms:
                    break
                continue
            parts = line.split()
            if len(parts) < len(cols):
                continue
            rec = parts[cols["group_PDB"]]
            if rec != "ATOM" and not (include_hetatm and rec == "HETATM"):
                continue
            resname = parts[cols["label_comp_id"]]
            if resname == "HOH":
                continue
            element = parts[cols["type_symbol"]].strip().upper()
            if element == "H":
                continue
            altloc = parts[cols["label_alt_id"]]
            if altloc not in (".", "?", "A"):
                continue
            try:
                xyz = (float(parts[cols["Cartn_x"]]),
                       float(parts[cols["Cartn_y"]]),
                       float(parts[cols["Cartn_z"]]))
            except ValueError:
                continue
            icode = parts[cols["pdbx_PDB_ins_code"]] if "pdbx_PDB_ins_code" in cols else "?"
            atoms.append(Atom(
                chain=parts[cols.get("auth_asym_id", cols["label_asym_id"])],
                resseq=parts[cols.get("auth_seq_id", cols["label_seq_id"])],
                resname=resname, element=element, xyz=xyz,
                icode="" if icode in (".", "?") else icode,
            ))
    return atoms


def _sphere_points(n):
    """n roughly uniform points on the unit sphere (golden-spiral)."""
    i = np.arange(n, dtype=float) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    theta = math.pi * (1 + 5 ** 0.5) * i
    return np.column_stack((np.cos(theta) * np.sin(phi),
                            np.sin(theta) * np.sin(phi),
                            np.cos(phi)))


def shrake_rupley(coords, radii, n_points=100, probe=PROBE_RADIUS):
    """Per-atom solvent-accessible surface area, in square Angstrom."""
    coords = np.asarray(coords, dtype=float)
    radii = np.asarray(radii, dtype=float)
    n = len(coords)
    if n == 0:
        return np.zeros(0)
    expanded = radii + probe
    sphere = _sphere_points(n_points)
    # pairwise distances; structures here are small enough for a dense matrix
    diff = coords[:, None, :] - coords[None, :, :]
    dist = np.sqrt((diff ** 2).sum(-1))
    cutoff = expanded[:, None] + expanded[None, :]
    neighbours = (dist < cutoff) & ~np.eye(n, dtype=bool)

    areas = np.empty(n)
    for i in range(n):
        idx = np.flatnonzero(neighbours[i])
        pts = coords[i] + sphere * expanded[i]
        if idx.size:
            d2 = ((pts[:, None, :] - coords[idx][None, :, :]) ** 2).sum(-1)
            accessible = ~(d2 < (expanded[idx] ** 2)[None, :]).any(axis=1)
        else:
            accessible = np.ones(n_points, dtype=bool)
        areas[i] = (4.0 * math.pi * expanded[i] ** 2) * accessible.mean()
    return areas


def residue_sasa(atoms, n_points=100):
    """Aggregate per-atom SASA into per-residue absolute and relative values.

    Returns a list of dicts ordered as in the file:
      {chain, resseq, resname, aa, sasa, rsa}
    RSA is SASA / MAX_ASA[resname], i.e. 0 = fully buried, ~1 = fully exposed.
    SASA is computed in the context of the whole file, so antigen chains
    correctly bury the paratope.
    """
    if not atoms:
        return []
    coords = [a.xyz for a in atoms]
    radii = [VDW_RADII.get(a.element, DEFAULT_RADIUS) for a in atoms]
    areas = shrake_rupley(coords, radii, n_points=n_points)

    out = []
    index = {}
    for atom, area in zip(atoms, areas):
        key = atom.residue_key  # includes the insertion code
        if key not in index:
            index[key] = len(out)
            out.append({"chain": atom.chain, "resseq": atom.resseq,
                        "icode": atom.icode, "resname": atom.resname,
                        "aa": THREE_TO_ONE.get(atom.resname, "X"),
                        "sasa": 0.0, "rsa": 0.0})
        out[index[key]]["sasa"] += float(area)
    for res in out:
        max_asa = MAX_ASA.get(res["resname"])
        res["rsa"] = res["sasa"] / max_asa if max_asa else float("nan")
    return out


def exposure_from_cif(path, chains=None, n_points=100):
    """Convenience: parse a .cif and return per-residue exposure.

    chains: iterable of auth chain ids to keep in the output. SASA is still
    computed with every chain present, so burial by the antigen is preserved.
    """
    atoms = parse_cif_atoms(path)
    residues = residue_sasa(atoms, n_points=n_points)
    if chains is not None:
        keep = set(chains)
        residues = [r for r in residues if r["chain"] in keep]
    return residues
