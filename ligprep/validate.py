"""Checks that run on every ligand before its PDBQT is written."""
from __future__ import annotations

import re
from dataclasses import dataclass, field

from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors

from .errors import LigprepError

RDLogger.DisableLog("rdApp.*")
HYDROGEN_TYPES = {"H", "HD", "HS"}
FLAT_TOLERANCE = 1e-3  # angstrom; a 2D structure has every z exactly 0
MANY_TORSIONS = 20
CHARGE_TOLERANCE = 0.1


@dataclass
class Report:
    formula: str = ""
    heavy_atoms: int = 0
    torsdof: int | None = None
    charge_sum: float = 0.0
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.errors


def parse_formula(formula: str) -> dict[str, int]:
    counts: dict[str, int] = {}
    for el, n in re.findall(r"([A-Z][a-z]?)(\d*)", formula):
        counts[el] = counts.get(el, 0) + (int(n) if n else 1)
    return counts


def read_sdf(sdf_text: str) -> Chem.Mol:
    mol = Chem.MolFromMolBlock(sdf_text, removeHs=False)
    if mol is None or mol.GetNumConformers() == 0:
        raise LigprepError("could not read a 3D structure from the SDF")
    return mol


def mol_formula(mol: Chem.Mol) -> str:
    return re.sub(r"[+-]\d*$", "", rdMolDescriptors.CalcMolFormula(mol))  # drop RDKit's charge suffix


def pdbqt_stats(pdbqt: str) -> dict:
    heavy = polar_h = 0
    charge = 0.0
    torsdof = None
    for line in pdbqt.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            atom_type = line[77:79].strip()
            if atom_type in HYDROGEN_TYPES:
                polar_h += 1
            else:
                heavy += 1
            try:
                charge += float(line[70:76])
            except ValueError:
                pass
        elif line.startswith("TORSDOF"):
            torsdof = int(line.split()[1])
    return {"heavy": heavy, "polar_h": polar_h, "charge": charge, "torsdof": torsdof}


def validate(mol: Chem.Mol, pubchem_formula: str, pdbqt: str, allow_flat: bool = False) -> Report:
    rep = Report(formula=mol_formula(mol))
    expected = parse_formula(pubchem_formula)
    rep.heavy_atoms = sum(n for el, n in expected.items() if el != "H")

    if rep.formula != pubchem_formula:
        rep.errors.append(f"formula mismatch: SDF {rep.formula} vs PubChem {pubchem_formula}")

    zs = [mol.GetConformer().GetAtomPosition(i).z for i in range(mol.GetNumAtoms())]
    if not allow_flat and max(abs(z) for z in zs) < FLAT_TOLERANCE and rep.heavy_atoms > 3:
        rep.errors.append("all z coordinates are 0: this is a 2D structure, not a 3D conformer")

    stats = pdbqt_stats(pdbqt)
    rep.torsdof, rep.charge_sum = stats["torsdof"], round(stats["charge"], 3)
    if stats["heavy"] != rep.heavy_atoms:
        rep.errors.append(f"PDBQT has {stats['heavy']} heavy atoms, formula says {rep.heavy_atoms}")
    if rep.torsdof is None:
        rep.errors.append("PDBQT has no TORSDOF line")
    elif rep.torsdof > MANY_TORSIONS:
        rep.warnings.append(f"TORSDOF {rep.torsdof}: many rotatable bonds, docking may be unreliable")

    formal = sum(a.GetFormalCharge() for a in mol.GetAtoms())
    if abs(stats["charge"] - formal) > CHARGE_TOLERANCE:
        rep.warnings.append(f"partial charges sum to {stats['charge']:+.2f}, formal charge is {formal:+d}")
    return rep
