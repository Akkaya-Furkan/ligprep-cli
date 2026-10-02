import pytest
from rdkit.Geometry import Point3D

from ligprep.errors import LigprepError
from ligprep.validate import MANY_TORSIONS, parse_formula, pdbqt_stats, read_sdf, validate


def fake_pdbqt(heavy: int = 13, charges: float = 0.0, torsdof: int | None = 4, polar_h: int = 1) -> str:
    """PDBQT text with the fixed columns Open Babel writes (charge 71-76, type 78-79)."""
    lines = []
    for i in range(heavy):
        q = charges if i == 0 else 0.0
        lines.append(f"ATOM  {i + 1:5d}  C   UNL     1    {i:8.3f}{1.0:8.3f}{0.5:8.3f}  0.00  0.00    {q:+6.3f} C ")
    for j in range(polar_h):
        lines.append(f"ATOM  {heavy + j + 1:5d}  H   UNL     1    {0.0:8.3f}{0.0:8.3f}{0.0:8.3f}  0.00  0.00    +0.000 HD")
    if torsdof is not None:
        lines.append(f"TORSDOF {torsdof}")
    return "\n".join(lines) + "\n"


@pytest.fixture
def mol(aspirin_sdf):
    return read_sdf(aspirin_sdf)


def test_parse_formula():
    assert parse_formula("C9H8O4") == {"C": 9, "H": 8, "O": 4}
    assert parse_formula("CH4N2Cl") == {"C": 1, "H": 4, "N": 2, "Cl": 1}


def test_pdbqt_stats_separates_polar_hydrogens():
    s = pdbqt_stats(fake_pdbqt(heavy=13, polar_h=2, torsdof=4))
    assert (s["heavy"], s["polar_h"], s["torsdof"]) == (13, 2, 4)


def test_read_sdf_rejects_garbage():
    with pytest.raises(LigprepError):
        read_sdf("this is not an sdf")


def test_a_good_ligand_passes(mol):
    rep = validate(mol, "C9H8O4", fake_pdbqt())
    assert rep.ok and not rep.warnings
    assert (rep.formula, rep.heavy_atoms, rep.torsdof) == ("C9H8O4", 13, 4)


def test_formula_mismatch_fails(mol):
    rep = validate(mol, "C9H8O5", fake_pdbqt(heavy=14))
    assert not rep.ok
    assert any("formula mismatch" in m for m in rep.errors)


def test_heavy_atom_count_mismatch_fails(mol):
    rep = validate(mol, "C9H8O4", fake_pdbqt(heavy=12))
    assert any("12 heavy atoms" in m for m in rep.errors)


def flatten(mol):
    conf = mol.GetConformer()
    for i in range(mol.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, Point3D(p.x, p.y, 0.0))
    return mol


def test_flat_structure_fails_unless_allowed(mol):
    flat = flatten(mol)
    assert any("2D structure" in m for m in validate(flat, "C9H8O4", fake_pdbqt()).errors)
    assert validate(flat, "C9H8O4", fake_pdbqt(), allow_flat=True).ok


def test_real_conformer_is_not_flat(mol):
    assert not any("2D" in m for m in validate(mol, "C9H8O4", fake_pdbqt()).errors)


def test_missing_torsdof_fails(mol):
    rep = validate(mol, "C9H8O4", fake_pdbqt(torsdof=None))
    assert any("TORSDOF" in m for m in rep.errors)


def test_many_torsions_only_warn(mol):
    rep = validate(mol, "C9H8O4", fake_pdbqt(torsdof=MANY_TORSIONS + 1))
    assert rep.ok and any("rotatable" in m for m in rep.warnings)


def test_charge_far_from_formal_charge_warns(mol):
    rep = validate(mol, "C9H8O4", fake_pdbqt(charges=0.5))
    assert rep.ok and any("partial charges" in m for m in rep.warnings)
