import csv

import responses

from ligprep import cli
from ligprep.pubchem import BASE

from .conftest import needs_obabel


def run(tmp_path, lines, *extra):
    lst = tmp_path / "list.txt"
    lst.write_text("\n".join(lines) + "\n")
    out = tmp_path / "out"
    code = cli.main(["--input", str(lst), "--out", str(out), "--delay", "0", *extra])
    return code, out


def manifest(out):
    with open(out / "manifest.csv") as fh:
        return list(csv.DictReader(fh))


@needs_obabel
def test_end_to_end_with_mocked_pubchem(tmp_path, mock_aspirin, capsys):
    mock_aspirin()
    code, out = run(tmp_path, ["aspirin"])
    assert code == 0
    row = manifest(out)[0]
    assert (row["status"], row["cid"], row["formula"]) == ("ok", "2244", "C9H8O4")
    assert row["smiles"] == "CC(=O)OC1=CC=CC=C1C(=O)O"
    assert (out / "aspirin_CID2244.pdbqt").exists()
    assert "1/1 ligands prepared" in capsys.readouterr().out


@needs_obabel
def test_formula_mismatch_writes_no_pdbqt_and_exits_1(tmp_path, mock_aspirin):
    mock_aspirin(formula="C9H8O5")
    code, out = run(tmp_path, ["aspirin"])
    assert code == 1
    row = manifest(out)[0]
    assert row["status"] == "failed" and "formula mismatch" in row["messages"]
    assert not list(out.glob("*.pdbqt"))


@needs_obabel
def test_smiles_input_resolves_to_the_same_cid(tmp_path, mock_aspirin):
    mock_aspirin()
    code, out = run(tmp_path, ["smiles:CC(=O)Oc1ccccc1C(=O)O"])
    assert code == 0 and (out / "smiles_CID2244.pdbqt").exists()


@needs_obabel
def test_keep_sdf(tmp_path, mock_aspirin):
    mock_aspirin()
    _, out = run(tmp_path, ["aspirin"], "--keep-sdf")
    assert (out / "aspirin_CID2244.sdf").exists()


@needs_obabel
def test_one_failure_does_not_stop_the_run(tmp_path, mock_aspirin, rsps):
    mock_aspirin()
    rsps.add(responses.GET, f"{BASE}/compound/name/zzz/cids/TXT", status=404)
    code, out = run(tmp_path, ["zzz", "aspirin"])
    assert code == 1
    assert [r["status"] for r in manifest(out)] == ["failed", "ok"]


def test_empty_list_exits_2(tmp_path, capsys):
    code, _ = run(tmp_path, ["# nothing"])
    assert code == 2 and "empty" in capsys.readouterr().err


def test_malformed_cid_exits_2_with_the_line_number(tmp_path, capsys):
    code, _ = run(tmp_path, ["aspirin", "cid:abc"])
    assert code == 2 and "list.txt:2" in capsys.readouterr().err


def test_missing_input_file_exits_2(tmp_path):
    assert cli.main(["--input", str(tmp_path / "nope.txt"), "--out", str(tmp_path / "o")]) == 2


def test_missing_obabel_exits_2(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr("ligprep.convert.shutil.which", lambda name: None)
    code, _ = run(tmp_path, ["aspirin"])
    assert code == 2 and "Open Babel" in capsys.readouterr().err
