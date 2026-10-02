import shutil
from pathlib import Path

import pytest
import responses

from ligprep.pubchem import BASE

DATA = Path(__file__).parent / "data"
ASPIRIN_FORMULA = "C9H8O4"
ASPIRIN_SMILES = "CC(=O)OC1=CC=CC=C1C(=O)O"

needs_obabel = pytest.mark.skipif(shutil.which("obabel") is None, reason="Open Babel not installed")


@pytest.fixture
def aspirin_sdf() -> str:
    return (DATA / "aspirin_3d.sdf").read_text()


@pytest.fixture
def rsps():
    with responses.RequestsMock(assert_all_requests_are_fired=False) as r:
        yield r


@pytest.fixture
def mock_aspirin(rsps, aspirin_sdf):
    """PubChem answers for aspirin (CID 2244); no network is ever touched."""
    def _mock(formula: str = ASPIRIN_FORMULA, sdf: str | None = None):
        rsps.add(responses.GET, f"{BASE}/compound/name/aspirin/cids/TXT", body="2244\n")
        rsps.add(responses.POST, f"{BASE}/compound/smiles/cids/TXT", body="2244\n")
        rsps.add(responses.GET, f"{BASE}/compound/cid/2244/property/MolecularFormula,SMILES/JSON",
                 json={"PropertyTable": {"Properties": [
                     {"CID": 2244, "MolecularFormula": formula, "SMILES": ASPIRIN_SMILES}]}})
        rsps.add(responses.GET, f"{BASE}/compound/cid/2244/SDF", body=sdf or aspirin_sdf)
        return rsps
    return _mock


@pytest.fixture(autouse=True)
def no_sleep(monkeypatch):
    monkeypatch.setattr("ligprep.pubchem.time.sleep", lambda s: None)
