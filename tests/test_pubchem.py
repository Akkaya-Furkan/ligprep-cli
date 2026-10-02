import pytest
import requests
import responses

from ligprep.errors import LigprepError
from ligprep.inputs import Entry
from ligprep.pubchem import BASE, USER_AGENT, PubChem

NAME_URL = f"{BASE}/compound/name/aspirin/cids/TXT"


def client(**kw):
    return PubChem(delay=0, **kw)


def test_name_lookup(rsps):
    rsps.add(responses.GET, NAME_URL, body="2244\n")
    assert client().cid_for(Entry("name", "aspirin", "aspirin")) == 2244


def test_cid_needs_no_request(rsps):
    assert client().cid_for(Entry("cid", "3672", "3672")) == 3672
    assert len(rsps.calls) == 0


def test_smiles_goes_in_the_post_body(rsps):
    rsps.add(responses.POST, f"{BASE}/compound/smiles/cids/TXT", body="2244\n")
    assert client().cid_for(Entry("smiles", "CC(=O)O", "smiles:CC(=O)O")) == 2244
    assert rsps.calls[0].request.body == "smiles=CC%28%3DO%29O"


def test_unknown_name_is_a_ligprep_error(rsps):
    rsps.add(responses.GET, f"{BASE}/compound/name/zzz/cids/TXT", status=404)
    with pytest.raises(LigprepError, match="no compound"):
        client().cid_for(Entry("name", "zzz", "zzz"))


def test_cid_zero_means_not_found(rsps):
    rsps.add(responses.GET, NAME_URL, body="0\n")
    with pytest.raises(LigprepError):
        client().cid_for(Entry("name", "aspirin", "aspirin"))


def test_retries_a_503_then_succeeds(rsps):
    rsps.add(responses.GET, NAME_URL, status=503)
    rsps.add(responses.GET, NAME_URL, body="2244\n")
    assert client().cid_for(Entry("name", "aspirin", "aspirin")) == 2244
    assert len(rsps.calls) == 2


def test_gives_up_after_the_retry_budget(rsps):
    rsps.add(responses.GET, NAME_URL, status=503)
    with pytest.raises(LigprepError, match="HTTP 503"):
        client(retries=2).cid_for(Entry("name", "aspirin", "aspirin"))
    assert len(rsps.calls) == 3


def test_network_error_becomes_a_ligprep_error(rsps):
    rsps.add(responses.GET, NAME_URL, body=requests.ConnectionError("down"))
    with pytest.raises(LigprepError, match="network error"):
        client(retries=1).cid_for(Entry("name", "aspirin", "aspirin"))


def test_properties(mock_aspirin):
    mock_aspirin()
    props = client().properties(2244)
    assert props["MolecularFormula"] == "C9H8O4"


def test_malformed_property_response(rsps):
    rsps.add(responses.GET, f"{BASE}/compound/cid/1/property/MolecularFormula,SMILES/JSON", json={})
    with pytest.raises(LigprepError, match="unexpected"):
        client().properties(1)


def test_missing_3d_conformer_has_a_clear_message(rsps):
    rsps.add(responses.GET, f"{BASE}/compound/cid/5462224/SDF", status=404)
    with pytest.raises(LigprepError, match="no 3D conformer"):
        client().sdf_3d(5462224)


def test_sdf_request_asks_for_3d_and_sends_a_user_agent(mock_aspirin):
    rsps = mock_aspirin()
    client().sdf_3d(2244)
    req = rsps.calls[0].request
    assert "record_type=3d" in req.url
    assert req.headers["User-Agent"] == USER_AGENT
