"""Minimal PubChem PUG REST client: resolve an entry to a CID, fetch properties and a 3D SDF."""
from __future__ import annotations

import time

import requests

from .errors import LigprepError
from .inputs import Entry

BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
USER_AGENT = "ligprep-cli/0.1 (+https://github.com/Akkaya-Furkan/ligprep-cli)"
RETRY_STATUS = {429, 500, 502, 503, 504}


class PubChem:
    def __init__(self, session: requests.Session | None = None, delay: float = 0.25,
                 retries: int = 3, timeout: float = 30.0):
        self.session = session or requests.Session()
        self.session.headers["User-Agent"] = USER_AGENT  # requests would otherwise send python-requests/x
        self.delay, self.retries, self.timeout = delay, retries, timeout
        self._last = 0.0

    def _request(self, method: str, url: str, **kw) -> requests.Response:
        for attempt in range(self.retries + 1):
            wait = self.delay - (time.monotonic() - self._last)  # PubChem asks for <= 5 requests/s
            if wait > 0:
                time.sleep(wait)
            self._last = time.monotonic()
            try:
                r = self.session.request(method, url, timeout=self.timeout, **kw)
            except requests.RequestException as e:
                if attempt == self.retries:
                    raise LigprepError(f"network error: {e}") from None
                time.sleep(2 ** attempt)
                continue
            if r.status_code in RETRY_STATUS and attempt < self.retries:
                time.sleep(2 ** attempt)
                continue
            return r
        raise LigprepError("unreachable")

    def cid_for(self, entry: Entry) -> int:
        if entry.kind == "cid":
            return int(entry.value)
        if entry.kind == "name":
            r = self._request("GET", f"{BASE}/compound/name/{requests.utils.quote(entry.value, safe='')}/cids/TXT")
        else:  # SMILES goes in the body: it contains characters that break URL paths
            r = self._request("POST", f"{BASE}/compound/smiles/cids/TXT", data={"smiles": entry.value})
        if r.status_code == 404:
            raise LigprepError(f"PubChem has no compound for {entry.kind} {entry.value!r}")
        if not r.ok:
            raise LigprepError(f"PubChem lookup failed (HTTP {r.status_code})")
        cids = [c for c in r.text.split() if c.isdigit() and c != "0"]
        if not cids:
            raise LigprepError(f"PubChem returned no CID for {entry.value!r}")
        return int(cids[0])

    def properties(self, cid: int) -> dict:
        r = self._request("GET", f"{BASE}/compound/cid/{cid}/property/MolecularFormula,SMILES/JSON")
        if not r.ok:
            raise LigprepError(f"could not read properties of CID {cid} (HTTP {r.status_code})")
        try:
            return r.json()["PropertyTable"]["Properties"][0]
        except (KeyError, IndexError, ValueError):
            raise LigprepError(f"unexpected property response for CID {cid}") from None

    def sdf_3d(self, cid: int) -> str:
        r = self._request("GET", f"{BASE}/compound/cid/{cid}/SDF", params={"record_type": "3d"})
        if r.status_code == 404:
            raise LigprepError(f"CID {cid} has no 3D conformer in PubChem")
        if not r.ok:
            raise LigprepError(f"could not download 3D SDF of CID {cid} (HTTP {r.status_code})")
        return r.text
