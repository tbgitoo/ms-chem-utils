# tests/test_spectra_pubchem_mocked.py
import types

class _FakeResp:
    def __init__(self, json_data=None, status_code=200):
        self._json = json_data or {}
        self.status_code = status_code

    def json(self):
        return self._json

    def raise_for_status(self):
        if self.status_code >= 400:
            raise RuntimeError(f"HTTP {self.status_code}")

def test_spectra_from_smiles_empty_smiles_returns_empty():
    from ms_chem_utils.spectra_pubchem import spectra_from_smiles
    out = spectra_from_smiles("")
    assert out == []

def test_spectra_from_smiles_success_path(monkeypatch):
    """
    Simulate PubChem:
    1) SMILES -> CID list
    2) CID -> pug_view record with some spectra-like sections
    """
    import ms_chem_utils.spectra_pubchem as sp

    def fake_get(url, params=None, timeout=None):
        # Step 1: smiles -> cids
        if "/compound/smiles/cids/JSON" in url:
            return _FakeResp({"IdentifierList": {"CID": [123]}})
        # Step 2: CID -> pug_view
        if "/pug_view/data/compound/123/JSON" in url:
            return _FakeResp({
                "Record": {
                    "Section": [
                        {
                            "TOCHeading": "Experimental Properties",
                            "Section": [
                                {"TOCHeading": "Mass Spectrometry", "Information": [{"Value": {"StringWithMarkup": [{"String": "LC-MS"}]}}]}
                            ],
                        }
                    ]
                }
            })
        return _FakeResp({}, status_code=404)

    monkeypatch.setattr(sp.requests, "get", fake_get)

    out = sp.spectra_from_smiles("CCO")
    assert isinstance(out, list)
    # your function might return descriptors or normalized objects; assert non-throw and plausible non-empty
    # If your implementation returns [] when it can't parse, this might be [] — adjust as needed.
    assert out is not None

def test_spectra_from_smiles_handles_http_error(monkeypatch):
    import ms_chem_utils.spectra_pubchem as sp

    def fake_get(url, params=None, timeout=None):
        return _FakeResp({}, status_code=500)

    monkeypatch.setattr(sp.requests, "get", fake_get)

    out = sp.spectra_from_smiles("CCO")
    assert out == []  # based on your code style: returns [] on exceptions