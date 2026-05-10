from __future__ import annotations

import time
from typing import Any, Dict, List, Optional
import requests

PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest"

class PubChemError(RuntimeError):
    pass

def _http_get_json(url: str, params: Optional[dict] = None, timeout: int = 20) -> Dict[str, Any]:
    r = requests.get(url, params=params, timeout=timeout, headers={"Accept": "application/json"})
    if r.status_code != 200:
        msg = r.text.strip()
        raise PubChemError(f"HTTP {r.status_code} from PubChem. Response: {msg[:500]}")
    try:
        return r.json()
    except ValueError:
        raise PubChemError(f"Non-JSON response from PubChem: {r.text[:500]}")

def pubchem_cids_from_smiles(smiles: str, timeout: int = 20, retries: int = 2, backoff_s: float = 0.6) -> List[int]:
    url = f"{PUBCHEM_BASE}/pug/compound/smiles/cids/JSON"
    params = {"smiles": (smiles or "").strip()}

    last_err: Optional[Exception] = None
    for attempt in range(retries + 1):
        try:
            data = _http_get_json(url, params=params, timeout=timeout)
            cids = data.get("IdentifierList", {}).get("CID", [])
            return [int(x) for x in cids]
        except PubChemError as e:
            last_err = e
            time.sleep(backoff_s * (2 ** attempt))
    raise last_err  # type: ignore

def pubchem_properties_from_cid(cid: int, properties: Optional[List[str]] = None, timeout: int = 20) -> Dict[str, Any]:
    if properties is None:
        properties = ["Title", "MolecularFormula", "MolecularWeight", "IsomericSMILES", "InChIKey", "XLogP"]
    prop_str = ",".join(properties)
    url = f"{PUBCHEM_BASE}/pug/compound/cid/{cid}/property/{prop_str}/JSON"
    data = _http_get_json(url, timeout=timeout)
    props = data.get("PropertyTable", {}).get("Properties", [])
    return props[0] if props else {}

def pubchem_lookup_smiles(smiles: str) -> Dict[str, Any]:
    smiles = (smiles or "").strip()
    if not smiles:
        return {"ok": False, "query": smiles, "reason": "empty SMILES string"}

    try:
        cids = pubchem_cids_from_smiles(smiles)
        if not cids:
            return {"ok": True, "query": smiles, "cids": [], "cid": None, "name": None, "properties": {}, "reason": "no CID found"}

        cid0 = cids[0]
        props = pubchem_properties_from_cid(cid0)
        name = props.get("Title")
        return {"ok": True, "query": smiles, "cids": cids, "cid": cid0, "name": name, "properties": props}
    except PubChemError as e:
        return {"ok": False, "query": smiles, "reason": str(e)}
    except Exception as e:
        return {"ok": False, "query": smiles, "reason": f"unexpected error: {e}"}