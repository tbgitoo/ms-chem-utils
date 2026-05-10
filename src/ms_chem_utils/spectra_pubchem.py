
from typing import List, Dict, Any


def spectra_from_smiles(smiles: str, timeout: int = 20) -> List[Dict[str, Any]]:
    """
    Given a SMILES, return available MS spectra metadata from PubChem (pug_view).

    This returns *spectral descriptors* (LC-MS / GC-MS), not full peak lists.
    Full spectra must be retrieved from the referenced source (e.g. MassBank).
    """

    
    smiles = (smiles or "").strip()
    if not smiles:
        return []

    # ---- 1) SMILES -> CID (lightweight, reliable) ----
    try:
        r = requests.get(
            f"{PUBCHEM_BASE}/pug/compound/smiles/cids/JSON",
            params={"smiles": smiles},
            timeout=timeout,
        )
        r.raise_for_status()
    except Exception:
        return []

    cids = r.json().get("IdentifierList", {}).get("CID", [])
    if not cids:
        return []

    cid = cids[0]

    # ---- 2) CID -> pug_view record (spectral metadata) ----
    try:
        r = requests.get(
            f"{PUBCHEM_BASE}/pug_view/data/compound/{cid}/JSON",
            timeout=timeout,
        )
        r.raise_for_status()
    except Exception:
        return []

    data = r.json()
    spectra: List[Dict[str, Any]] = []

    def walk(section):
        if not isinstance(section, dict):
            return

        toc = section.get("TOCHeading", "")
        if toc in {"LC-MS", "GC-MS", "Mass Spectrometry"}:
            # Group Information blocks by ReferenceNumber (= one spectrum)
            by_ref = {}
            for info in section.get("Information", []):
                ref = info.get("ReferenceNumber")
                if ref is None:
                    continue
                by_ref.setdefault(ref, []).append(info)

            for ref, items in by_ref.items():
                spec = {
                    "source": None,
                    "accession": None,
                    "ms_level": None,
                    "ionization": None,
                    "polarity": None,
                    "instrument_type": None,
                    "collision_energy": None,
                    "precursor_adduct": None,
                    "top_peaks": [],
                    "splash": None,
                    "external_url": None,
                }

                for it in items:
                    name = it.get("Name")
                    val = it.get("Value", {})
                    url = it.get("URL")

                    def sval():
                        swm = val.get("StringWithMarkup", [])
                        return swm[0]["String"] if swm else None

                    if name == "Accession ID":
                        spec["accession"] = sval()
                        spec["external_url"] = url
                        spec["source"] = "MassBank"
                    elif name == "MS Level":
                        spec["ms_level"] = sval()
                    elif name == "Ionization":
                        spec["ionization"] = sval()
                    elif name == "Ionization Mode":
                        spec["polarity"] = sval()
                    elif name == "Instrument Type":
                        spec["instrument_type"] = sval()
                    elif name == "Collision Energy":
                        spec["collision_energy"] = sval()
                    elif name == "Precursor Adduct":
                        spec["precursor_adduct"] = sval()
                    elif name == "Top 5 Peaks":
                        peaks = []
                        for p in val.get("StringWithMarkup", []):
                            try:
                                mz, inten = p["String"].split()
                                peaks.append((float(mz), float(inten)))
                            except Exception:
                                pass
                        spec["top_peaks"] = peaks
                    elif name == "SPLASH":
                        spec["splash"] = sval()

                spectra.append(spec)

        for s in section.get("Section", []):
            walk(s)

    for top in data.get("Record", {}).get("Section", []):
        walk(top)

    return spectra






def consensus_spectra_from_spectra(
    spectra: List[Dict[str, Any]],
    mz_round_digits: int = 0
) -> Dict[Tuple[str, str], List[Dict[str, float]]]:
    """
    Build consensus spectra aggregated by (GC/LC, polarity).

    - m/z values are grouped by rounded m/z (default: integer).
    - representative m/z = mean of contributing high-precision m/z values
    - mean intensity computed across *all spectra in group* (0 if absent)
    """

    # ---- Helper: GC vs LC classification ----
    def classify_gc_lc(instr: str | None) -> str:
        if not instr:
            return "UNKNOWN"
        s = instr.upper()
        if "EI" in s or "GC" in s:
            return "GC"
        if "LC" in s or "ESI" in s:
            return "LC"
        return "UNKNOWN"

    # ---- Group spectra ----
    groups: Dict[Tuple[str, str], List[Dict[str, Any]]] = defaultdict(list)

    for spec in spectra:
        if not spec.get("top_peaks"):
            continue

        mode = classify_gc_lc(spec.get("instrument_type"))
        polarity = (spec.get("polarity") or "UNKNOWN").lower()

        groups[(mode, polarity)].append(spec)

    # ---- Build consensus per group ----
    consensus = {}

    for key, specs in groups.items():
        # Collect all rounded m/z values
        mz_bins: Dict[float, List[Tuple[float, float]]] = defaultdict(list)

        for spec in specs:
            for mz, inten in spec["top_peaks"]:
                mz_round = round(mz, mz_round_digits)
                mz_bins[mz_round].append((mz, inten))

        # Compute mean intensity with zero fill
        peaks_out = []

        for mz_round, _ in mz_bins.items():
            mz_values = []
            intensities = []

            for spec in specs:
                found = False
                for mz, inten in spec["top_peaks"]:
                    if round(mz, mz_round_digits) == mz_round:
                        mz_values.append(mz)
                        intensities.append(inten)
                        found = True
                        break
                if not found:
                    intensities.append(0.0)

            peaks_out.append({
                "mz": sum(mz_values) / len(mz_values),
                "mz_round": mz_round,
                "mean_intensity": sum(intensities) / len(intensities)
            })

        # Sort by mean intensity descending
        peaks_out.sort(key=lambda x: x["mean_intensity"], reverse=True)
        consensus[key] = peaks_out

    return consensus







