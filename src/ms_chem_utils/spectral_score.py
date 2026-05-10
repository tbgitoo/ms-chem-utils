
from __future__ import annotations

import math
from typing import Any, Dict, Iterable, List, Mapping, Optional, Tuple, Union


def _safe_float(x: Any) -> Optional[float]:
    try:
        if x is None:
            return None
        return float(x)
    except Exception:
        return None


def _build_binned_spectrum(
    peaks: Iterable[Mapping[str, Any]],
    mz_key: str = "mz",
    intensity_key: str = "intensity",
    binning: str = "round",   # "round" only for now
    min_mz: float = 0.0,
) -> Dict[int, Dict[str, float]]:
    """
    Build integer-binned spectrum:
      bin -> { 'intensity': sum_intensity, 'mz_rep': intensity-weighted mean m/z }

    If multiple peaks fall into the same integer bin, intensities are summed and mz_rep
    becomes intensity-weighted mean.
    """
    bins: Dict[int, Dict[str, float]] = {}

    for p in peaks:
        mz = _safe_float(p.get(mz_key))
        inten = _safe_float(p.get(intensity_key))
        if mz is None or inten is None:
            continue
        if mz < min_mz or inten < 0:
            continue

        if binning == "round":
            b = int(round(mz))
        else:
            raise ValueError(f"Unsupported binning={binning}")

        if b not in bins:
            bins[b] = {"intensity": 0.0, "mz_num": 0.0, "mz_den": 0.0}

        bins[b]["intensity"] += inten
        bins[b]["mz_num"] += mz * inten
        bins[b]["mz_den"] += inten

    # finalize mz_rep
    for b, d in list(bins.items()):
        den = d.get("mz_den", 0.0)
        if den > 0:
            d["mz_rep"] = d["mz_num"] / den
        else:
            d["mz_rep"] = float(b)  # fallback
        # clean internals
        d.pop("mz_num", None)
        d.pop("mz_den", None)

    return bins


def _pearson_r(x: List[float], y: List[float]) -> Optional[float]:
    """
    Pearson correlation coefficient.
    Returns None if undefined (e.g., zero variance).
    """
    if len(x) != len(y) or len(x) < 2:
        return None

    mx = sum(x) / len(x)
    my = sum(y) / len(y)

    sx2 = 0.0
    sy2 = 0.0
    sxy = 0.0
    for xi, yi in zip(x, y):
        dx = xi - mx
        dy = yi - my
        sx2 += dx * dx
        sy2 += dy * dy
        sxy += dx * dy

    if sx2 <= 0.0 or sy2 <= 0.0:
        return None

    return sxy / math.sqrt(sx2 * sy2)


def score_spectral_match(
    experimental: List[Dict[str, Any]],
    consensus: Union[
        List[Dict[str, Any]],
        Dict[Tuple[str, str], List[Dict[str, Any]]]
    ],
    *,
    # consensus selection
    mode: str = "LC",
    polarity: str = "positive",
    # keying
    exp_mz_key: str = "mz",
    exp_intensity_key: str = "intensity",
    cons_mz_key: str = "mz",
    cons_intensity_key: str = "mean_intensity",
    # scoring params
    intensity_power: float = 0.5,      # sqrt-like compress; set 1.0 for raw
    ppm_scale: float = 30.0,           # smaller = stricter mass accuracy
    w_intensity: float = 0.75,
    w_ppm: float = 0.75,
    min_mz: float = 0.0,

    # control how much "agreeing zeros" contribute (log-scaled)
    zero_agreement_weight: float = 0,  # 0..1, typical small (0.05..0.25)
) -> Dict[str, Any]:
    """
    Score match between an experimental spectrum and a consensus spectrum.

    Intensity scoring:
      - Bin by integer m/z (round).
      - Build vectors from lowest bin in either spectrum to highest bin in either spectrum.
      - Missing bins assumed 0 (so extra peaks are penalized).
      - Pearson correlation on normalized/compressed intensities.
      - Plus a small log-scaled "zero-agreement" term to reflect shared absences
        without letting large zero regions dominate.

    Mass scoring:
      - Mean absolute ppm error across overlapping bins (both nonzero).
      - Convert to score via exp(-mean_abs_ppm / ppm_scale).

    Returns a structured dict with diagnostics and 'ok' flag.
    """

    # 1) Extract the consensus list if consensus is keyed by (mode, polarity)
    cons_list: Optional[List[Dict[str, Any]]] = None
    if isinstance(consensus, dict):
        cons_list = consensus.get((mode, polarity))
    else:
        cons_list = consensus

    if not cons_list:
        return {
            "ok": False,
            "reason": f"No consensus spectrum available for {(mode, polarity)}",
            "score": 0.0,
            "intensity_score": 0.0,
            "mass_score": 0.0,
            "pearson_r": None,
            "mean_abs_ppm": None,
            "n_overlap_bins": 0,
            "mode": mode,
            "polarity": polarity,
        }

    # 2) Bin both spectra
    exp_bins = _build_binned_spectrum(
        experimental,
        mz_key=exp_mz_key,
        intensity_key=exp_intensity_key,
        binning="round",
        min_mz=min_mz,
    )
    cons_bins = _build_binned_spectrum(
        cons_list,
        mz_key=cons_mz_key,
        intensity_key=cons_intensity_key,
        binning="round",
        min_mz=min_mz,
    )

    if not exp_bins or not cons_bins:
        return {
            "ok": False,
            "reason": "Empty binned spectrum (experimental or consensus)",
            "score": 0.0,
            "intensity_score": 0.0,
            "mass_score": 0.0,
            "pearson_r": None,
            "mean_abs_ppm": None,
            "n_overlap_bins": 0,
            "mode": mode,
            "polarity": polarity,
        }

    # Optional: treat spectra with zero total intensity as empty
    exp_total = sum(d.get("intensity", 0.0) for d in exp_bins.values())
    cons_total = sum(d.get("intensity", 0.0) for d in cons_bins.values())
    if exp_total <= 0.0 or cons_total <= 0.0:
        return {
            "ok": False,
            "reason": "Zero total intensity (experimental or consensus)",
            "score": 0.0,
            "intensity_score": 0.0,
            "mass_score": 0.0,
            "pearson_r": None,
            "mean_abs_ppm": None,
            "n_overlap_bins": 0,
            "mode": mode,
            "polarity": polarity,
        }

    # 3) Determine bin range: from lowest in either to highest in either
    min_bin = min(min(exp_bins.keys()), min(cons_bins.keys()))
    max_bin = max(max(exp_bins.keys()), max(cons_bins.keys()))
    if max_bin < min_bin:
        max_bin = min_bin

    n_bins_total = (max_bin - min_bin + 1)

    # 4) Build intensity vectors with zeros for missing bins across [min_bin..max_bin]
    x_raw = [0.0] * n_bins_total
    y_raw = [0.0] * n_bins_total

    for b, d in exp_bins.items():
        if min_bin <= b <= max_bin:
            x_raw[b - min_bin] = d["intensity"]
    for b, d in cons_bins.items():
        if min_bin <= b <= max_bin:
            y_raw[b - min_bin] = d["intensity"]

    # 5) Normalize to base peak to focus on shape
    max_x = max(x_raw) if x_raw else 0.0
    max_y = max(y_raw) if y_raw else 0.0

    x = x_raw[:]  # work copies
    y = y_raw[:]

    if max_x > 0:
        x = [v / max_x for v in x]
    if max_y > 0:
        y = [v / max_y for v in y]

    # 6) Optional intensity compression (reduce dominance of base peak)
    if intensity_power != 1.0:
        x = [v ** intensity_power for v in x]
        y = [v ** intensity_power for v in y]

    # 7) Pearson correlation as intensity similarity (clamp negatives to 0)
    r = _pearson_r(x, y)
    if r is None:
        pearson_component = 0.0
    else:
        pearson_component = max(0.0, min(1.0, r))

    # 7b) Log-scaled "agreeing zeros" term:
    # Count bins where BOTH are zero (after normalization/compression still zero iff raw was zero).
    n_zero_zero = 0
    for xv, yv in zip(x_raw, y_raw):
        if xv == 0.0 and yv == 0.0:
            n_zero_zero += 1

    # Score in [0,1], but grows slowly as zeros accumulate
    # If n_bins_total is 1, avoid divide-by-zero by forcing denom >= log1p(1)=~0.693
    denom = math.log1p(max(n_bins_total, 1))
    zero_agreement_score = (math.log1p(n_zero_zero) / denom) if denom > 0 else 0.0
    zero_agreement_score = max(0.0, min(1.0, zero_agreement_score))

    # Combine intensity score components
    zaw = max(0.0, min(1.0, float(zero_agreement_weight)))
    intensity_score = (1.0 - zaw) * pearson_component + zaw * zero_agreement_score

    # 8) Mass accuracy score over overlapping bins (both nonzero in the bin)
    overlaps = []
    abs_ppms = []

    for b in range(min_bin, max_bin + 1):
        if b in exp_bins and b in cons_bins:
            if exp_bins[b]["intensity"] > 0 and cons_bins[b]["intensity"] > 0:
                mz_e = exp_bins[b]["mz_rep"]
                mz_c = cons_bins[b]["mz_rep"]
                if mz_c and mz_c != 0:
                    ppm = (mz_e - mz_c) / mz_c * 1e6
                    overlaps.append({
                        "bin": b,
                        "mz_exp": mz_e,
                        "mz_cons": mz_c,
                        "ppm": ppm,
                        # provide normalized intensities at that bin (for debugging)
                        "int_exp_norm": x[b - min_bin] if (0 <= (b - min_bin) < n_bins_total) else None,
                        "int_cons_norm": y[b - min_bin] if (0 <= (b - min_bin) < n_bins_total) else None,
                    })
                    abs_ppms.append(abs(ppm))

    n_overlap = len(abs_ppms)
    if n_overlap == 0:
        mean_abs_ppm = None
        mass_score = 0.0
    else:
        mean_abs_ppm = sum(abs_ppms) / n_overlap
        mass_score = math.exp(-mean_abs_ppm / max(ppm_scale, 1e-9))

    # 9) Combine scores with normalized weights
    

    combined = (w_intensity * intensity_score) + (w_ppm* mass_score)

    return {
        "ok": True,
        "score": float(combined),
        "intensity_score": float(intensity_score),
        "mass_score": float(mass_score),

        # intensity diagnostics
        "pearson_r": None if r is None else float(r),
        "pearson_component": float(pearson_component),
        "zero_agreement_score": float(zero_agreement_score),
        "n_zero_zero_bins": int(n_zero_zero),
        "zero_agreement_weight": float(zaw),

        # mass diagnostics
        "mean_abs_ppm": None if mean_abs_ppm is None else float(mean_abs_ppm),
        "n_overlap_bins": int(n_overlap),
        "overlaps": overlaps,

        # range diagnostics
        "min_bin_used": int(min_bin),
        "max_bin_used": int(max_bin),
        "n_bins_total": int(n_bins_total),

        "mode": mode,
        "polarity": polarity,
        "params": {
            "intensity_power": float(intensity_power),
            "ppm_scale": float(ppm_scale),
            "w_intensity": float(w_intensity),
            "w_ppm": float(w_ppm),
            "zero_agreement_weight": float(zaw),
        },
    }


    
    

