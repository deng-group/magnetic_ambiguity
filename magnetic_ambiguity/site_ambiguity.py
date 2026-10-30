"""Per-site ambiguity from the stored neighbour lists.

    a_i = median_{j in twins(i)} |m_i - m_j|          signed moments (default; the T8 definition)
    a_i^u = median_{j in twins(i)} | |m_i| - |m_j| |  unsigned (decision 2026-10-01 branch; kept as a column)

Moments are the final collinear DFT moments in muB (float64).  Both values are always computed and
stored; `moments` in the settings selects which one is the headline column `ambiguity_muB`.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def prefix_median(x: np.ndarray, counts: np.ndarray) -> np.ndarray:
    """median of x[i, :counts[i]] for every row (vectorised over equal counts)."""
    counts = np.asarray(counts, int)
    out = np.empty(x.shape[0], dtype=np.float64)
    for c in np.unique(counts):
        rows = counts == c
        out[rows] = np.median(x[rows, :c], axis=1)
    return out


def site_ambiguity(moments: np.ndarray, nb_idx: np.ndarray, k_eff: np.ndarray, unsigned: bool = False) -> np.ndarray:
    m = np.abs(np.asarray(moments, float)) if unsigned else np.asarray(moments, float)
    delta = np.abs(m[nb_idx] - m[:, None])
    return prefix_median(delta, k_eff)


def site_table(sites: pd.DataFrame, nb_idx: np.ndarray, nb_dist: np.ndarray, k_eff: np.ndarray,
               headline: str = "signed") -> pd.DataFrame:
    """Attach ambiguity columns to the site table (one row per central-element site, source order)."""
    m = sites["magmom_muB"].to_numpy(float)
    frames = sites["frame_id"].to_numpy()
    out = sites.copy()
    a_signed = site_ambiguity(m, nb_idx, k_eff, unsigned=False)
    a_unsigned = site_ambiguity(m, nb_idx, k_eff, unsigned=True)
    out["adaptive_k"] = np.asarray(k_eff, np.int16)
    out["ambiguity_signed_muB"] = a_signed
    out["ambiguity_unsigned_muB"] = a_unsigned
    out["ambiguity_muB"] = a_signed if headline == "signed" else a_unsigned
    out["ambiguity_exact_k_signed_muB"] = np.median(np.abs(m[nb_idx] - m[:, None]), axis=1)
    out["neighbor_median_magmom_muB"] = prefix_median(m[nb_idx], k_eff)
    out["nearest_neighbor_soap_distance"] = nb_dist[:, 0]
    out["band_radius_soap_distance"] = np.float32(1.0) * nb_dist[np.arange(len(k_eff)), np.asarray(k_eff, int) - 1]
    same = frames[nb_idx] == frames[:, None]
    out["n_same_frame_twins"] = np.asarray([same[i, :k].sum() for i, k in enumerate(k_eff)], np.int16)
    out["n_same_frame_exact_k"] = same.sum(axis=1).astype(np.int16)
    return out
