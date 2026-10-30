"""Dataset-level statistics of the per-site ambiguity a_i.

Headline (since 2026-10-28):
    T_p = mean of the upper-p mass of the a_i distribution, p = 0.08  ("T8")
          exact fractional tail: n*p observations, with a fractional weight on the boundary one.
Also reported for continuity: T20, A = mean a_i, median (sigma), C = A / median,
P = fraction of sites with a_i > tau (tau = 0.5 muB), k_eff statistics, and the legacy
(sigma, p) verdict.  Force errors never enter any statistic.

`tail` is the function of the source project (scan_tail_descriptors_20261011.py) verbatim.
"""

from __future__ import annotations

from typing import Iterable, Optional

import numpy as np

P_FINE, P_AFFORDABLE, SIGMA_SEPARABLE = 0.01, 0.15, 0.01


def tail(a, q: float) -> tuple[float, float]:
    """Exact upper-q mass with a fractional observation at the boundary.

    Returns (H_q, T_q): H_q = share of the total sum carried by the upper-q mass,
    T_q = mean ambiguity of that mass (T_q = A * H_q / q).
    """
    x = np.sort(np.asarray(a, float))[::-1]
    mass = len(x) * q
    k = int(np.floor(mass))
    total = x[:k].sum() + ((mass - k) * x[k] if k < len(x) else 0.0)
    return total / x.sum(), total / mass


def tail_mean(a, q: float) -> float:
    """T_q: mean of the upper-q mass of a (fractional boundary observation).  T_0.08 is 'T8'."""
    x = np.asarray(a, float)
    if x.size == 0:
        return float("nan")
    return tail(x, q)[1]


def tail_flag(a, q: float) -> tuple[np.ndarray, float]:
    """Integer-count tail membership: n_tail = ceil(q n), threshold = n_tail-th largest a_i.
    Used to mark the top-p sites in figures; T_p itself uses the fractional definition."""
    a = np.asarray(a, float)
    n_tail = int(np.ceil(len(a) * q))
    thr = np.sort(a)[::-1][n_tail - 1]
    return a >= thr, float(thr)


def frame_force_mask(frame_fmax: np.ndarray, fmax: Optional[float]) -> np.ndarray:
    """Query-site mask of the frame-force filter: keep sites whose frame has max|F_DFT| <= fmax.
    Frames without a force record (NaN) are kept and should be reported.  None -> keep all."""
    frame_fmax = np.asarray(frame_fmax, float)
    if fmax is None:
        return np.ones(len(frame_fmax), bool)
    return ~(frame_fmax > fmax)


def legacy_verdict(p_other_state: float, sigma: float) -> str:
    """(sigma, p) rule of 2026-10-02: p < 1 % fine; p <= 15 % and sigma <= 0.01 muB down-weight; else magnetic model."""
    if p_other_state < P_FINE:
        return "dataset is fine"
    if p_other_state <= P_AFFORDABLE and sigma <= SIGMA_SEPARABLE:
        return "down-weight with d_i"
    return "needs a magnetic model"


def dataset_statistics(a, k_eff=None, *, tail_fraction: float = 0.08,
                       extra_tail_fractions: Iterable[float] = (0.20,), tau: float = 0.5,
                       k_cap: Optional[int] = None) -> dict:
    """One row of dataset statistics from the (already filtered) per-site ambiguities."""
    a = np.asarray(a, float)
    if a.size == 0 or not np.isfinite(a).all() or (a < 0).any():
        raise ValueError("ambiguities must be a non-empty, finite, non-negative array")
    median = float(np.median(a))
    p = float(np.mean(a > tau))
    row = {
        "n_sites": int(a.size),
        f"T{round(tail_fraction * 100):d}_muB": tail_mean(a, tail_fraction),
        "tail_fraction": float(tail_fraction),
        "A_muB": float(a.mean()),
        "median_muB": median,
        "C": float(a.mean() / max(median, 1e-9)),
        f"P_percent_a_gt_{tau:g}": 100.0 * p,
        "p_other_state": p,
        "sigma_muB": median,
        "frac_zero": float(np.mean(a == 0)),
        "max_muB": float(a.max()),
        "legacy_verdict_sigma_p": legacy_verdict(p, median),
    }
    for q in extra_tail_fractions:
        row[f"T{round(q * 100):d}_muB"] = tail_mean(a, q)
    if k_eff is not None:
        k_eff = np.asarray(k_eff)
        row.update({"k_eff_median": float(np.median(k_eff)), "k_eff_mean": float(k_eff.mean()),
                    "k_eff_min": int(k_eff.min()), "k_eff_max": int(k_eff.max())})
        if k_cap is not None:
            row["k_eff_at_cap_fraction"] = float(np.mean(k_eff == k_cap))
    return row
