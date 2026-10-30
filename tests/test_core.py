"""Unit tests that need only numpy/pandas (no SOAP, no data): run with `python -m pytest tests`
or `python tests/test_core.py`."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from magnetic_ambiguity.dataset_ambiguity import (dataset_statistics, frame_force_mask, legacy_verdict,  # noqa: E402
                                                  tail, tail_flag, tail_mean)
from magnetic_ambiguity.neighbors import adaptive_k, exact_knn, fixed_k  # noqa: E402
from magnetic_ambiguity.settings import Settings, NeighborSettings, corpus_settings, inhouse_settings  # noqa: E402
from magnetic_ambiguity.site_ambiguity import prefix_median, site_ambiguity, site_table  # noqa: E402
from magnetic_ambiguity.soap import BlockSparseDescriptors, descriptors_from_dense  # noqa: E402


def test_tail_analytic():
    # definition checks of the source project (scan_tail_descriptors_20261011.py)
    assert np.allclose(tail(np.ones(10), .15), (.15, 1.))
    assert np.allclose(tail(np.r_[10., np.zeros(9)], .10), (1., 10.))
    a = np.arange(1, 101, dtype=float)              # 1..100: upper 8 % = 8 observations 100..93
    assert np.isclose(tail_mean(a, .08), np.mean(np.arange(93, 101)))
    # fractional boundary: n = 10, q = 0.25 -> 2.5 observations
    a = np.array([5., 4., 3., 2., 1., 0, 0, 0, 0, 0])
    assert np.isclose(tail(a, .25)[1], (5 + 4 + 0.5 * 3) / 2.5)
    assert np.isclose(tail(a, .25)[1], a.mean() * tail(a, .25)[0] / .25)     # T_q = A H_q / q
    flag, thr = tail_flag(np.arange(1, 101, dtype=float), .08)
    assert flag.sum() == 8 and thr == 93


def test_adaptive_rule():
    d = np.array([[0.0, 0.1, 0.2, 0.3, 0.4, 0.45, 0.49, 0.5, 0.51, 1.0]], dtype=np.float32)
    k = adaptive_k(d, alpha=1.25, m_ref=5)          # radius 1.25 * 0.4 = 0.5 -> 8 neighbours (<= 0.5)
    assert k.tolist() == [8]
    assert adaptive_k(d, alpha=1.0, m_ref=5).tolist() == [5]       # minimum m
    assert adaptive_k(d, alpha=10.0, m_ref=5).tolist() == [10]     # cap K
    assert fixed_k(d, 3).tolist() == [3]


def test_prefix_median_and_site_ambiguity():
    m = np.array([1.0, -1.0, 1.1, 3.0, 0.0])
    nb = np.array([[1, 2, 3, 4], [0, 2, 3, 4], [0, 1, 3, 4], [0, 1, 2, 4], [0, 1, 2, 3]])
    keff = np.array([2, 4, 3, 2, 4])
    a = site_ambiguity(m, nb, keff)                 # signed
    assert np.isclose(a[0], np.median([2.0, 0.1]))
    assert np.isclose(a[1], np.median([2.0, 2.1, 4.0, 1.0]))
    u = site_ambiguity(m, nb, keff, unsigned=True)
    assert np.isclose(u[0], np.median([0.0, 0.1]))
    assert np.allclose(prefix_median(np.abs(m[nb] - m[:, None]), keff), a)


def test_exact_knn_matches_brute_force_and_excludes_self():
    rng = np.random.default_rng(0)
    x = rng.normal(size=(60, 7)).astype(np.float32)
    x /= np.linalg.norm(x, axis=1, keepdims=True)
    frames = np.repeat(np.arange(12), 5)
    desc = descriptors_from_dense(x)
    nb, dist = exact_knn(desc, frames, k=10, query_chunk=17)
    full = np.sqrt(np.clip(2 - 2 * (x @ x.T), 0, None))
    for i in range(60):
        assert i not in nb[i]
        others = np.delete(np.arange(60), i)
        expect = others[np.argsort(full[i, others], kind="stable")[:10]]
        assert np.array_equal(nb[i], expect)
        assert np.allclose(dist[i], full[i, expect], atol=2e-6)
    nb2, _ = exact_knn(desc, frames, k=10, exclude_same_frame=True)
    assert not np.any(frames[nb2] == frames[:, None])


def test_block_sparse_dot_equals_dense():
    """Two species bases with a shared block: the block-sparse dot equals the zero-padded dense dot."""
    rng = np.random.default_rng(1)
    n = 8
    xa = rng.normal(size=(n, 3)).astype(np.float32)            # block A|A (all sites)
    xb = rng.normal(size=(4, 2)).astype(np.float32)            # block A|B (sites 0-3 only)
    dense = np.zeros((n, 5), np.float32); dense[:, :3] = xa; dense[:4, 3:] = xb
    dense /= np.linalg.norm(dense, axis=1, keepdims=True)
    desc = BlockSparseDescriptors(n_sites=n, basis_of_site=["A|B"] * 4 + ["A"] * 4)
    desc.blocks["A|A"] = (np.arange(n), dense[:, :3].copy())
    desc.blocks["A|B"] = (np.arange(4), dense[:4, 3:].copy())
    out = np.empty((n, n), np.float32)
    assert np.allclose(desc.dots(0, n, out), dense @ dense.T, atol=1e-6)


def test_dataset_statistics_and_filters():
    a = np.r_[np.zeros(90), np.full(10, 2.0)]
    row = dataset_statistics(a, k_eff=np.full(100, 7), tail_fraction=0.08, tau=0.5, k_cap=100)
    assert np.isclose(row["T8_muB"], 2.0) and np.isclose(row["T20_muB"], 1.0)
    assert np.isclose(row["P_percent_a_gt_0.5"], 10.0) and row["k_eff_median"] == 7
    assert legacy_verdict(0.001, 0.1) == "dataset is fine"
    assert legacy_verdict(0.06, 0.004) == "down-weight with d_i"
    assert legacy_verdict(0.3, 0.1) == "needs a magnetic model"
    fm = np.array([1.0, 25.0, np.nan, 20.0])
    assert frame_force_mask(fm, 20.0).tolist() == [True, False, True, True]
    assert frame_force_mask(fm, None).all()


def test_site_table_same_frame_counts():
    sites = pd.DataFrame({"site_id": range(4), "frame_id": [0, 0, 1, 1], "atom_index": [0, 1, 0, 1],
                          "element": "Co", "magmom_muB": [1.0, 1.0, -1.0, 2.0]})
    nb = np.array([[1, 2, 3], [0, 3, 2], [3, 0, 1], [2, 1, 0]])
    dist = np.tile(np.array([0.1, 0.2, 0.3], np.float32), (4, 1))
    t = site_table(sites, nb, dist, np.array([2, 2, 3, 1]))
    assert t["n_same_frame_twins"].tolist() == [1, 1, 1, 1]
    assert t["n_same_frame_exact_k"].tolist() == [1, 1, 1, 1]
    assert np.isclose(t.loc[2, "ambiguity_signed_muB"], np.median([3.0, 2.0, 2.0]))
    assert np.isclose(t.loc[2, "ambiguity_unsigned_muB"], np.median([1.0, 0.0, 0.0]))


def test_settings_version_is_stable_and_sensitive():
    s = Settings()
    assert len(s.version()) == 16 and s.version() == Settings().version()
    assert s.version() != s.with_(neighbors=NeighborSettings(k_exact=20)).version()
    assert inhouse_settings().version() != corpus_settings().version()          # species basis / |m| filter
    # the dataset section (tail fraction, force filter) must NOT change the cache version
    from magnetic_ambiguity.settings import DatasetSettings
    assert s.version() == s.with_(dataset=DatasetSettings(tail_fraction=0.2)).version()


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn(); print(f"ok  {name}")
