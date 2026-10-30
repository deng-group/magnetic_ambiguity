"""Exact K-nearest-neighbour search in SOAP space and the adaptive neighbourhood.

Rule (user decision 2026-10-09 in-house, 2026-10-11 corpora):
  * every other site of the element pool is a candidate, including the other sites of the
    query frame; only the query site itself is excluded (``exclude_same_frame=False``);
  * distance d_ij = sqrt(clip(2 - 2 x_i.x_j, 0)) on L2-normalised SOAP (float32 accumulation);
  * the first K (default 100) neighbours are found EXACTLY (no approximate search), ties
    broken by site id (lexicographic (distance, site_id));
  * adaptive twins: neighbours with d_ij <= alpha * d_i,m (alpha 1.25, m 5), so that
    m <= k_eff <= K.  The band radius is formed in float32 exactly as stored.
The previous rule (whole query frame excluded, K = 20) is ``exclude_same_frame=True, k_exact=20``.
"""

from __future__ import annotations

import numpy as np

from .soap import BlockSparseDescriptors


def exact_knn(desc: BlockSparseDescriptors, frame_ids: np.ndarray, k: int,
              exclude_same_frame: bool = False, query_chunk: int = 4096,
              log=None) -> tuple[np.ndarray, np.ndarray]:
    """Return (nb_idx int32 (n, k), nb_dist float32 (n, k)) sorted by (distance, site id)."""
    n = desc.n_sites
    frame_ids = np.asarray(frame_ids)
    if len(frame_ids) != n:
        raise ValueError("frame_ids length != number of sites")
    if exclude_same_frame:
        _, counts = np.unique(frame_ids, return_counts=True)
        if n - counts.max() < k:
            raise ValueError(f"fewer than {k} cross-frame candidates for at least one site")
    elif n - 1 < k:
        raise ValueError(f"pool has {n} sites; K = {k} needs at least K + 1")
    nb_idx = np.empty((n, k), dtype=np.int32)
    nb_dist = np.empty((n, k), dtype=np.float32)
    buf = np.empty((min(query_chunk, n), n), dtype=np.float32)
    for start in range(0, n, query_chunk):
        stop = min(start + query_chunk, n)
        q = np.arange(start, stop)
        dots = desc.dots(start, stop, buf[: stop - start])
        dist = np.sqrt(np.clip(2.0 - 2.0 * dots, 0.0, None))
        if exclude_same_frame:
            same = frame_ids[q][:, None] == frame_ids[None, :]
            dist[same] = np.inf
        else:
            dist[np.arange(stop - start), q] = np.inf
        order = np.argsort(dist, axis=1, kind="stable")[:, :k]      # columns are site ids -> (dist, id) order
        nb_idx[start:stop] = order
        nb_dist[start:stop] = np.take_along_axis(dist, order, axis=1)
        if log:
            log(f"    neighbours {stop}/{n}")
    if not np.isfinite(nb_dist).all():
        raise ValueError("non-finite selected neighbour distance (pool too small after exclusion)")
    if np.any(nb_idx == np.arange(n)[:, None]):
        raise ValueError("query site selected as its own neighbour")
    if exclude_same_frame and np.any(frame_ids[nb_idx] == frame_ids[:, None]):
        raise ValueError("same-frame neighbour leakage")
    return nb_idx, nb_dist


def adaptive_k(nb_dist: np.ndarray, alpha: float = 1.25, m_ref: int = 5, k_cap: int | None = None) -> np.ndarray:
    """k_eff per site: number of stored neighbours with d <= alpha * d_(m_ref), clipped to [m_ref, K]."""
    nb_dist = np.asarray(nb_dist, dtype=np.float32)
    k = nb_dist.shape[1] if k_cap is None else min(k_cap, nb_dist.shape[1])
    if m_ref < 1 or m_ref > k:
        raise ValueError("m_ref must be within 1..K")
    radius = np.float32(alpha) * nb_dist[:, m_ref - 1, None]
    keff = np.sum(nb_dist[:, :k] <= radius, axis=1)
    return np.clip(keff, m_ref, k).astype(np.int16)


def fixed_k(nb_dist: np.ndarray, k: int) -> np.ndarray:
    """k_eff for a fixed-k rule (sensitivity comparisons)."""
    if k > nb_dist.shape[1]:
        raise ValueError("fixed k exceeds the stored neighbour count")
    return np.full(nb_dist.shape[0], k, dtype=np.int16)
