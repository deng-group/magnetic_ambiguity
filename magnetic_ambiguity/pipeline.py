"""One element pool in, one site table + one summary row out."""

from __future__ import annotations

import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from . import __version__
from .dataset_ambiguity import dataset_statistics, frame_force_mask
from .io import Pool, central_site_force_errors, read_pool
from .neighbors import adaptive_k, exact_knn
from .settings import Settings
from .site_ambiguity import site_table
from .soap import soap_descriptors


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(16 << 20), b""):
            h.update(block)
    return h.hexdigest()


def _atomic_write_text(path: Path, text: str) -> None:
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(text)
    os.replace(tmp, path)


def summarize_sites(sites: pd.DataFrame, settings: Settings, column: str = "ambiguity_muB") -> dict:
    """Dataset statistics with the frame-force filter applied to QUERY sites only."""
    ds = settings.dataset
    keep = frame_force_mask(sites["frame_fmax_eV_A"].to_numpy(), ds.fmax_filter_eV_A)
    a = sites.loc[keep, column].to_numpy(float)
    row = dataset_statistics(a, sites.loc[keep, "adaptive_k"].to_numpy(), tail_fraction=ds.tail_fraction,
                             extra_tail_fractions=ds.extra_tail_fractions, tau=ds.tau_muB,
                             k_cap=settings.neighbors.k_exact)
    n_frames = sites["frame_id"].nunique()
    n_frames_kept = sites.loc[keep, "frame_id"].nunique()
    row.update({
        "n_sites_pool": int(len(sites)),
        "n_frames_pool": int(n_frames),
        "n_frames_dropped_fmax": int(n_frames - n_frames_kept),
        "fmax_filter_eV_A": ds.fmax_filter_eV_A,
        "n_frames_without_force_record": int(sites.loc[sites["frame_fmax_eV_A"].isna(), "frame_id"].nunique()),
        "moments": settings.ambiguity.moments,
        "mean_abs_moment_muB": float(np.abs(sites.loc[keep, "magmom_muB"]).mean()),
        "same_frame_twins_mean": float(sites.loc[keep, "n_same_frame_twins"].mean()),
        "fraction_sites_with_same_frame_twin": float((sites.loc[keep, "n_same_frame_twins"] > 0).mean()),
    })
    if ds.fmax_filter_eV_A is not None:
        row[f"T{round(ds.tail_fraction * 100):d}_unfiltered_muB"] = dataset_statistics(
            sites[column].to_numpy(float), tail_fraction=ds.tail_fraction, extra_tail_fractions=())[
            f"T{round(ds.tail_fraction * 100):d}_muB"]
    if "error_fx_eV_A" in sites:
        res = sites.loc[keep, ["error_fx_eV_A", "error_fy_eV_A", "error_fz_eV_A"]].to_numpy(float)
        row["force_component_MAE_eV_A"] = float(np.mean(np.abs(res)))
    return row


def run_pool(pool: Pool, settings: Settings, n_jobs: int = 1, query_chunk: int = 4096,
             log=print) -> tuple[pd.DataFrame, dict, dict]:
    """SOAP -> exact kNN -> adaptive twins -> a_i -> dataset statistics.  Returns (sites, summary, neighbours)."""
    ns = settings.neighbors
    log(f"{pool.element}: {len(pool.frames)} frames, {len(pool.sites)} sites; settings version {settings.version()}")
    desc = soap_descriptors(pool.frames, pool.centers, settings.soap, n_jobs=n_jobs, log=log)
    log(f"  SOAP done: {len(desc.blocks)} pair blocks, {len(set(desc.basis_of_site))} species bases")
    nb_idx, nb_dist = exact_knn(desc, pool.sites["frame_id"].to_numpy(), ns.k_exact,
                                exclude_same_frame=ns.exclude_same_frame, query_chunk=query_chunk, log=log)
    k_eff = adaptive_k(nb_dist, ns.alpha, ns.m_ref, ns.k_exact)
    sites = site_table(pool.sites, nb_idx, nb_dist, k_eff, headline=settings.ambiguity.moments)
    errors = central_site_force_errors(pool)
    if errors is not None:
        sites = pd.concat([sites, errors], axis=1)
    summary = {"element": pool.element, "source": pool.source, "settings_version": settings.version(),
               "n_frames_read": pool.n_frames_read, "n_frames_dropped_abs_moment": pool.n_frames_dropped_abs_moment,
               "n_frames_without_element": pool.n_frames_without_element,
               **summarize_sites(sites, settings)}
    log(f"  T{round(settings.dataset.tail_fraction * 100)} = {summary[f'T{round(settings.dataset.tail_fraction * 100)}_muB']:.4f} muB"
        f"  (A {summary['A_muB']:.4f}, median {summary['median_muB']:.4f}, k_eff median {summary['k_eff_median']:.0f})")
    return sites, summary, {"neighbor_site_ids": nb_idx, "neighbor_distances": nb_dist}


def write_outputs(out_dir: Path, pool: Pool, settings: Settings, sites: pd.DataFrame, summary: dict,
                  neighbors: dict, save_neighbors: bool = True) -> None:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    sites.to_parquet(out_dir / "sites.parquet", index=False, compression="zstd")
    if save_neighbors:
        np.savez_compressed(out_dir / "neighbors.npz", **neighbors)
    pd.DataFrame([summary]).to_csv(out_dir / "summary.csv", index=False, float_format="%.15g")
    meta = {"status": "complete", "created_utc": datetime.now(timezone.utc).isoformat(),
            "package_version": __version__, "settings_version": settings.version(),
            "soap_version": settings.soap_version(), "settings": settings.to_dict(),
            "source": {"path": pool.source, "sha256": sha256_file(Path(pool.source)),
                       "magmom_key": pool.magmom_key, "force_key": pool.force_key},
            "outputs": {"sites": str((out_dir / "sites.parquet").resolve()),
                        "neighbors": str((out_dir / "neighbors.npz").resolve()) if save_neighbors else None}}
    _atomic_write_text(out_dir / "metadata.json", json.dumps(meta, indent=2, default=str) + "\n")


def rescore(out_dir: Path, settings: Settings, alpha: Optional[float] = None, m_ref: Optional[int] = None,
            fixed_k: Optional[int] = None) -> tuple[pd.DataFrame, dict]:
    """Recompute a_i and statistics from stored neighbour lists with another band rule (no new kNN)."""
    from .neighbors import fixed_k as fixed_rule
    out_dir = Path(out_dir)
    sites = pd.read_parquet(out_dir / "sites.parquet")
    nb = np.load(out_dir / "neighbors.npz")
    nb_idx, nb_dist = nb["neighbor_site_ids"], nb["neighbor_distances"]
    base_cols = ["site_id", "frame_id", "source_frame_index", "atom_index", "element", "magmom_muB",
                 "frame_fmax_eV_A", "n_atoms_frame"]
    base = sites[[c for c in base_cols if c in sites]]
    if fixed_k is not None:
        k_eff = fixed_rule(nb_dist, fixed_k)
    else:
        k_eff = adaptive_k(nb_dist, alpha if alpha is not None else settings.neighbors.alpha,
                           m_ref if m_ref is not None else settings.neighbors.m_ref, nb_dist.shape[1])
    new = site_table(base, nb_idx, nb_dist, k_eff, headline=settings.ambiguity.moments)
    for c in ("error_fx_eV_A", "error_fy_eV_A", "error_fz_eV_A"):
        if c in sites:
            new[c] = sites[c].to_numpy()
    summary = summarize_sites(new, settings)
    summary.update({"alpha": None if fixed_k is not None else (alpha if alpha is not None else settings.neighbors.alpha),
                    "m_ref": None if fixed_k is not None else (m_ref if m_ref is not None else settings.neighbors.m_ref),
                    "fixed_k": fixed_k})
    return new, summary


def run_file(path: Path, element: str, settings: Settings, out_dir: Path, n_jobs: int = 1,
             query_chunk: int = 4096, save_neighbors: bool = True, log=print) -> dict:
    pool = read_pool(path, element, settings.population)
    sites, summary, neighbors = run_pool(pool, settings, n_jobs=n_jobs, query_chunk=query_chunk, log=log)
    write_outputs(out_dir, pool, settings, sites, summary, neighbors, save_neighbors=save_neighbors)
    return summary
