"""Read an element pool from extxyz: central-element sites, moments, frame ids, frame max |F_DFT|.

Per-atom arrays looked up (first present wins unless fixed in PopulationSettings):
    moments : initial_magmoms, magmoms, dft_magmom, magmom, magmoms_dft        [muB, collinear]
    forces  : DFT_force, REF_forces, REF_force, forces, dft_forces             [eV/A]
Duplicate frames are retained (the metric counts them; decision 2026-08-20).  Frames with any
|m| > abs_moment_max_muB are dropped whole when that option is set (8 muB for MatPES / MP-MM).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence

import numpy as np
import pandas as pd

from .settings import FORCE_KEYS, MAGMOM_KEYS, PopulationSettings


@dataclass
class Pool:
    element: str
    frames: list                 # ASE Atoms, frames that survive the population rules, source order
    sites: pd.DataFrame          # one row per central-element site
    centers: list[list[int]]     # per frame, atom indices of the central-element sites
    magmom_key: str
    force_key: Optional[str]
    n_frames_read: int
    n_frames_dropped_abs_moment: int
    n_frames_without_element: int
    source: str


def _pick(arrays: dict, fixed: Optional[str], candidates: Sequence[str], what: str) -> Optional[str]:
    if fixed is not None:
        if fixed not in arrays:
            raise KeyError(f"{what} array '{fixed}' not in frame arrays {sorted(arrays)}")
        return fixed
    for key in candidates:
        if key in arrays:
            return key
    return None


def read_pool(path: str | Path, element: str, population: PopulationSettings = PopulationSettings(),
              frames: Optional[list] = None, require_forces: bool = False) -> Pool:
    from ase.io import read

    path = Path(path)
    if frames is None:
        frames = read(path, index=":", format="extxyz")
    kept, centers, rows = [], [], []
    n_drop_m, n_no_el = 0, 0
    magmom_key = force_key = None
    sid = 0
    for frame_id, atoms in enumerate(frames):
        if not atoms.pbc.all():
            raise ValueError(f"{path.name} frame {frame_id}: non-periodic cell")
        if atoms.info.get("config_type") == "IsolatedAtom":
            continue
        mk = _pick(atoms.arrays, population.magmom_key, MAGMOM_KEYS, "magmom")
        if mk is None:
            raise KeyError(f"{path.name} frame {frame_id}: no magnetic-moment array among {MAGMOM_KEYS}")
        fk = _pick(atoms.arrays, population.force_key, FORCE_KEYS, "force")
        if fk is None and require_forces:
            raise KeyError(f"{path.name} frame {frame_id}: no DFT force array among {FORCE_KEYS}")
        magmom_key, force_key = mk, fk
        moments = np.asarray(atoms.arrays[mk], float).reshape(len(atoms), -1)
        if moments.shape[1] != 1:
            raise ValueError(f"{path.name} frame {frame_id}: moments must be collinear scalars, got shape {moments.shape}")
        moments = moments[:, 0]
        if not np.isfinite(moments).all():
            raise ValueError(f"{path.name} frame {frame_id}: non-finite moment")
        if population.abs_moment_max_muB is not None and np.any(np.abs(moments) > population.abs_moment_max_muB):
            n_drop_m += 1
            continue
        symbols = np.asarray(atoms.get_chemical_symbols())
        selected = np.flatnonzero(symbols == element)
        if len(selected) == 0:
            n_no_el += 1
            continue
        fmax = float(np.linalg.norm(np.asarray(atoms.arrays[fk], float), axis=1).max()) if fk else np.nan
        kept_id = len(kept)
        kept.append(atoms)
        centers.append(selected.tolist())
        for atom_index in selected:
            rows.append({"site_id": sid, "frame_id": kept_id, "source_frame_index": frame_id,
                         "atom_index": int(atom_index), "element": element,
                         "magmom_muB": float(moments[atom_index]), "frame_fmax_eV_A": fmax,
                         "n_atoms_frame": len(atoms)})
            sid += 1
    if not rows:
        raise ValueError(f"{path.name}: no {element} sites after the population rules")
    sites = pd.DataFrame(rows)
    return Pool(element=element, frames=kept, sites=sites, centers=centers, magmom_key=magmom_key,
                force_key=force_key, n_frames_read=len(frames), n_frames_dropped_abs_moment=n_drop_m,
                n_frames_without_element=n_no_el, source=str(path.resolve()))


def central_site_force_errors(pool: Pool, predicted_key: str = "MACE_force") -> Optional[pd.DataFrame]:
    """Optional: per-site force error of a model (predicted - DFT) when both arrays exist (in-house files)."""
    if pool.force_key is None or any(predicted_key not in a.arrays for a in pool.frames):
        return None
    err = []
    for fid, atoms in enumerate(pool.frames):
        d = np.asarray(atoms.arrays[predicted_key], float) - np.asarray(atoms.arrays[pool.force_key], float)
        for atom_index in pool.centers[fid]:
            e = d[atom_index]
            err.append({"error_fx_eV_A": e[0], "error_fy_eV_A": e[1], "error_fz_eV_A": e[2],
                        "force_vector_error_eV_A": float(np.linalg.norm(e))})
    return pd.DataFrame(err)
