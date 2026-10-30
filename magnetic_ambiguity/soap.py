"""SOAP descriptors of the central-element sites, L2-normalised, stored block-sparse.

Convention (source project, decisions 2026-08-11 / 2026-08-20): dscribe SOAP, periodic, r_cut 6 A,
n_max 8, l_max 3, all other dscribe defaults; float64 -> sklearn L2 row normalisation -> float32.
Each structure may be described in its OWN species basis (species = sorted unique symbols of the
structure).  Because the descriptor blocks of species pairs absent from a structure are exactly
zero, this is identical to a global basis after normalisation (validated bit-exact in the source
project, gate V2), and the dot product between two sites is the sum over their SHARED pair blocks.
Storage is therefore block-sparse: {pair -> (site rows, float32 matrix)}.  A fixed global basis
(SoapSettings.species) is the special case where every site has every block.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Sequence

import numpy as np

from .settings import SoapSettings


@dataclass
class BlockSparseDescriptors:
    n_sites: int
    # pair key "A|B" with A <= B (sorted symbols) -> (row indices int64, X float32 (n_rows, width))
    blocks: dict[str, tuple[np.ndarray, np.ndarray]] = field(default_factory=dict)
    basis_of_site: list[str] = field(default_factory=list)   # "A|B|C" basis key per site
    # per basis key: dscribe feature layout {pair -> (start, stop)} and total width, so that a
    # dense matrix can be rebuilt in dscribe's native column order (bit-exact to a direct call)
    native_slices: dict[str, dict[str, tuple[int, int]]] = field(default_factory=dict)
    native_dim: dict[str, int] = field(default_factory=dict)

    def pairs_of_site(self, site: int) -> list[str]:
        species = self.basis_of_site[site].split("|")
        return [f"{a}|{b}" for i, a in enumerate(species) for b in species[i:]]

    def to_dense(self) -> np.ndarray:
        """Dense (n_sites, dim) matrix in dscribe's native feature order; valid only when every
        site shares one basis."""
        if len(set(self.basis_of_site)) != 1:
            raise ValueError("to_dense needs one common species basis for all sites")
        key = self.basis_of_site[0]
        pairs = self.pairs_of_site(0)
        for p in pairs:
            if not np.array_equal(self.blocks[p][0], np.arange(self.n_sites)):
                raise ValueError("block rows are not the full site range")
        if key in self.native_slices:
            dense = np.empty((self.n_sites, self.native_dim[key]), dtype=np.float32)
            for p, (lo, hi) in self.native_slices[key].items():
                dense[:, lo:hi] = self.blocks[p][1]
            return dense
        return np.concatenate([self.blocks[p][1] for p in pairs], axis=1)

    def dense_if_uniform(self):
        """Dense matrix when every site shares one basis (cached); None otherwise."""
        if not hasattr(self, "_dense_cache"):
            self._dense_cache = self.to_dense() if len(set(self.basis_of_site)) == 1 else None
        return self._dense_cache

    def dots(self, start: int, stop: int, out: np.ndarray) -> np.ndarray:
        """out[q - start, j] = <x_q, x_j> for query sites q in [start, stop) and all pool sites j.
        One basis for all sites: a single dense float32 matmul (the in-house reference convention).
        Several bases: accumulated block by block over the pair blocks shared by q and j into the
        float32 buffer `out` (zeroed here).  Blocks a site does not have are exactly zero in a
        global basis, so the sum over shared blocks IS the full dot product; only the float32
        summation order differs between the two paths (distances agree to ~1e-4 near 0)."""
        dense = self.dense_if_uniform()
        if dense is not None:
            np.matmul(dense[start:stop], dense.T, out=out)
            return out
        out[:] = 0.0
        for rows, x in self.blocks.values():
            in_query = (rows >= start) & (rows < stop)
            if not in_query.any():
                continue
            out[rows[in_query, None] - start, rows[None, :]] += x[in_query] @ x.T
        return out


def pair_key(a: str, b: str) -> str:
    return "|".join(sorted((a, b)))


def _soap_for(species: Sequence[str], settings: SoapSettings):
    from dscribe.descriptors import SOAP
    return SOAP(species=list(species), periodic=settings.periodic, r_cut=settings.r_cut,
                n_max=settings.n_max, l_max=settings.l_max)


def soap_descriptors(frames: Sequence, centers: Sequence[Sequence[int]], settings: SoapSettings,
                     n_jobs: int = 1, max_centers_per_call: int = 256,
                     log=None) -> BlockSparseDescriptors:
    """Normalised SOAP rows for `centers[f]` of every frame, in frame order, block-sparse.

    frames  : ASE Atoms (periodic)
    centers : per frame, the atom indices of the central-element sites (site order = row order)
    """
    from sklearn.preprocessing import normalize

    n_sites = int(sum(len(c) for c in centers))
    desc = BlockSparseDescriptors(n_sites=n_sites)
    soaps: dict[tuple[str, ...], object] = {}
    slices: dict[tuple[str, ...], dict[str, tuple[int, int]]] = {}
    block_rows: dict[str, list[np.ndarray]] = {}
    block_x: dict[str, list[np.ndarray]] = {}
    row = 0
    global_species = tuple(settings.species) if settings.species else None
    for f, (atoms, cen) in enumerate(zip(frames, centers)):
        cen = [int(c) for c in cen]
        if not cen:
            continue
        basis = global_species or tuple(sorted(set(atoms.get_chemical_symbols())))
        if global_species is not None:
            missing = set(atoms.get_chemical_symbols()) - set(global_species)
            if missing:
                raise ValueError(f"frame {f}: species {sorted(missing)} not in the SOAP species basis")
        if basis not in soaps:
            soap = _soap_for(basis, settings)
            soaps[basis] = soap
            slices[basis] = {}
            for i, a in enumerate(basis):
                for b in basis[i:]:
                    loc = soap.get_location((a, b))
                    slices[basis][pair_key(a, b)] = (int(loc.start), int(loc.stop))
            desc.native_slices["|".join(basis)] = dict(slices[basis])
            desc.native_dim["|".join(basis)] = int(soap.get_number_of_features())
        soap = soaps[basis]
        parts = []
        for start in range(0, len(cen), max_centers_per_call):
            parts.append(np.asarray(soap.create(atoms, centers=cen[start:start + max_centers_per_call],
                                                n_jobs=n_jobs), dtype=np.float64))
        x64 = np.concatenate(parts, axis=0) if len(parts) > 1 else parts[0]
        norms = np.linalg.norm(x64, axis=1)
        if (norms <= 1e-12).any():
            raise ValueError(f"frame {f}: {int((norms <= 1e-12).sum())} SOAP rows with zero norm")
        x32 = normalize(x64, norm="l2").astype(np.float32)
        rows = np.arange(row, row + len(cen), dtype=np.int64)
        key = "|".join(basis)
        desc.basis_of_site.extend([key] * len(cen))
        for pair, (lo, hi) in slices[basis].items():
            block_rows.setdefault(pair, []).append(rows)
            block_x.setdefault(pair, []).append(x32[:, lo:hi])
        row += len(cen)
        if log and (f + 1) % 50 == 0:
            log(f"    SOAP frame {f + 1}/{len(frames)}")
    if row != n_sites:
        raise ValueError("SOAP/site count mismatch")
    for pair in block_rows:
        desc.blocks[pair] = (np.concatenate(block_rows[pair]), np.concatenate(block_x[pair], axis=0))
    return desc


def descriptors_from_dense(x: np.ndarray) -> BlockSparseDescriptors:
    """Wrap an already L2-normalised dense matrix as a one-block descriptor set (tests, external descriptors)."""
    x = np.ascontiguousarray(x, dtype=np.float32)
    desc = BlockSparseDescriptors(n_sites=len(x), basis_of_site=["X"] * len(x))
    desc.blocks["X|X"] = (np.arange(len(x), dtype=np.int64), x)
    return desc
