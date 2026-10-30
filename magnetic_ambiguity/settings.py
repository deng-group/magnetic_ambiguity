"""All output-affecting settings in one place, with a stable 16-hex version hash.

Defaults are the parameters in use since 2026-10-28.  Infrastructure knobs (chunk sizes,
thread counts) are deliberately NOT part of the settings so that they never change the
version hash; they were shown to be output-identical in the source project (V4/V5 gates).

Legacy version hashes of the source project, for cross-reference only (their payloads had a
different key layout, so this module does not reproduce them):
    a6afa070e5c40baa  in-house, same-frame K=100, alpha 1.25 (2026-10-09)
    6cf4308f08ed0600  MatPES / MP-MM corpora, same rule (2026-10-11)
    294407a65722bf16  in-house, cross-frame K=20 (superseded 2026-10-09)
    a917fc084b83a8e6  corpora, cross-frame K=20 (superseded 2026-10-11)
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass, field, replace
from typing import Optional

LEGACY_VERSIONS = {
    "inhouse_sameframe_k100": "a6afa070e5c40baa",
    "corpus_sameframe_k100": "6cf4308f08ed0600",
    "inhouse_crossframe_k20": "294407a65722bf16",
    "corpus_crossframe_k20": "a917fc084b83a8e6",
}

# Species basis used by every in-house run; a per-structure basis (species=None) is
# numerically identical after L2 normalisation (zero blocks for absent pairs).
INHOUSE_SOAP_SPECIES = ("Co", "Cr", "Li", "Na", "O", "P", "V")


@dataclass(frozen=True)
class SoapSettings:
    r_cut: float = 6.0
    n_max: int = 8
    l_max: int = 3
    periodic: bool = True
    # None -> each structure is described in its own sorted-unique species basis
    # (per-compound basis, block-sparse storage).  A tuple fixes one global basis.
    species: Optional[tuple[str, ...]] = None
    normalization: str = "sklearn l2 on float64, stored float32"


@dataclass(frozen=True)
class NeighborSettings:
    k_exact: int = 100          # exact neighbours stored (candidate cap of the adaptive band)
    alpha: float = 1.25         # twins: d_ij <= alpha * d_i,m_ref
    m_ref: int = 5              # reference neighbour of the adaptive radius; also the minimum k_eff
    exclude_same_frame: bool = False   # False: only the query site is excluded (rule since 2026-10-09/11)
    distance: str = "sqrt(clip(2 - 2*dot, 0)) on L2-normalised SOAP, float32"
    tie_order: str = "lexicographic (distance, site_id)"


@dataclass(frozen=True)
class AmbiguitySettings:
    moments: str = "signed"     # "signed": |m_i - m_j|;  "unsigned": ||m_i| - |m_j||
    statistic: str = "median over adaptive twins"


@dataclass(frozen=True)
class DatasetSettings:
    tail_fraction: float = 0.08         # T_p: mean of the upper-p mass of a_i (T8)
    extra_tail_fractions: tuple[float, ...] = (0.20,)
    fmax_filter_eV_A: Optional[float] = 20.0   # drop QUERY sites of frames with any |F_DFT| above this; None disables
    tau_muB: float = 0.5                # threshold of P(a_i > tau) and of the legacy (sigma, p) rule


@dataclass(frozen=True)
class PopulationSettings:
    abs_moment_max_muB: Optional[float] = None  # whole-frame drop if any |m| above (8.0 for MatPES/MP-MM); None keeps all
    duplicates: str = "retained"
    magmom_key: Optional[str] = None            # extxyz per-atom array; None -> first of MAGMOM_KEYS present
    force_key: Optional[str] = None             # extxyz per-atom DFT force array; None -> first of FORCE_KEYS present


@dataclass(frozen=True)
class Settings:
    soap: SoapSettings = field(default_factory=SoapSettings)
    neighbors: NeighborSettings = field(default_factory=NeighborSettings)
    ambiguity: AmbiguitySettings = field(default_factory=AmbiguitySettings)
    dataset: DatasetSettings = field(default_factory=DatasetSettings)
    population: PopulationSettings = field(default_factory=PopulationSettings)
    schema_version: int = 1

    def to_dict(self) -> dict:
        return asdict(self)

    def version(self) -> str:
        """16-hex sha256 of the output-affecting settings (SOAP, neighbours, ambiguity, population)."""
        payload = {k: v for k, v in self.to_dict().items() if k != "dataset"}
        return hashlib.sha256(json.dumps(payload, sort_keys=True, default=str).encode()).hexdigest()[:16]

    def soap_version(self) -> str:
        return hashlib.sha256(json.dumps(asdict(self.soap), sort_keys=True, default=str).encode()).hexdigest()[:16]

    def with_(self, **sections) -> "Settings":
        """Return a copy with whole sections replaced, e.g. s.with_(neighbors=NeighborSettings(k_exact=20))."""
        return replace(self, **sections)


MAGMOM_KEYS = ("initial_magmoms", "magmoms", "dft_magmom", "magmom", "magmoms_dft")
FORCE_KEYS = ("DFT_force", "REF_forces", "REF_force", "forces", "dft_forces")


def inhouse_settings(**overrides) -> Settings:
    """Settings of the in-house runs: global 7-species SOAP basis, no |m| filter."""
    s = Settings(soap=SoapSettings(species=INHOUSE_SOAP_SPECIES))
    return s.with_(**overrides) if overrides else s


def corpus_settings(**overrides) -> Settings:
    """Settings of the MatPES / MP-MM runs: per-compound SOAP basis, whole-frame |m| > 8 muB drop."""
    s = Settings(population=PopulationSettings(abs_moment_max_muB=8.0))
    return s.with_(**overrides) if overrides else s
