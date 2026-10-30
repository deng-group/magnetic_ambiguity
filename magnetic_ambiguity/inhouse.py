"""Registry of the in-house DFT datasets shipped in data/inhouse/ (29 ambiguity populations).

Every file carries per-atom `initial_magmoms` (final collinear DFT moments, muB), `DFT_force` and a
model force `MACE_force` (eV/A).  Two files duplicate the DFT content of another entry and are
therefore excluded from the ambiguity populations (kept for force comparisons only):
    NVP_all_U31_MagMACE  same DFT as NVP_all_U31_MACE (only MACE_force differs)
    N2VP_800             80-frame structural subset of NVP_x2_U31_x2
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

DATA_DIR = Path(__file__).resolve().parents[1] / "data" / "inhouse"


@dataclass(frozen=True)
class InhouseDataset:
    key: str
    filename: str
    element: str
    material: str
    analysis_group: str          # "static" or "temperature"
    hubbard_u: str               # "yes", "no", "yes (U=3.1 eV)"
    temperature_K: int | None = None
    label: str = ""

    @property
    def path(self) -> Path:
        return DATA_DIR / self.filename


STATIC = (
    InhouseDataset("cro2_u", "CrO2_all_U_MACE_atomic_forces.extxyz", "Cr", "CrO2", "static", "yes", label="CrO2 +U"),
    InhouseDataset("cro2_nou", "CrO2_all_noU_MACE_atomic_forces.extxyz", "Cr", "CrO2", "static", "no", label="CrO2 no U"),
    InhouseDataset("licoo2_u", "LiCoO2_all_U_MACE_atomic_forces.extxyz", "Co", "LiCoO2", "static", "yes", label="LiCoO2 +U"),
    InhouseDataset("licoo2_nou", "LiCoO2_all_noU_MACE_atomic_forces.extxyz", "Co", "LiCoO2", "static", "no", label="LiCoO2 no U"),
    InhouseDataset("li05coo2_u", "Li0.5CoO2_U_MACE_atomic_forces.extxyz", "Co", "Li0.5CoO2", "static", "yes", label="Li0.5CoO2 +U"),
    InhouseDataset("nvp_all_u31", "NVP_all_U31_MACE_atomic_forces.extxyz", "V", "NVP", "static", "yes (U=3.1 eV)", label="NVP +U"),
    InhouseDataset("nvp_all_u0", "NVP_all_U0_MACE_atomic_forces_matched.extxyz", "V", "NVP", "static", "no", label="NVP no U"),
    InhouseDataset("nvp_x2_u31", "NVP_x2_U31_x2_MACE_atomic_forces.extxyz", "V", "NVP-x2", "static", "yes (U=3.1 eV)", 800, "NVP x2 +U"),
)

TEMPERATURE = tuple(
    [InhouseDataset(f"co3o4_{t}K", f"Co3O4_{t}_MACE_atomic_forces.extxyz", "Co", "Co3O4", "temperature", "yes", t, f"Co3O4 {t} K")
     for t in (50, 100, 150, 200, 250, 400, 500, 600, 700, 800, 900)]
    + [InhouseDataset(f"n2vp_{t}K", f"N2VP_{t}_MACE_atomic_forces.extxyz", "V", "N2VP", "temperature", "yes", t, f"N2VP {t} K")
       for t in (50, 100, 200, 300, 400, 500, 600, 700, 800, 900)]
)

DATASETS = STATIC + TEMPERATURE
BY_KEY = {d.key: d for d in DATASETS}

# force-comparison only (same DFT structures as an ambiguity population above)
FORCE_ONLY = (
    InhouseDataset("nvp_all_u31_magmace", "NVP_all_U31_MagMACE_atomic_forces.extxyz", "V", "NVP", "force_only", "yes (U=3.1 eV)", label="NVP +U (MagMACE forces)"),
    InhouseDataset("n2vp_800_subset", "N2VP_800_MACE_atomic_forces.extxyz", "V", "N2VP", "force_only", "yes", 800, "N2VP 800 K subset"),
)
