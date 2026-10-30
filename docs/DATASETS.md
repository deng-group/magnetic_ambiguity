# In-house datasets (`data/inhouse/`)

30 extxyz files (179 MB), the force/moment files of the source project's `data/atomic_force_errors/`.
Every frame is periodic and carries per-atom `initial_magmoms` (final collinear DFT moments, muB,
3 decimals), `DFT_force` and `MACE_force` (eV/A; `REF_forces` duplicates `DFT_force`).  The signed
model error is `MACE_force - DFT_force`.  Integrity: `cd data/inhouse && sha256sum -c checksums.sha256`
(hashes are those of the source project's `checksums.sha256`).

## Ambiguity populations (29)

Registry: `magnetic_ambiguity/inhouse.py`.  Each file is one neighbour pool.

| key | file | element | frames | sites | Hubbard U | note |
|---|---|---|---:|---:|---|---|
| cro2_u | CrO2_all_U_MACE_atomic_forces.extxyz | Cr | 160 | 5,760 | yes | |
| cro2_nou | CrO2_all_noU_MACE_atomic_forces.extxyz | Cr | 160 | 5,760 | no | |
| licoo2_u | LiCoO2_all_U_MACE_atomic_forces.extxyz | Co | 152 | 3,648 | yes | |
| licoo2_nou | LiCoO2_all_noU_MACE_atomic_forces.extxyz | Co | 160 | 3,840 | no | |
| li05coo2_u | Li0.5CoO2_U_MACE_atomic_forces.extxyz | Co | 283 | 4,528 | yes | Li0.5CoO2 subset of the mixed-composition LiCoO2 source; 74 duplicate frames retained |
| nvp_all_u31 | NVP_all_U31_MACE_atomic_forces.extxyz | V | 2,600 | 20,680 | U = 3.1 eV | all compositions, ordinary MACE forces |
| nvp_all_u0 | NVP_all_U0_MACE_atomic_forces_matched.extxyz | V | 2,562 | 20,496 | no | 38 unconverged frames removed (`NVP_all_noU_excluded_unconverged_structures.csv`) |
| nvp_x2_u31 | NVP_x2_U31_x2_MACE_atomic_forces.extxyz | V | 270 | 2,160 | U = 3.1 eV | 270-frame x2 subset, x2-specific model |
| co3o4_{50..900}K | Co3O4_{T}_MACE_atomic_forces.extxyz | Co | 51–186 | 1,224–4,464 | yes | T = 50 100 150 200 250 400 500 600 700 800 900 K |
| n2vp_{50..900}K | N2VP_{T}_MACE_atomic_forces.extxyz | V | 80 | 640 | yes | T = 50 100 200 300 400 500 600 700 800 900 K |

Frame counts from `manifest.csv`; site counts from the reference tables (sites = frames x magnetic
atoms per cell).  Two further files share DFT content with the above and are **force-comparison
only** (not ambiguity populations): `NVP_all_U31_MagMACE_atomic_forces.extxyz` (same DFT as
`NVP_all_U31_MACE`, MagMACE forces) and `N2VP_800_MACE_atomic_forces.extxyz` (80-frame structural
subset of `NVP_x2_U31_x2`, predictions from the temperature-series model).

## Settings used for the reference values

Global SOAP basis `Co Cr Li Na O P V` (`--species inhouse`), K = 100, alpha 1.25, m 5, same-frame
candidates, signed moments, T8 over **all** sites (the in-house reference tables were produced before the
frame-force filter existed and never applied it: pass `--fmax none` to reproduce them), no |m| filter.
Source version `a6afa070e5c40baa`.

Reference T8 (muB), static populations (`results/reference/inhouse_static_T3_T100_a6afa070e5c40baa.csv`,
column `T8`): CrO2 +U (Cr) 0.231, Li0.5CoO2 +U (Co) 1.713, LiCoO2 +U (Co) 1.420, NVP no U (V) 0.601, NVP +U (V) 0.856.  Temperature series
(`results/reference/inhouse_temperature_series_T8_a6afa070e5c40baa.csv`): Co3O4 0.008 (50 K) ->
1.486 (900 K), monotonic; N2VP 0.026 (50 K) -> 0.835 (800 K), 0.634 at 900 K.

## Not copied from the source project

`data/n3vp/`, `data/cross_material/`, `data/temperature_series/` of the source repository hold
atom-wise **MACE descriptors** (and the N3VP subset structures) used by the earlier MACE-descriptor
metrics; they are not inputs of the SOAP ambiguity and were left in the internal project.  The NVP mMACE model (84 MB) was also left there.
