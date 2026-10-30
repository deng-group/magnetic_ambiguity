# magnetic_ambiguity

How much of a DFT dataset's magnetism is **not determined by its geometry**?  For every
magnetic-ion site we look for its geometric twins (nearest neighbours in SOAP space) and
ask how far their DFT moments disagree; the dataset number is the mean ambiguity of the
most ambiguous 8 % of sites.

```
a_i  = median_{j in twins(i)} | m_i - m_j |              site ambiguity   [muB]  (signed moments)
twins(i): exact K = 100 SOAP neighbours of i in the element pool, kept if d_ij <= 1.25 d_i,5
          (5 <= k_eff <= 100); same-frame sites are candidates, only i itself is excluded
T8   = mean of the upper 8 % mass of {a_i}                 dataset ambiguity [muB]
       (query sites in frames with any |F_DFT| > 20 eV/A excluded; neighbours from the full pool)
```

This is the definition in force since 2026-10-28 ("ambiguity" = T8).  Every parameter is a
command-line option; the defaults are the current ones.  Full derivation, alternatives and the
decision history: [docs/ALGORITHM.md](docs/ALGORITHM.md), [docs/HISTORY.md](docs/HISTORY.md).

## Layout

```
magnetic_ambiguity/      the package (SOAP -> exact kNN -> adaptive twins -> a_i -> T_p)
  settings.py            all parameters + 16-hex settings version hash (defaults = 2026-10-28 definition)
  soap.py                dscribe SOAP, L2-normalised, per-structure species basis, block-sparse store
  neighbors.py           exact blockwise kNN (no approximate search), adaptive band rule
  site_ambiguity.py      a_i signed / unsigned, k_eff, same-frame twin counts
  dataset_ambiguity.py   tail(), T_p, A, median, C, P(>tau), frame-force filter, legacy (sigma, p) rule
  io.py / pipeline.py    extxyz pool reader, run / rescore / summarize, provenance metadata
  inhouse.py             registry of the 29 in-house populations
  cli.py                 python -m magnetic_ambiguity {run,rescore,summarize,settings}
scripts/run_inhouse_all.py      all in-house datasets with the in-house settings
scripts/build_matpes_pool.py    MatPES JSON -> per-element extxyz pools (population rules applied)
data/inhouse/                   the 30 in-house extxyz files (moments + DFT + MACE forces), manifest, checksums
data/matpes/                    element list and population rules of the MatPES pools (data rebuilt from the public JSON)
docs/                           ALGORITHM, DATASETS, MATPES_ELEMENTS, HISTORY
results/reference/              reference T8 tables of the source project (regression targets)
tests/                          numpy-only unit tests (tail, band rule, kNN vs brute force, block-sparse dot)
```

## Install

```bash
python -m venv .venv && .venv/bin/python -m pip install -r requirements.txt   # pinned metric stack
.venv/bin/python -m pytest tests          # or: .venv/bin/python tests/test_core.py
```
Python 3.12 was used throughout.

## Run

```bash
# one pool, current definition (T8, K 100, alpha 1.25, m 5, same-frame candidates, signed, fmax 20)
python -m magnetic_ambiguity run --input data/inhouse/LiCoO2_all_U_MACE_atomic_forces.extxyz \
    --element Co --out results/licoo2_u --species inhouse --fmax none --n-jobs 4

# the in-house convention explicitly, and the options behind the defaults
python -m magnetic_ambiguity run --input POOL.extxyz --element Fe --out OUT \
    --k 100 --alpha 1.25 --m-ref 5 --moments signed --tail 0.08 --fmax 20 [--max-abs-moment 8]

# superseded rules, for comparison
python -m magnetic_ambiguity run ... --k 20 --exclude-same-frame            # cross-frame K = 20 (2026-08/09)
python -m magnetic_ambiguity run ... --moments unsigned                     # unsigned a_i (2026-10-01 branch)
python -m magnetic_ambiguity rescore --dir OUT --alpha 1.5                  # other band from stored neighbours
python -m magnetic_ambiguity rescore --dir OUT --fixed-k 20
python -m magnetic_ambiguity summarize --sites OUT/sites.parquet --tail 0.20 --fmax none

# all 29 in-house populations (resumable)
python scripts/run_inhouse_all.py --out results/inhouse --n-jobs 12

# a MatPES element pool
python scripts/build_matpes_pool.py --source /path/to/MatPES-PBE-2025.2-charges.json \
    --out data/matpes/pools/pbe --elements V Co Fe Cr Mn Ni Cu Ti Nb Mo Tc Ru Rh Pd
python -m magnetic_ambiguity run --input data/matpes/pools/pbe/Fe.extxyz --element Fe --out results/matpes_pbe/Fe
```

Outputs per run: `sites.parquet` (one row per site: moment, `adaptive_k`, `ambiguity_signed_muB`,
`ambiguity_unsigned_muB`, `ambiguity_muB` = headline, neighbour-median moment, same-frame twin
counts, frame max |F|, central-site force error when a model force array exists), `neighbors.npz`
(exact K neighbour ids and distances, for `rescore`), `summary.csv` (T8, T20, A, median, C,
P(a_i > 0.5), k_eff statistics, filter counts), `metadata.json` (settings, version hashes, source sha256).

A pool of 10^5 sites needs a few GB of memory (`--query-chunk` x pool size x 4 bytes for the distance
block) and minutes to an hour on 12 cores; the in-house datasets run in seconds to minutes each.

## Input format

extxyz, periodic frames, per-atom arrays: a collinear moment array (`initial_magmoms`, `magmoms`,
`dft_magmom`, `magmom`; the in-house files use `initial_magmoms`, the MatPES pools `dft_magmom`),
a DFT force array for the frame-force filter (`DFT_force`, `REF_forces`, `forces`), optionally a
model force array `MACE_force` for the per-site error columns.  Keys can be fixed with
`--magmom-key` / `--force-key`.  Frames with `config_type=IsolatedAtom` are skipped.

## Reference values (regression targets, `results/reference/`)

Validation 2026-10-30: `run --species inhouse --fmax none` on Co3O4 50 K and N2VP 50 K reproduces the
stored reference site tables (descriptors, 100-neighbour lists and distances, k_eff, signed and
unsigned a_i, T8) bit for bit.

| population | n sites | T8 (muB) | source version |
|---|---:|---:|---|
| in-house CrO2 +U (Cr) | 5,760 | 0.231 | a6afa070e5c40baa, no force filter |
| in-house Li0.5CoO2 +U (Co) | 4,528 | 1.713 | a6afa070e5c40baa, no force filter |
| in-house LiCoO2 +U (Co) | 3,648 | 1.420 | a6afa070e5c40baa, no force filter |
| in-house NVP no U (V) | 20,496 | 0.601 | a6afa070e5c40baa, no force filter |
| in-house NVP +U (V) | 20,680 | 0.856 | a6afa070e5c40baa, no force filter |
| in-house Co3O4 900 K (Co) | 1,920 | 1.486 | a6afa070e5c40baa, no force filter |
| MatPES-PBE Fe | 101,637 | 3.939 | 6cf4308f08ed0600, fmax 20 |
| MatPES-r2SCAN Cr | 31,679 | 5.595 | 6cf4308f08ed0600, fmax 20 |
| MP-MM Mn | 126,398 | 8.396 | 6cf4308f08ed0600, fmax 20 |

The in-house values were computed with a global 7-species SOAP basis, no |m| filter and no
frame-force filter (`--species inhouse --fmax none`); the corpora with the per-structure basis, the
whole-frame |m| > 8 muB drop and the 20 eV/A frame-force filter (`--max-abs-moment 8`, default `--fmax 20`).
Both bases give identical descriptors after L2 normalisation.

## Provenance

Consolidated on 2026-10-30 from the internal `geometry-magnetic-mapping-benchmark` project
(`analysis/dataset_ambiguity`, `matpes`, `mp-mm`, `analysis/ambiguity_vs_error_decrease`).  The
algorithms are reproduced function by function; file layouts and names are new.  The source
project's version hashes (`a6afa070e5c40baa`, `6cf4308f08ed0600`, ...) are listed in
`settings.LEGACY_VERSIONS` for cross-reference; this package computes its own hash from a
cleaner settings layout.
