# MatPES element pools

The MatPES pools are not stored in this repository (tens of GB of derived shards); they are
rebuilt from the public releases with `scripts/build_matpes_pool.py`.  This document records
which elements were selected, the population rules, the sources and the reference T8 values.

## Elements (14)

| group | elements | first built | note |
|---|---|---|---|
| 3d, first set | **V, Co, Fe, Cr** | 2026-08-25 (PBE), 2026-08-27 (r2SCAN) | the four elements of the in-house materials (NVP, LiCoO2 / Co3O4, CrO2) |
| 3d, extension | **Mn, Ni, Cu** | 2026-09 | completes the magnetic late-3d row |
| early-3d / 4d | **Ti, Nb, Mo, Tc, Ru, Rh, Pd** | 2026-10-11 (both functionals) | validation of T20/T8 on elements not used in its development |

Order used in all tables: `V Co Fe Cr Mn Ni Cu Ti Nb Mo Tc Ru Rh Pd`.
The periodic-table figures cover periods 4–5, groups 3–12; the six remaining elements of that block
(**Sc, Zn, Y, Zr, Ag, Cd**) were shown from all atoms of the frames but never built as pools (their
moments are essentially zero).  MP-MM pools for the same 14 elements exist in the source project
(`mp-mm/`), scored with the same definition on 2026-10-28.

Machine-readable list: [`../data/matpes/elements.json`](../data/matpes/elements.json).

## Sources

| corpus | file | size (bytes) | sha256 |
|---|---|---:|---|
| MatPES-PBE 2025.2 (with charges) | `MatPES-PBE-2025.2-charges.json` | 2,358,718,427 | `fa3bb54eb05bd336251a130e02839cc7193fd45c92a4569755e7a194f01bf6fe` |
| MatPES-r2SCAN 2025.2 | `MatPES-R2SCAN-2025.2.json` | 1,946,549,525 | `d17878d1b5c33821b49a00412ccb7aabcdec3a67cea85552594954469680947f` |

Records: 433,189 (PBE), 386,544 (r2SCAN).  JSON layout per record: `structure` (pymatgen dict with
per-site `properties.magmom`), `forces`, `energy`, `stress`, `matpes_id`, `formula_pretty`,
`functional`, optional `bader_charges` / `bader_magmoms` (not used).

## Population rules (unchanged since 2026-08-25)

1. One pool per element: every frame with >= 1 atom of the element; all those atoms are sites.
2. Whole-frame drop if any atomic |magmom| > 8 muB (PBE: 0 frames, r2SCAN: 1 frame).
3. The separate isolated-atom reference file is **not** appended.
4. Duplicates retained; signed site moments retained (source final collinear moments).
5. Missing / non-finite moments and partial occupancies are fatal (none occurred).
6. Parent (material) id = MatPES id without its final underscore step (`..._999484_73` -> `..._999484`).
7. Charges stay in the source and never enter the metric.

Neighbour rule and statistics: [ALGORITHM.md](ALGORITHM.md) (per-structure SOAP basis, K = 100,
alpha 1.25, m 5, same-frame candidates, signed a_i, T8 with the 20 eV/A frame-force filter).

## Pool sizes

Sites per element (pool = after rules 1–5; "T8 population" = query sites after the 20 eV/A
frame-force filter).  Source: build manifests and `results/reference/matpes_T8_frameforce20_*.csv`.

| element | PBE pool | PBE T8 pop. | r2SCAN pool | r2SCAN T8 pop. |
|---|---:|---:|---:|---:|
| V  | 45,010 | 44,533 | 34,940 | 34,512 |
| Co | 67,334 | 67,169 | 49,639 | 49,524 |
| Fe | 101,889 | 101,637 | 79,380 | 79,205 |
| Cr | 39,367 | 39,095 | 31,862 | 31,679 |
| Mn | 72,440 | 72,179 | 55,886 | 55,672 |
| Ni | 78,596 | 78,266 | 64,544 | 64,360 |
| Cu | 51,390 | 51,261 | 42,118 | 42,034 |
| Ti | 38,938 | 38,810 | 32,153 | 32,065 |
| Nb | 28,587 | 28,466 | 24,067 | 23,955 |
| Mo | 29,028 | 28,896 | 25,176 | 25,085 |
| Tc | 12,596 | 12,592 | 11,896 | 11,892 |
| Ru | 23,982 | 23,977 | 21,139 | 21,134 |
| Rh | 24,543 | 24,532 | 21,397 | 21,388 |
| Pd | 26,146 | 26,110 | 22,409 | 22,375 |

## Reference T8 (muB, frame-force filter 20 eV/A, version 6cf4308f08ed0600)

| element | PBE | r2SCAN |   | element | PBE | r2SCAN |
|---|---:|---:|---|---|---:|---:|
| V  | 2.449 | 3.172 | | Ti | 0.850 | 1.316 |
| Co | 2.565 | 3.303 | | Nb | 0.812 | 1.367 |
| Fe | 3.939 | 4.948 | | Mo | 1.441 | 2.050 |
| Cr | 4.095 | 5.595 | | Tc | 1.659 | 2.581 |
| Mn | 5.604 | 6.792 | | Ru | 1.028 | 1.645 |
| Ni | 1.333 | 1.809 | | Rh | 0.456 | 0.789 |
| Cu | 0.559 | 0.799 | | Pd | 0.267 | 0.407 |

Full table with T20 and the mMACE force-MAE reduction y: `results/reference/matpes_T8_frameforce20_6cf4308f08ed0600.csv`.
MP-MM T8 of the same elements: `results/reference/mpmm_T8_14_elements_6cf4308f08ed0600_20261028.csv`.

## Rebuilding

```bash
# streams the JSON once (ijson), minutes per GB
python scripts/build_matpes_pool.py --source /path/to/MatPES-PBE-2025.2-charges.json \
    --out data/matpes/pools/pbe
python scripts/build_matpes_pool.py --source /path/to/MatPES-R2SCAN-2025.2.json \
    --out data/matpes/pools/r2scan
# then per element:
python -m magnetic_ambiguity run --input data/matpes/pools/pbe/Fe.extxyz --element Fe \
    --out results/matpes_pbe/Fe --max-abs-moment 8 --n-jobs 12
```
The builder writes `dft_magmom` and `REF_forces` per atom and `matpes_id`, `parent_id`,
`source_record` per frame; the frame-force filter reads `REF_forces`.  Differences to the source
project's builder: structures are converted directly (lattice matrix + Cartesian coordinates) instead
of through pymatgen, and one extxyz per element replaces the 2048-frame shards + parquet site tables;
sites, moments and forces are the same.
