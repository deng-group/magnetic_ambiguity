# Reference tables (regression targets)

Copied unchanged from the source project on 2026-10-30.

| file | content | source path (relative to the internal geometry-magnetic-mapping-benchmark project) |
|---|---|---|
| `inhouse_static_T3_T100_a6afa070e5c40baa.csv` | 5 in-house static populations: n_sites, A, median, C, T3..T100 (signed a_i), MACE / mMACE component MAE | `analysis/ambiguity_vs_error_decrease/results/Tn_methods_and_best_20261011/inhouse_T3_T100_descriptors.csv` |
| `inhouse_temperature_series_T8_a6afa070e5c40baa.csv` | Co3O4 and N2VP series: T8, T20, A, median, top-8 % count / threshold, k_eff | `../temperature/T8_ambiguity_20261020/tseries_T8_signed_a6afa070e5c40baa.csv` (first 14 columns) |
| `matpes_T8_frameforce20_6cf4308f08ed0600.csv` | 14 elements x PBE / r2SCAN: n_sites (fmax 20), T8, T20, y = mMACE force-MAE reduction (%) | `analysis/ambiguity_vs_error_decrease/results/p_convergence_T8_20261012/matpes_T8_T20_y_frameforce20.csv` |
| `mpmm_T8_14_elements_6cf4308f08ed0600_20261028.csv` | MP-MM 14 elements: pool / filtered sites, T8, T20, A, median, P, k_eff, unfiltered T8, MatPES T8 | `mp-mm/results/element_T8_ambiguity_6cf4308f08ed0600_20261028.csv` |

A rerun of this package on `data/inhouse/*.extxyz` with `--species inhouse --fmax none` reproduces the
in-house site tables bit for bit (verified 2026-10-30 on co3o4_50K, T8 0.00811124, and n2vp_50K,
T8 0.0255371; neighbour lists, distances, k_eff and a_i identical).  Spot-check with
`python -m magnetic_ambiguity run --input data/inhouse/Co3O4_900_MACE_atomic_forces.extxyz --element Co --out results/check/co3o4_900K --species inhouse --fmax none`
and compare `summary.csv:T8_muB` with 1.48618.
