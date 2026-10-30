# Decision history of the ambiguity definition

Dates are those of the internal source project (`geometry-magnetic-mapping-benchmark`).  Version hashes are
the source project's settings hashes (listed in `settings.LEGACY_VERSIONS`).

| date | change | version |
|---|---|---|
| 2026-08-11 | Training-free SOAP metric: SOAP(6 A, n 8, l 3), L2, global basis Co Cr Li Na O P V; a_i = median over the 20 nearest cross-frame neighbours of the signed |m_i − m_j|; A = mean a_i, C = A / median, P = frac(a_i > 0.5); 4,000-site subsample, duplicate frames removed | (fixed k 20) |
| 2026-08-20 | MP-MM corpus: per-compound SOAP basis + zero padding (bit-exact to global), |m| > 8 muB whole-frame drop, duplicates retained, true self-exclusion, exact blockwise kNN, settings-version caches | |
| 2026-08-21 / 09-04 | Adaptive band d_ij <= 1.25 d_i,5 with 5–20 neighbours becomes the only primary rule (selected by `14_neighbor_rule_selection`); fixed k 20 kept as provenance | a917fc084b83a8e6 (corpora), 4389b7138e1cc474 (in-house signed) |
| 2026-08-25 / 27 | MatPES-PBE and MatPES-r2SCAN pools V Co Fe Cr | |
| 2026-09-04 | Review: C fragile (median -> 0), AFM sign caveat; near-duplicate control; 6 A cutoff kept | |
| 2026-10-01 / 02 | Unsigned a_i = median ||m_i| − |m_j||; headline p = frac(a_i > 0.5 muB), sigma = median; (sigma, p) rule replaces (A, C) | 294407a65722bf16 (in-house) |
| 2026-10-09 | In-house: same-frame sites are candidates (only the query site excluded); K 20 -> 100 (band was capped for 30–86 % of sites) | a6afa070e5c40baa |
| 2026-10-11 | Corpora on the same rule; T20 (upper-20 % tail mean, signed) as the force-reduction predictor; Mn Ni Cu already present, Ti Nb Mo Tc Ru Rh Pd added | 6cf4308f08ed0600 |
| 2026-10-12 | Frame-force filter |F_DFT| > 20 eV/A on query sites; p-convergence scan fixes q = 8 % (T8) | |
| 2026-10-22 | alpha-rule convergence: alpha 1.25 on the rho(T8, MAE reduction) plateau for in-house and both MatPES functionals | |
| 2026-10-28 | **"Ambiguity" = T8** (signed-median a_i, same-frame K 100, alpha 1.25, fmax 20) for every table and figure; MP-MM scored for the 14 elements | |
| 2026-10-30 | This repository: algorithm consolidated as a package with these defaults; every earlier choice remains an option | |
