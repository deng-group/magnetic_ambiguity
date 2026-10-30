# Geometry-magnetic ambiguity: the algorithm

Status: definition in force since 2026-10-28.  Defaults of the package = this document.
Alternatives that were used before are kept as options and are listed in §7 and in
[HISTORY.md](HISTORY.md).

## 0. Question

A geometry-only interatomic potential (MACE, ...) predicts forces from positions and species.
If two sites with indistinguishable local geometry carry different DFT magnetic moments, the
force label of one of them cannot be learned from geometry alone.  "Ambiguity" measures how
much moment is undetermined by geometry in a dataset, with DFT data only (no model, no fit).

## 1. Population: element pools

* One **pool per central element** (V, Co, Fe, ...).  A site is every atom of that element in
  every frame of the dataset.  Neighbours are drawn from the same pool only (same element).
* **Duplicates are retained** (exact duplicate frames count; decision 2026-08-20).
* **Corpora** (MatPES, MP-MM): a whole frame is dropped if any atom has |m| > 8 muB
  (`--max-abs-moment 8`; paper Methods of MP-MM, 0–1 frames in MatPES); the separate
  isolated-atom reference file is not appended; missing / non-finite moments and partial
  occupancies abort the build.  **In-house** datasets: no |m| filter.
* Frames with `config_type = IsolatedAtom` are skipped by the reader.
* Moments are the **final collinear DFT moments** in muB, float64 (stored to 3 decimals in the
  in-house files, so a_i lands on a 0.001 muB grid).

## 2. Descriptor: SOAP, L2-normalised

dscribe `SOAP(periodic=True, r_cut=6.0, n_max=8, l_max=3)`, all other dscribe defaults
(sigma 1.0, crossover on, no averaging).  Centres = the pool sites.  Rows are computed in
float64, L2-normalised with `sklearn.preprocessing.normalize`, stored float32.

**Species basis.**  Each structure is described in its own basis `species = sorted(set(symbols))`
("per-compound basis", decision 2026-08-20).  The descriptor blocks of species pairs absent from a
structure are exactly zero in a global basis, so after L2 normalisation the per-compound vector is
the global-basis vector with its zero blocks removed; the dot product of two sites is the sum over
their **shared pair blocks**.  This was verified bit-exact against the global basis (source gate V2,
max |delta| = 0.0).  The in-house runs used the fixed global basis `Co Cr Li Na O P V`
(`--species inhouse`); both give the same numbers.  Block widths for n_max 8, l_max 3 are 144
(same-species) and 256 (cross) features.

The 6 A cutoff is shorter than a MACE receptive field (2 layers x 5–6 A); kept deliberately
(decision 2026-09-04) and the measured sensitivity is small (§6).

## 3. Neighbours: exact kNN, then the adaptive band

Distance between L2-normalised vectors: `d_ij = sqrt(clip(2 - 2 x_i . x_j, 0))`, dot products
accumulated in float32 over shared blocks (disjoint chemistry -> dot exactly 0 -> d = sqrt 2, never
computed).  **No approximate search anywhere.**

1. **Candidates**: every other site of the pool, **including the other sites of the query's own
   frame**; only the query site itself is excluded (its distance set to +inf; an exact duplicate
   of the site in another frame stays at distance 0).  Rule since 2026-10-09 (in-house) /
   2026-10-11 (corpora).  Previous rule: all sites of the query frame excluded (`--exclude-same-frame`).
2. **Exact K = 100 nearest** neighbours are stored, ties broken by site id (stable sort on distance
   with columns in site-id order = lexicographic (distance, site_id)).  K was 20 until 2026-10-09;
   with K = 20 the adaptive band was capped for 30–86 % of in-house sites.  K is only the candidate
   cap of the band; at alpha 1.25 it caps ~0.05 % of in-house sites.
3. **Adaptive band ("twins")**: `d_ij <= alpha * d_i,m` with alpha = 1.25 and m = 5, i.e. all stored
   neighbours not farther than 1.25 times the distance to the 5th neighbour; `k_eff` is clipped to
   `[m, K] = [5, 100]`.  The radius is formed in float32 exactly as stored.  Typical k_eff: in-house
   means 14–23, MatPES medians 9–10, MP-MM medians 7–8 (MP-MM is dense in near-duplicate relaxation frames).

Why a data-adaptive band and not an absolute SOAP radius: within one thermal dataset the DFT force
difference between two same-element sites is ~1 eV/A at *every* SOAP distance 0.002–0.3 and for every
cutoff 2.5–6 A (Spearman rho(d, dF) = 0.00); no force-equivalence radius exists in thermal data.  SOAP
distance does resolve the moment (median |dm| of the 5 nearest twins vs random pairs: NVP +U 0.03 vs
0.46 muB).  The band therefore measures moment disagreement among the *nearest sampled* environments,
which is the resolution a model trained on the same data has (analysis 2026-10-11).

## 4. Site ambiguity

```
a_i      = median_{j in twins(i)} | m_i - m_j |          signed   (headline since 2026-10-28)
a_i^u    = median_{j in twins(i)} | |m_i| - |m_j| |      unsigned (option --moments unsigned)
```
Both are always computed and stored (`ambiguity_signed_muB`, `ambiguity_unsigned_muB`); the headline
column `ambiguity_muB` follows the `moments` setting.  Also stored: the median over the full K list
(`ambiguity_exact_k_signed_muB`), the twin-median moment, the number of same-frame twins, the
nearest-neighbour distance and the band radius.

Signed vs unsigned: the signed value counts FM/AFM ordering disagreements and sign flips of
time-reversed sublattices as ambiguity; the unsigned value ignores them.  For predicting the force
error that a *magnetic* model removes from a geometry-only model, the signed statistic correlated
better on the corpora (ordering ambiguity also costs forces there), which is why T8 is signed; the
unsigned branch (2026-10-01/02) is kept as an option.  `median` over the twins is robust to one
outlier twin; the mean over the twins was never used.

## 5. Dataset ambiguity: T_p

```
tail(a, q):  sort a descending; mass = n q; k = floor(mass)
             total = sum(a[:k]) + (mass - k) a[k]            fractional boundary observation
             H_q = total / sum(a)      T_q = total / mass     (T_q = A H_q / q)
T8 = T_0.08 = mean ambiguity of the most ambiguous 8 % of (filtered) sites        [muB]
```
`tail` is the source project's function verbatim (`scan_tail_descriptors_20261011.py`).  For figures
the top-8 % *sites* are flagged with the integer rule n_tail = ceil(0.08 n), threshold = n_tail-th
largest a_i (`tail_flag`); T8 itself uses the fractional mass.

**Frame-force filter.**  Query sites in frames with any |F_DFT| > 20 eV/A are excluded from the
statistic; their sites still serve as neighbours of other sites (`--fmax 20`, `--fmax none` disables).
Introduced 2026-10-12 for the MatPES periodic-table tables (sub-angstrom collision frames dominate
the force error) and part of the definition since 2026-10-28.  The in-house reference tables predate it
and were computed over all sites (`--fmax none`).

**Why the upper tail and not the mean.**  A = mean a_i is floor + tail: the floor (median a_i,
0.003–0.1 muB) is diffuse same-state moment noise that forces barely feel; the tail is the
other-state / other-ordering part that a magnetic model removes.  The floor can outweigh the tail
(NVP no-U: floor 0.10 muB on every V but 0.07 % other-state sites), so A mis-ranks datasets.  Over
q = 3..20 %, T_q was the descriptor with the highest and most stable Spearman correlation with the
force-MAE reduction of mMACE over MACE across in-house (5 populations), MatPES-PBE and
MatPES-r2SCAN (7, later 14 elements each); q = 8 % was fixed on 2026-10-12 (p-convergence scan) and
re-confirmed in the alpha-convergence scan of 2026-10-22.  Also reported for continuity: T20, A,
median (sigma), C = A / median (ill-conditioned when median -> 0, retired), P = fraction of sites
with a_i > 0.5 muB, and the legacy (sigma, p) verdict (p < 1 % fine; p <= 15 % and sigma <= 0.01 muB
down-weight; else needs a magnetic model).  Force errors never enter any statistic.

## 6. Sensitivity (measured, source project)

* Band rule, in-house (28 adaptive rules m 2–20 x alpha 1.1–2.0, and fixed k 5–100): the (sigma, p)
  verdict changes in 1 of 29 x 33 cases; p moves by <= 2 pp; alpha >= 1.5 fills the K = 100 list for
  most sites.  Corpora: A within −8…+25 % of baseline for every pool; fixed k 5…100 moves A by up to 50 %.
* rho(T8, MAE reduction) vs alpha (1.0–3.0) and fixed k (5–100), 2026-10-22: alpha 1.25 is on the
  plateau for in-house and both MatPES functionals (`alpha_rule_convergence_20261022`).
* Same-frame candidates vs cross-frame (2026-10-11): A −0.09…+0.04 muB; P 1–9 points lower; C rises
  because same-frame twins are near-duplicates (61–71 % MatPES, 84–94 % MP-MM of sites have one).
* Near-duplicates (MP-MM): 40–46 % of retained edges at SOAP d < 0.01 pin the median near 0 and inflate
  C 3–7x; T8 is insensitive to this floor, which is one more reason for the tail statistic.

## 7. Options and what they reproduce

| option | default | alternatives used before |
|---|---|---|
| `--species` | per-structure basis | `inhouse` = Co Cr Li Na O P V global basis (identical values) |
| `--r-cut --n-max --l-max` | 6.0 8 3 | never changed |
| `--k` | 100 | 20 (until 2026-10-09/11) |
| `--alpha --m-ref` | 1.25 5 | scans only |
| `--exclude-same-frame` | off | on (cross-frame rule, until 2026-10-09/11) |
| `--moments` | signed | unsigned (2026-10-01 … 10-27 headline for p / sigma) |
| `--tail` | 0.08 | 0.20 (T20, 2026-10-11), 0.10 |
| `--fmax` | 20 eV/A | none |
| `--max-abs-moment` | none (in-house) | 8 (MatPES / MP-MM) |

## 8. Implementation notes

* `soap.BlockSparseDescriptors.dots(start, stop, out)`: when every site shares one species basis
  (global basis, in-house) the dense float32 matrix is rebuilt in dscribe's native column order and
  one matmul is used; this reproduces the in-house reference descriptors, neighbour lists, distances,
  k_eff and a_i **bit for bit** (checked on Co3O4 50 K and N2VP 50 K, 2026-10-30).  With several bases
  (corpora) `X_q X_p^T` is accumulated block by block into the float32 `(n_query, n_pool)` buffer;
  the different float32 summation order moves distances by up to ~2e-4 near 0 and can reorder
  near-equidistant neighbours (observed effect on T8: 5e-6 muB).  `neighbors.exact_knn` then takes the
  stable argsort per row.  Memory: `query_chunk x n_pool x 4` bytes (4096 x 137k = 2.2 GB).
* `neighbors.npz` stores the exact K lists so that `rescore` can evaluate any alpha / m / fixed k
  without a new kNN (that is how the sensitivity scans were done).
* The settings hash (`Settings.version()`) covers SOAP, neighbour rule, moment convention and
  population rules; the dataset section (tail fraction, force filter, tau) does not change it because
  it is a post-processing choice on the same site table.
* All distances, k_eff and a_i are reproducible bit-for-bit from the stored descriptors; float32
  accumulation order can move a distance by ~1e-4 near 0 and swap equidistant duplicates (source
  gates V3–V5 quantify this as delta A ~ 2e-5 muB).
