# Correlation relevance networks — TCGA-CESC, HPV clades A7 and A9

Pearson and Spearman co-expression networks built on the same preprocessed
matrices and the same subsamples that the ARACNe analysis uses, so the three
inference methods can be compared on equal footing.

Method: relevance network (Butte, Tamayo, Slonim, Golub & Kohane, *PNAS*
97:12182–12186, 2000). All pairwise correlations are computed and everything
below a significance threshold is dropped; what survives is the network.

## Input

`03_preprocessing/2_7_Full_counts_A7_annot.tsv` and `2_8_..._A9_...`, i.e. the
output of `03_preprocessing/01_normalise_rnaseq.R`: **15,282 gene symbols**,
primary tumour samples only (A7 = 65, A9 = 202), GC/length + TMM normalised,
no ARSyN.

Nothing here renormalises. Subsamples are column cuts of those two files.

## Scripts

| script | what it does |
|---|---|
| `01_generate_submatrices.sh` | 100 leave-x-out subsamples per clade, n = 60 each (A7 65→60, A9 202→60), fixed seeds |
| `02_pearson_networks.R` | the 200 Pearson networks |
| `03_spearman_networks.R` | the 200 Spearman networks |
| `04_pairwise_metrics.R` | the four pairwise metrics over all 19,900 network pairs per method |
| `05_spearman_full_networks.R` | the two full-cohort Spearman networks |
| `06_qap_contrasts.R` | **superseded** — label permutation over the finished networks; kept for the record, see *Why the QAP was wrong* |
| `07_overlap_control.R` | regresses the patient-overlap confound out of the A9 contrast; supporting evidence |
| `08_figure_internal_validation.R` | the publication figure |
| `09_patient_label_null.R` | **the test the figures report** — permutes the clade label of the patients and rebuilds the networks |
| `10_pearson_full_networks.R` | the two full-cohort Pearson networks |
| `11_figure_null_distribution.R` | the permutation-null figure |

## Edge criterion (identical for both methods and both clades)

```
t     = r * sqrt(n - 2) / sqrt(1 - r^2),   df = n - 2
alpha = 0.05 / NP,   NP = 15282 * 15281 / 2 = 116,762,121 pairs
```

The threshold is a function of `n` and `NP` only — nothing is tuned to the
data — so every network's edge count is a result, not a setting. For the 200
subsample networks `n = 60` gives **|r| ≥ 0.7015**.

Per-network edge and node counts: `pairwise_metrics/network_stats_*.csv`.

## Full-cohort networks — `full_networks/`

`Gene1`, `Gene2`, `rho`, sorted by |rho| descending.

| method | clade | n | threshold | edges | connected nodes | % of 15,282 | mean degree | % negative |
|---|---|---|---|---|---|---|---|---|
| Spearman | A7 | 65 | 0.681088 | 91,278 | 8,387 | 54.9 | 21.8 | 6.9 |
| Spearman | A9 | 202 | 0.421280 | 1,064,581 | 14,813 | 96.9 | 143.7 | 18.6 |
| Pearson | A7 | 65 | 0.681088 | 178,529 | 9,022 | 59.0 | 39.6 | 23.3 |
| Pearson | A9 | 202 | 0.421280 | 1,248,435 | 14,811 | 96.9 | 168.6 | 22.8 |

Pearson files carry `r` where the Spearman ones carry `rho`; the threshold
formula is the same, so at equal `n` the two methods are cut at the same
number. Pearson keeps more edges in both clades and three times as many
negative ones — it is the raw linear correlation, with no rank transform to
damp the influence of a few extreme samples.

### Do not compare these two by edge count or density

A Bonferroni threshold on a correlation is a function of `n`, so the larger
cohort is cut at a lower point and will always look denser. Counting edges in
both clades at **matched** thresholds shows the effect reverses:

| threshold \|rho\| | A7 (n=65) | A9 (n=202) | A9/A7 |
|---|---|---|---|
| 0.4238 (A9's own cut) | 4,742,992 | 1,255,126 | 0.26 |
| 0.55 | 777,814 | 265,073 | 0.34 |
| 0.6844 (A7's own cut) | 95,781 | 60,862 | 0.64 |
| 0.80 | 16,326 | 12,884 | 0.79 |

At any common cut A7 has **more** edges than A9 — the opposite of what the
per-clade thresholds suggest. (Those counts were measured on the 17,481-gene
variant of the matrix; the ratios carry over.)

A common threshold is not a fix either. Under H0 the s.d. of `r` is about
`1/sqrt(n-1)`: 0.125 for A7 against 0.071 for A9, so the smaller cohort's
correlations are intrinsically more dispersed and its edge count is inflated
on its own.

**Valid between-clade comparisons** are either at a matched *number of edges*,
or on the 200 leave-x-out networks, which are all n = 60 and therefore share a
single null. These two files describe each clade **separately**.

## Pairwise metrics — `pairwise_metrics/`

`pairwise_metrics_long.csv.gz` — one row per unordered pair of networks per
method: 19,900 pairs × 3 methods = 59,700 rows. This is the input for the
violin panels.

| column | meaning |
|---|---|
| `method`, `method_key` | ARACNe / Pearson / Spearman |
| `group` | `intra-A7`, `intra-A9`, `inter A7-A9` |
| `Jaccard` | edge overlap — **similarity**, higher = more alike |
| `Spearman` | rho between the two degree sequences, over the fixed 15,282-gene universe with zeros kept so absent genes still count — **similarity** |
| `ARI` | adjusted Rand index of the two Louvain partitions on shared nodes — **similarity** |
| `KS` | two-sample KS statistic between the degree vectors — **distance**, lower = more alike |

`KS` is the D statistic used as an effect size, not its p-value: degrees are
not independent draws, so the p-value is meaningless at this many nodes.

`pairwise_metrics_summary.csv` — group means, the QAP p-value and Cohen's d
per method and metric:

| method | metric | intra-A7 | intra-A9 | inter | QAP p | d |
|---|---|---|---|---|---|---|
| Pearson | Jaccard | 0.661 | 0.382 | 0.204 | 0.0002 | −2.82 |
| Pearson | Spearman | 0.921 | 0.675 | 0.553 | 0.0002 | −2.66 |
| Pearson | ARI | 0.571 | 0.520 | 0.360 | 0.0002 | −3.58 |
| Pearson | KS | 0.018 | 0.027 | 0.102 | 0.0002 | +4.09 |
| Spearman | Jaccard | 0.699 | 0.401 | 0.300 | 0.0002 | −2.28 |
| Spearman | Spearman | 0.896 | 0.651 | 0.529 | 0.0002 | −2.75 |
| Spearman | ARI | 0.603 | 0.504 | 0.387 | 0.0002 | −2.95 |
| Spearman | KS | 0.015 | 0.026 | 0.099 | 0.0002 | +4.25 |

The QAP statistic is always `inter − intra`, so a clade effect is negative for
the three similarities and positive for KS D. `p = 0.0002` is the floor for
B = 5000 permutations.

### The test — `patient_null/patient_null_pvalues.csv`

Each intra group is contrasted against inter; `A7:A7` against `A9:A9` is
**not** admissible, because two A7 subsamples share 55–58 of 60 patients (A7
has only 65) while two A9 subsamples share 8–29. The statistic is
`mean(inter) − mean(intra)`, as before.

The null is built by permuting the clade label of the **patients**
(`09_patient_label_null.R`): the 267 tumours are split at random into a
pseudo-A7 of 65 and a pseudo-A9 of 202 — the real sizes — and then the
leave-x-out subsampling, the networks and the four metrics are rebuilt from
scratch. Overlap structure is preserved exactly; only clade identity is
destroyed. B = 100, 20 + 20 networks per permutation, so the P floor is
1/101 = 0.0099. `z` is `(observed − null mean) / null s.d.`

| metric | A7:A7 vs A7:A9 | | A9:A9 vs A7:A9 | |
|---|---|---|---|---|
| | *P* | *z* | *P* | *z* |
| Jaccard | 0.33 | −0.55 | 0.14 | −1.16 |
| Degree ρ | **0.0099** | −4.33 | **0.0099** | −3.57 |
| ARI | 0.55 | +0.15 | **0.0099** | −2.80 |
| KS D | **0.0099** | +4.37 | **0.0099** | +4.48 |

`0.0099` is the floor: no permutation out of 100 reached the observation.

**Edge identity does not carry a clade signal.** Jaccard fails in both clades,
and the ARI of A7 fails too. With clades assigned at random the intra-A7
Jaccard averages **0.7197**, slightly *above* the observed 0.6953: that 0.70 is
the 55–58 shared patients, not the clade.

**Topology does.** Degree-sequence correlation and the KS distance between
degree distributions survive in both clades, and module agreement survives in
A9. The clade changes which genes are hubs and how the degree mass is
distributed, not the exact list of surviving edges — which is the part the
subsampling design inflates.

### Why the QAP was wrong — `figure/qap_contrasts.csv`

The earlier test (`06_qap_contrasts.R`) permuted the clade label of the 200
finished networks. That assumes the 200 are exchangeable under H0. They are
not: the overlap structure above is baked into them, and permuting their
labels destroys it. Its null was therefore far too flat — for the Jaccard of
A7, mean 0.00004 ± 0.011 against a correct null of −0.3888 ± 0.0159 — so every
contrast hit the 1/(B+1) floor regardless of the data. It reported
*P* = 5 × 10⁻⁵ and *d* = −14.8 for the Jaccard contrast that the correct null
cannot distinguish from chance. The file is kept so the difference is on the
record; do not quote its P values.

### The overlap control — `figure/overlap_control.csv`

Inter-clade pairs share **zero** patients by construction: A7 and A9 are
disjoint patient sets. So every intra-vs-inter contrast starts with an overlap
advantage, and the permutation test does not remove it — permuting labels
destroys clade structure but not the fact that same-clade networks were drawn
from one patient pool.

Within A9 the shared count varies enough (8–29 of 60) to regress it out and
extrapolate to zero shared patients. It survives on all four metrics:

| metric | intra-A9 observed | at 0 shared [95% CI] | inter | survives |
|---|---|---|---|---|
| Jaccard | 0.4014 | 0.3325 [0.3261, 0.3390] | 0.3003 | yes |
| Degree ρ | 0.6508 | 0.5767 [0.5732, 0.5802] | 0.5289 | yes |
| ARI | 0.5037 | 0.4627 [0.4545, 0.4708] | 0.3871 | yes |
| KS D | 0.0262 | 0.0299 [0.0273, 0.0325] | 0.0990 | yes |

**A7 cannot be controlled this way.** Its pairs share 55–58 of 60, so
extrapolating to zero is a 55-patient extrapolation outside the observed
range. Its slope is reported instead (0.0297 Jaccard per shared patient for
Spearman).

This regression is now supporting evidence rather than the main control: the
patient-label null above handles the same confound for both clades at once,
without extrapolating. Where the two can be compared they agree — A9 survives
on all four metrics under both.

For contrast, the same control applied to the ARACNe networks
(`05_topology`/leave-x-out analysis) removes the effect entirely: intra-A9
Jaccard 0.0486 observed, 0.0303 extrapolated, against 0.0324 inter — a ratio
of 0.93. The clade signal is present at the co-expression level and is lost to
ARACNe's sparsification.

### Never read A7:A7 against A9:A9

`A7:A7` exceeds `A9:A9` on every metric and that is the sampling design, not
biology — see the two sections above. The valid contrasts are each intra group
against inter, controlled for overlap where the design permits it.

## Figure — `figure/`

Two figures, one question each.

`Figure_Spearman_internal_validation.pdf` — **descriptive**: how alike are
networks within a clade and between clades. Violins of every pair among the
100 + 100 networks, with the two brackets and their P values.

`Figure_Spearman_permutation_null.pdf` — **inferential**: where those P values
come from. One histogram per contrast — the 100 patient-label permutations —
with the observed statistic drawn on it. A line inside its own histogram is a
contrast the subsampling design produces by itself; panel A is the Jaccard,
and both its lines sit inside.

An earlier version put the null as a band inside the violins. One panel, two
questions, and the band read as a second boxplot.

Palette `#2E6DA4` / `#CC0000` / `#8C8C8C`, validated for the lightness band,
CVD separation (worst adjacent pair ΔE 14.5 deuteranopia, 24.4 tritanopia),
normal-vision ΔE 24.4 and contrast against the surface. The third slot is a
deliberate neutral — `A7 : A9` is the reference group, not a third category —
and group identity is carried by the x axis, never by colour alone. No legend,
because the groups *are* the x categories.

## Known limitation of the preprocessing

Inherited from `03_preprocessing/01_normalise_rnaseq.R`. EDASeq's
`withinLaneNormalization` always normalises `counts(x)`, never `normCounts(x)`,
so the two-step call

```r
gcFull <- withinLaneNormalization(mydataEDA, "gc_percent", which = "full")
lFull  <- withinLaneNormalization(gcFull,    "length",     which = "full")
```

does **not** compose: the second call restarts from the raw counts and the GC
correction is discarded. What reaches TMM is length-corrected but not
GC-corrected. It applies equally to both clades, so it does not bias the
A7-vs-A9 contrast, but it should be reported as a limitation.

## References

- Butte, Tamayo, Slonim, Golub & Kohane. *PNAS* 97(22):12182–12186 (2000) — relevance networks.
- Song, Langfelder & Horvath. *BMC Bioinformatics* 13:328 (2012) — correlation vs mutual information as a co-expression measure.
- Vandenbon. *PLoS ONE* 17(1):e0263344 (2022) — Pearson vs Spearman as a function of dataset size.
- van Dam, Võsa, van der Graaf, Franke & de Magalhães. *Briefings in Bioinformatics* 19(4):575–592 (2018) — review.
