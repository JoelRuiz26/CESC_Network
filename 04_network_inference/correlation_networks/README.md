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
| `04_pairwise_metrics.R` | the four pairwise metrics over all 19,900 network pairs per method, plus the QAP test |
| `05_spearman_full_networks.R` | the two full-cohort Spearman networks |

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

| clade | n | threshold \|rho\| | edges | connected nodes | % of 15,282 | mean degree | % negative |
|---|---|---|---|---|---|---|---|
| A7 | 65 | 0.681088 | 91,278 | 8,387 | 54.9 | 21.8 | 6.9 |
| A9 | 202 | 0.421280 | 1,064,581 | 14,813 | 96.9 | 143.7 | 18.6 |

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

### Significance test

The clade label of each of the 200 networks is permuted B = 5000 times keeping
the 100/100 split, and the group means are recomputed from the same fixed
pairwise matrix. **The permutation unit is the network, not the pair** — the
19,900 pairwise values are not independent, and treating them as such is the
pseudo-replication the test exists to avoid.

### Read intra-A7 vs intra-A9 with care

`intra-A7` exceeds `intra-A9` on every metric, and that is a sampling artefact,
not biology: A7 has 65 patients and each subsample takes 60, so any two A7
subsamples share at least 55 of 60 patients (92%). A9 has 202, so two A9
subsamples share only 12–29. The valid contrast is **intra vs inter**.

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
