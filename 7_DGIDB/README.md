# 7_DGIDB — Drug Proximity to DepMap Interest Genes

Network-based drug repurposing using the closest proximity metric (Guney et al. 2016, *Nat Commun*).  
Interest genes = DEG ∩ Hub genes from DepMap analysis (A7: 604 genes, A9: 422 genes).  
Networks: HPV-A7 and HPV-A9 co-expression networks (17,481 nodes, ~595,000 edges each).  
Drugs tested: 2,814 FDA-approved drugs from DGIdb.

---

## Scripts

| Script | Descripción |
|--------|-------------|
| `6_1_DGIBD_metrics.R` | Descarga y filtra interacciones DGIdb FDA-approved; construye grafo bipartito droga–gen |
| `6_2_bipartite_cosmos.R` | Visualización bipartita del grafo droga–gen (tema cosmos y blanco) |
| `6_3_proximity_interest_genes.R` | Calcula proximidad de red (d_c), z-score por permutaciones, selectividad y significancia |
| `6_4_proximity_plots.R` | Genera las 7 figuras de publicación a partir de los resultados de proximidad |

---

## Método de proximidad

**Closest proximity** (Guney 2016):

```
d_c(T, S) = (1/|T|) * Σ_{t∈T} min_{s∈S} dist(t, s)
```

**Selectividad:**

```
S = d_c(T, non-interest) / d_c(T, interest)    [S > 1 → droga preferentemente cerca de genes de interés]
```

**Significancia:** z-score por permutación degree-matched de targets (500 permutaciones).  
Umbral: **z < −1.96** (convención Guney 2016; BH también reportado como referencia).

> Nota: BH sobre 2,814 drogas exige z < −4.1 para significancia, lo cual es excesivamente estricto para drug repurposing. La literatura (Guney 2016, Cheng 2019, Zhou 2020) usa umbrales de z directamente.

---

## Outputs principales

| Archivo | Contenido |
|---------|-----------|
| `6_3_proximity_A7.tsv` | Scores por droga en red A7 |
| `6_3_proximity_A9.tsv` | Scores por droga en red A9 |
| `6_3_proximity_combined.tsv` | Tabla combinada A7 + A9 |
| `6_3_proximity_objects.rds` | Objetos R (resultados + drug_map + interest genes) |
| `6_4_significant_drugs_summary.tsv` | Solo drogas con z < −1.96 |

---

## Figuras (6_4_P*.pdf / .png)

### P1 · Volcano — z-score vs selectividad
Vista panorámica de las 2,814 drogas. Eje X: z-score (más negativo = más cercana a genes de interés que por azar). Eje Y: selectividad (d_c no-interés / d_c interés; >1 = droga prefiere genes de interés sobre el resto de la red). El cuadrante ideal es **izquierda + arriba**. El tamaño del punto refleja el número de targets en la red.

### P2 · Lollipop top 20
Las 20 drogas con z-score más negativo por cohorte, rankeadas de mayor a menor proximidad. Asterisco rojo = también pasa BH < 0.05. Es la figura más directa para identificar candidatos y comparar A7 vs A9.

### P3 · Scatter d_c(interest) vs d_c(no-interest)
Cada punto es una droga. Eje X: distancia media a genes de interés. Eje Y: distancia media al resto de la red. Puntos **sobre la diagonal** = droga más cerca de genes de interés que del fondo (selectiva). Color = z-score. Las drogas ideales están abajo-izquierda y sobre la diagonal.

### P4 · Violin selectividad
Compara la distribución de selectividad entre drogas significativas (z < −1.96) y no significativas. Si las significativas tienen mediana > 1 y distribución desplazada hacia arriba, confirma que la proximidad estadística va acompañada de especificidad biológica real.

### P5 · Cross-cohort A7 vs A9
Cada punto es una droga presente en ambas redes. Eje X = z-score en A7, eje Y = z-score en A9. El cuadrante sombreado (ambos z < −1.96) contiene drogas candidatas en **ambas cohortes HPV**. Las que solo caen en un eje son específicas de ese clade. Útil para priorizar candidatos transversales.

### P6 · Bubble top 30
Top 30 por z-score, con color y segmento codificando selectividad (azul oscuro = baja, naranja = alta). Combina ranking, selectividad y número de targets en una sola figura. La más informativa para publicación.

### P7 · Null distributions (ridge plot)
Para las top 5 drogas significativas de cada cohorte: distribución nula gaussiana (500 permutaciones, azul) vs d_c observado (línea roja). Cuanto más a la izquierda cae la línea roja respecto a la campana, más extrema y confiable es la proximidad. Validación visual del z-score.

---

## Resultados (z < −1.96)

| Red | z < −1.96 | z < −2 | z < −3 | BH < 0.05 |
|-----|-----------|--------|--------|-----------|
| A7  | 39        | 38     | 4      | 2         |
| A9  | 25        | 24     | 3      | 0         |

**Top candidatos A7 (z < −3):** TRIPTORELIN (z=−4.14), ANAGRELIDE HCl (z=−3.96), ESTRONE (z=−3.40), OXALIPLATIN (z=−3.13)  
**Top candidatos A9 (z < −3):** GEMCITABINE (z=−3.60), CYTARABINE (z=−3.24), CLOFARABINE (z=−3.06)

---

## Referencias

- Guney et al. 2016 — https://doi.org/10.1038/ncomms10331
- Cheng et al. 2019 — https://doi.org/10.1038/s41467-019-10744-6
- Zhou et al. 2020 — https://doi.org/10.1038/s41421-020-0153-3
