# 6_4_hubgene_essentiality_table.tsv

Lista completa de genes hub de cada red (HPV-A7, HPV-A9) con su esencialidad DepMap.

| Columna | Significado |
|---|---|
| `network` | Red del gen hub (HPV-A7 / HPV-A9). |
| `gene` | Gen hub (top 10 % por grado). |
| `depmap_dependency` | Dependencia CRISPR media (más negativa = más dependiente; `NA` si no hay dato). |
| `essential` | 1 si dependencia ≤ −0.5. |
| `pan_essential` | 1 si es esencial universal (flag DepMap). |
| `important` | 1 si esencial **y no** pan-esencial. |

## Conteos

| Red | Total hubs | essential | pan_essential | important |
|---|---|---|---|---|
| HPV-A7 | 1609 | 265 | 296 | 23 |
| HPV-A9 | 1586 | 337 | 391 | 25 |
