#!/bin/bash
# Uso: bash run.sh <matriz.tsv>
# Ejemplo: bash run.sh 2_7_Full_counts_A7_annot_pvalue.tsv | tee -a salida &

set -euo pipefail

partools="$(pwd)/../parallel"

# ---------- Args ----------
ftsv="${1:-}"
if [[ -z "$ftsv" || ! -f "$ftsv" ]]; then
  echo "[run] ERROR: Debes pasar un .tsv válido. Ej: bash run.sh 2_7_Full_counts_A7_annot_pvalue.tsv" >&2
  exit 15
fi

echo "Processing MI calculations for: $ftsv ..."
nom="$(basename "$ftsv" .tsv)"

# ---------- Arranque en limpio ----------
rm -rf adj log mat.adj 2>/dev/null || true

# ---------- node.list ----------
# 1a columna completa (incluye encabezado, p.ej. 'Gene')
awk '{print $1}' "$ftsv" > node.list
cname="$(head -1 node.list)"
echo "Column Index Name: $cname"

# ---------- Parámetros ----------
CORES=70        # fijos
ALPHA="1e-2"    # menos estricto para recuperar densidad; ajusta a 1e-7 si quieres
EPS="0.1"         # 1 = SIN DPI (diagnóstico). Luego puedes poner 0.2
BOOTS="0"       # sin bootstraps (mínimos cambios)

# ---------- ARACNE paralelizado ----------
SECONDS=0
python "${partools}/aracne-par.py" "$ftsv" node.list "$cname" $CORES "$ALPHA" "$EPS" "$BOOTS" &> aracne.log
echo "ARACNe time: $(echo "$SECONDS/60" | bc -l) minutes."

# ---------- Unir .adj en mat.adj ----------
SECONDS=0
echo "Parameters to join: $nom node.list $cname"
python "${partools}/joinadj.py" "$nom" node.list "$cname" | tee -a salida
echo "join ADJ matriz time: $(echo "$SECONDS/60" | bc -l) minutes."

# ---------- Validar y mover mat.adj ----------
if [[ ! -s adj/mat.adj ]]; then
  echo "[run] ERROR: adj/mat.adj no existe o está vacío. Revisa 'aracne.log' y el contenido de 'adj/'." >&2
  exit 2
fi

echo "Moving adjacency matrix"
mv -f adj/mat.adj .

# ---------- SIF ----------
SECONDS=0
python "${partools}/adj2sif.py" > "${nom}.sif"
echo "Creating SIF: $(echo "$SECONDS/60" | bc -l) minutes."

# ---------- Sort ----------
SECONDS=0
sort -r -k3,3 "${nom}.sif" > "${nom}.sort"
echo "Sorting: $(echo "$SECONDS/60" | bc -l) minutes."

# ---------- Conteo ----------
echo "[run] Edges en SIF:  $(wc -l < "${nom}.sif")"
echo "[run] Edges en SORT: $(wc -l < "${nom}.sort")"

# ---------- No limpiar durante depuración ----------
# rm -rf adj log mat.adj node.list
