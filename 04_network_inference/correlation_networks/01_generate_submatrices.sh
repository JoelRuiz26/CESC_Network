#!/bin/bash
# Genera 100 submatrices de A7 (65->60, leave-5-out) y 100 de A9 (202->60, leave-142-out)
# usando las matrices de 03_preprocessing (15,282 genes, solo tumor primario).
set -e

REPO="${HOME}/CESC_Network"
MATDIR="${REPO}/03_preprocessing"          # matrices 2_7 / 2_8, 15,282 genes
BASE="${REPO}/04_network_inference/correlation_networks"
N_REPS=100
mkdir -p "${BASE}/A7/submatrices" "${BASE}/A9/submatrices"

subsample_matrix() {
    local infile=$1 outfile=$2 n_total_cols=$3 n_remove=$4 seed=$5 label=$6 rep=$7
    local repstr=$(printf '%02d' ${rep})

    local remove_cols=$(python3 -c "
import random
random.seed(${seed})
cols = random.sample(range(2, ${n_total_cols}+1), ${n_remove})
cols.sort()
print(' '.join(map(str, cols)))
")
    local keep_cols=$(python3 -c "
remove = set(map(int, '${remove_cols}'.split()))
print(','.join(str(i) for i in range(1, ${n_total_cols}+1) if i not in remove))
")
    cut -f "${keep_cols}" "${infile}" > "${outfile}"
    local out_n=$(head -1 "${outfile}" | awk '{print NF-1}')
    echo "[${label} r${repstr}] seed=${seed} | $((n_total_cols-1)) -> ${out_n} muestras"
}

echo "--- A7 (65 -> 60 muestras, leave-5-out) ---"
for rep in $(seq 1 ${N_REPS}); do
    subsample_matrix \
        "${MATDIR}/2_7_Full_counts_A7_annot.tsv" \
        "${BASE}/A7/submatrices/2_7_Full_counts_A7_r$(printf '%02d' ${rep}).tsv" \
        66 5 ${rep} "A7" ${rep}
done

echo "--- A9 (202 -> 60 muestras, leave-142-out) ---"
for rep in $(seq 1 ${N_REPS}); do
    subsample_matrix \
        "${MATDIR}/2_8_Full_counts_A9_annot.tsv" \
        "${BASE}/A9/submatrices/2_8_Full_counts_A9_r$(printf '%02d' ${rep}).tsv" \
        203 142 $((rep + 1000)) "A9" ${rep}
done

echo "200 submatrices listas (100xA7 + 100xA9) en ${BASE}/{A7,A9}/submatrices/"
