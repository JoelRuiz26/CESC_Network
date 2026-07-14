#!/bin/bash
# Resume A9 ARACNe run after interruption.
# 90 complete bootstraps present -> continue seeds 91..100 with 8 threads, then consolidate.

source /home/joelr/miniconda3/etc/profile.d/conda.sh
conda activate aracne

JAR=/home/joelr/ARACNe-AP/dist/aracne.jar
XMX=-Xmx10G
BASE=/home/joelr/CESC_Network
MATRIX=$BASE/2_Prepro_TCGA_GTEx/2_8_Full_counts_A9_annot.tsv
OUT=$BASE/4_Network/A9
TFS=$OUT/genes_A9.txt
PVAL=1E-8
THREADS=8
NBOOT=100
START=91

echo "############################################################"
echo "[$(date '+%F %T')] RESUME A9  (threads=$THREADS, seeds $START..$NBOOT)"
echo "############################################################"

# MI threshold already present -> skip; recompute only if missing
if ! ls "$OUT"/miThreshold_*.txt >/dev/null 2>&1; then
  echo "[$(date '+%F %T')] A9: MI threshold missing, calculating"
  java $XMX -jar $JAR -e "$MATRIX" -o "$OUT" --tfs "$TFS" \
       --pvalue $PVAL --seed 1 --calculateThreshold >> "$OUT/threshold.log" 2>&1
else
  echo "[$(date '+%F %T')] A9: MI threshold present, skipping"
fi

# Remaining bootstraps, 8 threads each, sequential
for i in $(seq $START $NBOOT); do
  java $XMX -jar $JAR -e "$MATRIX" -o "$OUT" --tfs "$TFS" \
       --pvalue $PVAL --seed $i --threads $THREADS \
       >> "$OUT/bootstrap.log" 2>> "$OUT/bootstrap_errors.log"
  echo "[$(date '+%F %T')] A9: bootstrap $i/$NBOOT done"
done

# Consolidate
echo "[$(date '+%F %T')] A9: consolidating"
java $XMX -jar $JAR -o "$OUT" --consolidate >> "$OUT/consolidate.log" 2>&1

nfiles=$(ls "$OUT"/bootstrapNetwork_*.txt 2>/dev/null | wc -l)
edges=$(( $(wc -l < "$OUT/network.txt") - 1 ))
echo "[$(date '+%F %T')] DONE A9 -> $OUT/network.txt ($edges edges, from $nfiles bootstraps)"
echo "[$(date '+%F %T')] ALL NETWORKS COMPLETE"
