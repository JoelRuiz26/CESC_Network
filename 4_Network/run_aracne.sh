#!/bin/bash
# ARACNe-AP full run for CESC networks A7 and A9
# All genes as regulators, 100 bootstraps, 9 threads, run sequentially.
source /home/joelr/miniconda3/etc/profile.d/conda.sh
conda activate aracne

JAR=/home/joelr/ARACNe-AP/dist/aracne.jar
XMX=-Xmx10G
BASE=/home/joelr/CESC_Network
MAT_DIR=$BASE/2_Prepro_TCGA_GTEx
OUT=$BASE/4_Network
PVAL=1E-8
THREADS=7
NBOOT=100

run_network () {
  local tag=$1 matrix=$2 outdir=$3 tfs=$4
  echo "############################################################"
  echo "[$(date '+%F %T')] START network $tag"
  echo "############################################################"

  # Step 1: MI threshold (skip if already present from timing runs)
  echo "[$(date '+%F %T')] $tag: calculating MI threshold"
  java $XMX -jar $JAR -e "$matrix" -o "$outdir" --tfs "$tfs" \
       --pvalue $PVAL --seed 1 --calculateThreshold \
       >> "$outdir/threshold.log" 2>&1

  # Step 2: 100 bootstraps, seeds 1..100, 9 threads each, sequential
  for i in $(seq 1 $NBOOT); do
    java $XMX -jar $JAR -e "$matrix" -o "$outdir" --tfs "$tfs" \
         --pvalue $PVAL --seed $i --threads $THREADS \
         >> "$outdir/bootstrap.log" 2>> "$outdir/bootstrap_errors.log"
    echo "[$(date '+%F %T')] $tag: bootstrap $i/$NBOOT done"
  done

  # Step 3: consolidate
  echo "[$(date '+%F %T')] $tag: consolidating"
  java $XMX -jar $JAR -o "$outdir" --consolidate \
       >> "$outdir/consolidate.log" 2>&1

  local edges=$(( $(wc -l < "$outdir/network.txt") - 1 ))
  echo "[$(date '+%F %T')] DONE network $tag -> $outdir/network.txt ($edges edges)"
}

run_network A7 "$MAT_DIR/2_7_Full_counts_A7_annot.tsv" "$OUT/A7" "$OUT/A7/genes_A7.txt"
run_network A9 "$MAT_DIR/2_8_Full_counts_A9_annot.tsv" "$OUT/A9" "$OUT/A9/genes_A9.txt"

echo "[$(date '+%F %T')] ALL NETWORKS COMPLETE"
