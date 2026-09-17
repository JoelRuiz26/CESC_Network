#!/bin/bash
# ARACNe-AP full run for CESC networks A7 and A9
# All genes as regulators, 100 bootstraps, 9 threads, run sequentially.
# ARACNe-AP is installed outside this repository; point CONDA_SH and JAR at
# your own installation, or drop the conda lines if aracne.jar is already on
# a reachable path.
CONDA_SH=${CONDA_SH:-$HOME/miniconda3/etc/profile.d/conda.sh}
source "$CONDA_SH"
conda activate aracne

JAR=${ARACNE_JAR:-$HOME/ARACNe-AP/dist/aracne.jar}
XMX=-Xmx10G
BASE=${BASE:-$HOME/CESC_Network}
MAT_DIR=$BASE/03_preprocessing
OUT=$BASE/04_network_inference
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
