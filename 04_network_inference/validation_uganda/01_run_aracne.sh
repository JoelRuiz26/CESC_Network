#!/bin/bash
# ARACNe-AP full run for the Uganda (HTMCP-CC) validation networks A7 and A9.
# Mirrors ../01_run_aracne.sh exactly (same JAR, pvalue, threads, bootstraps),
# pointed at the Uganda preprocessed matrices instead of the CESC ones.
# All genes as regulators, 100 bootstraps, 7 threads, run sequentially.
CONDA_SH=${CONDA_SH:-$HOME/miniconda3/etc/profile.d/conda.sh}
source "$CONDA_SH"
conda activate aracne

JAR=${ARACNE_JAR:-$HOME/ARACNe-AP/dist/aracne.jar}
XMX=-Xmx10G
BASE=${BASE:-$HOME/CESC_Network}
MAT_DIR=$BASE/03_preprocessing/validation_uganda
OUT=$BASE/04_network_inference/validation_uganda
PVAL=1E-8
THREADS=7
NBOOT=100

mkdir -p "$OUT/A7" "$OUT/A9"

# Regulator lists (all genes), derived from the Uganda matrices themselves
tail -n +2 "$MAT_DIR/2_7_Full_counts_A7_annot.tsv" | cut -f1 > "$OUT/A7/genes_A7.txt"
tail -n +2 "$MAT_DIR/2_8_Full_counts_A9_annot.tsv" | cut -f1 > "$OUT/A9/genes_A9.txt"

run_network () {
  local tag=$1 matrix=$2 outdir=$3 tfs=$4
  echo "############################################################"
  echo "[$(date '+%F %T')] START network $tag"
  echo "############################################################"

  # Step 1: MI threshold
  echo "[$(date '+%F %T')] $tag: calculating MI threshold"
  java $XMX -jar $JAR -e "$matrix" -o "$outdir" --tfs "$tfs" \
       --pvalue $PVAL --seed 1 --calculateThreshold \
       >> "$outdir/threshold.log" 2>&1

  # Step 2: 100 bootstraps, seeds 1..100, 7 threads each, sequential
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

  # ARACNe writes the consolidated network as network.txt; rename with a
  # Uganda_ prefix so it can never be confused with the CESC network_<tag>_AP.txt
  # files (same tag, different cohort), and drop the 100 per-bootstrap
  # intermediates + the now-redundant network.txt (CESC's own dirs don't keep them).
  mv "$outdir/network.txt" "$outdir/Uganda_${tag}_AP.txt"
  rm -f "$outdir"/bootstrapNetwork_*.txt

  local edges=$(( $(wc -l < "$outdir/Uganda_${tag}_AP.txt") - 1 ))
  echo "[$(date '+%F %T')] DONE network $tag -> $outdir/Uganda_${tag}_AP.txt ($edges edges)"
}

run_network A7 "$MAT_DIR/2_7_Full_counts_A7_annot.tsv" "$OUT/A7" "$OUT/A7/genes_A7.txt"
run_network A9 "$MAT_DIR/2_8_Full_counts_A9_annot.tsv" "$OUT/A9" "$OUT/A9/genes_A9.txt"

echo "[$(date '+%F %T')] ALL UGANDA NETWORKS COMPLETE"
