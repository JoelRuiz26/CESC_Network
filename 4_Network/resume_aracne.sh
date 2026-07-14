#!/bin/bash
# Resume ARACNe-AP: finish A7 at 7 threads (already has 17 bootstraps),
# then run A9 fully at 9 threads. All genes as regulators, 100 bootstraps total each.

source /home/joelr/miniconda3/etc/profile.d/conda.sh
conda activate aracne

JAR=/home/joelr/ARACNe-AP/dist/aracne.jar
XMX=-Xmx10G
BASE=/home/joelr/CESC_Network
MAT_DIR=$BASE/2_Prepro_TCGA_GTEx
OUT=$BASE/4_Network
PVAL=1E-8
NBOOT=100

# run_network <tag> <matrix> <outdir> <tfs> <threads> <start_seed>
run_network () {
  local tag=$1 matrix=$2 outdir=$3 tfs=$4 threads=$5 start=$6
  echo "############################################################"
  echo "[$(date '+%F %T')] START/RESUME network $tag  (threads=$threads, seeds $start..$NBOOT)"
  echo "############################################################"

  # Step 1: MI threshold only if missing
  if ls "$outdir"/miThreshold_*.txt >/dev/null 2>&1; then
    echo "[$(date '+%F %T')] $tag: MI threshold already present, skipping"
  else
    echo "[$(date '+%F %T')] $tag: calculating MI threshold"
    java $XMX -jar $JAR -e "$matrix" -o "$outdir" --tfs "$tfs" \
         --pvalue $PVAL --seed 1 --calculateThreshold \
         >> "$outdir/threshold.log" 2>&1
  fi

  # Step 2: bootstraps from start_seed..NBOOT, N threads each, sequential
  for i in $(seq "$start" $NBOOT); do
    java $XMX -jar $JAR -e "$matrix" -o "$outdir" --tfs "$tfs" \
         --pvalue $PVAL --seed $i --threads $threads \
         >> "$outdir/bootstrap.log" 2>> "$outdir/bootstrap_errors.log"
    echo "[$(date '+%F %T')] $tag: bootstrap $i/$NBOOT done"
  done

  # Step 3: consolidate
  echo "[$(date '+%F %T')] $tag: consolidating"
  java $XMX -jar $JAR -o "$outdir" --consolidate \
       >> "$outdir/consolidate.log" 2>&1

  local nfiles=$(ls "$outdir"/bootstrapNetwork_*.txt 2>/dev/null | wc -l)
  local edges=$(( $(wc -l < "$outdir/network.txt") - 1 ))
  echo "[$(date '+%F %T')] DONE network $tag -> $outdir/network.txt ($edges edges, from $nfiles bootstraps)"
}

# A7: already has 17 complete bootstraps -> resume at seed 18, 7 threads
run_network A7 "$MAT_DIR/2_7_Full_counts_A7_annot.tsv" "$OUT/A7" "$OUT/A7/genes_A7.txt" 7 18

# A9: fresh, full 100 bootstraps at 9 threads
run_network A9 "$MAT_DIR/2_8_Full_counts_A9_annot.tsv" "$OUT/A9" "$OUT/A9/genes_A9.txt" 9 1

echo "[$(date '+%F %T')] ALL NETWORKS COMPLETE"
