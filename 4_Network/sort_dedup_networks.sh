#!/bin/bash
# Sort each network by MI descending and remove duplicate undirected pairs (A-B / B-A),
# keeping the first occurrence (highest MI, since sorted descending).
# Input : <DIR>/<NET>/network_<NET>_AP.txt   (cols: Regulator  Target  MI  pvalue)
# Output: <DIR>/<NET>/network_<NET>_AP_sort_unique.txt  (same directory)

DIR="/home/joelr/CESC_Network/4_Network"

for NET in A7 A9; do
    INPUT="$DIR/$NET/network_${NET}_AP.txt"
    OUTPUT="$DIR/$NET/network_${NET}_AP_sort_unique.txt"

    echo "[$(date '+%H:%M:%S')] Sorting $NET ..."
    SORTED=$(mktemp)
    head -1 "$INPUT" > "$SORTED"
    tail -n +2 "$INPUT" | sort -t$'\t' -k3,3rn >> "$SORTED"

    echo "[$(date '+%H:%M:%S')] Deduplicating $NET ..."
    awk 'BEGIN { FS=OFS="\t" }
         NR == 1 { print; next }
         {
             key = ($1 < $2) ? $1 "\t" $2 : $2 "\t" $1
             if (!seen[key]++) print
         }' "$SORTED" > "$OUTPUT"

    rm -f "$SORTED"

    IN_ROWS=$(tail -n +2 "$INPUT"  | wc -l)
    OUT_ROWS=$(tail -n +2 "$OUTPUT" | wc -l)
    echo "  Input: $IN_ROWS  ->  After sort+dedup: $OUT_ROWS  (removed $(( IN_ROWS - OUT_ROWS )))"
done

echo "Done. Files written:"
ls -lh "$DIR"/*/network_*_AP_sort_unique.txt
