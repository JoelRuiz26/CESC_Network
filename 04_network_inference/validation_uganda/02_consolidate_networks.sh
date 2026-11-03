#!/bin/bash
# Sort each Uganda network by MI descending and remove duplicate undirected pairs
# (A-B / B-A), keeping the first occurrence (highest MI, since sorted descending).
# Mirrors ../02_consolidate_networks.sh exactly, pointed at validation_uganda.
# Input : <DIR>/<NET>/Uganda_<NET>_AP.txt   (cols: Regulator  Target  MI  pvalue)
# Output: <DIR>/<NET>/Uganda_<NET>_AP_sort_unique.txt  (same directory)
# Uganda_ prefix keeps these from ever being confused with CESC's own
# network_<NET>_AP*.txt files one directory up, which use the same A7/A9 tags.

DIR="${HOME}/CESC_Network/04_network_inference/validation_uganda"

for NET in A7 A9; do
    INPUT="$DIR/$NET/Uganda_${NET}_AP.txt"
    OUTPUT="$DIR/$NET/Uganda_${NET}_AP_sort_unique.txt"

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
ls -lh "$DIR"/*/Uganda_*_AP_sort_unique.txt
