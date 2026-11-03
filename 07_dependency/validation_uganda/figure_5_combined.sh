#!/bin/bash
# figure_5_combined.sh
#
# CESC + Uganda combined view of Figure 5 (leading-edge circos, both the
# outer-diamond and inner-tick variants). Does NOT touch the original
# CESC-only figures (figures/Figure_5_leading_edge_circos*.png/.pdf).
#
# Unlike every other figure, 10_/10b_figure_5_leading_edge_circos*.R draw
# with base-R graphics (circlize), not ggplot/patchwork, and each already
# renders BOTH clades of one cohort into a single huge (10800x25920 px)
# raster. Recomposing that at the pixel level (not re-running circlize) is
# the practical way to combine cohorts here: resize each cohort's already-
# rendered image, label it, and place the two side by side under one tag.
#
# Prerequisite: a pure Uganda-only render of each variant must exist first.
# The 07_dependency/10*.R scripts in this repo are the CESC-only originals
# (left untouched, per instructions); to (re)produce a Uganda-only render,
# temporarily add the same COHORT-argument branching used elsewhere in this
# project to a scratch copy of the script (do not edit 10_/10b_ in place),
# or restore it from an earlier parametrized version if one is kept
# elsewhere. This script only does the pixel-level combination step.
#
#   A) [CESC circos | Uganda circos]   (each already A7 | A9 internally)
#
# Usage:
#   figure_5_combined.sh <cesc_png> <uganda_png> <out_name>
# Example:
#   figure_5_combined.sh \
#     "$HOME/CESC_Network/figures/Figure_5_leading_edge_circos.png" \
#     "$HOME/CESC_Network/figures/validation_uganda/Figure_5_leading_edge_circos_uganda_ONLY.png" \
#     Figure_5_leading_edge_circos_uganda
set -e
CESC_SRC=$1
UGANDA_SRC=$2
OUT_NAME=$3
[[ -z "$OUT_NAME" ]] && { echo "usage: figure_5_combined.sh <cesc_png> <uganda_png> <out_name>"; exit 1; }

export MAGICK_CONFIGURE_PATH=/tmp/im_policy
mkdir -p /tmp/im_policy
cat > /tmp/im_policy/policy.xml << 'EOF'
<policymap>
  <policy domain="resource" name="memory" value="4GiB"/>
  <policy domain="resource" name="map" value="8GiB"/>
  <policy domain="resource" name="width" value="60KP"/>
  <policy domain="resource" name="height" value="60KP"/>
  <policy domain="resource" name="area" value="1GP"/>
  <policy domain="resource" name="disk" value="8GiB"/>
  <policy domain="resource" name="time" value="600"/>
</policymap>
EOF

OUT="$HOME/CESC_Network/figures/validation_uganda"
mkdir -p "$OUT"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

convert "$CESC_SRC"   -resize 2400x "$TMP/cesc.png"
convert "$UGANDA_SRC" -resize 2400x "$TMP/uganda.png"
convert "$TMP/cesc.png"   -gravity North -background white -splice 0x140 \
  -gravity North -pointsize 90 -font DejaVu-Sans-Bold -fill black -annotate +0+30 "CESC" \
  "$TMP/cesc_lab.png"
convert "$TMP/uganda.png" -gravity North -background white -splice 0x140 \
  -gravity North -pointsize 90 -font DejaVu-Sans-Bold -fill black -annotate +0+30 "Uganda" \
  "$TMP/uganda_lab.png"
convert "$TMP/cesc_lab.png" "$TMP/uganda_lab.png" +append "$TMP/combined.png"
convert "$TMP/combined.png" -gravity NorthWest -pointsize 110 -font DejaVu-Sans-Bold \
  -fill black -annotate +20+10 "A" "$OUT/${OUT_NAME}.png"
convert "$OUT/${OUT_NAME}.png" "$OUT/${OUT_NAME}.pdf"

echo "[OK] saved: $OUT/${OUT_NAME}.png / .pdf"
