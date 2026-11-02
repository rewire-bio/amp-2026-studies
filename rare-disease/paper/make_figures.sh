#!/usr/bin/env bash
# Make the paper's figures in figures/: the charts from ../figures/rd_charts.py as PDF, and the Mermaid diagrams in
# diagrams/*.mmd as PDF and PNG.
#
#   bash make_figures.sh            # charts and diagrams
#   bash make_figures.sh charts     # charts only
#   bash make_figures.sh diagrams   # diagrams only
#
# Charts: charts_pdf.py runs rd_charts.py unchanged and saves a PDF next to each PNG in ../figures/out/ (not
# committed); the ones the paper uses are copied to figures/, and their PNGs to ../figures/ for COMPANION.md.
# Diagrams: mmdc with diagrams/config.json (neutral theme) and -f (page fitted to the diagram). mmdc's bundled Chrome
# is not installed, so diagrams/puppeteer.json points it at Google Chrome; set MMDC_PUPPETEER to override.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
STUDY="$(dirname "$HERE")"
OUT="$HERE/figures"
WHAT="${1:-all}"
mkdir -p "$OUT"
CHARTS="rd_by_type rd_combined rd_gains rd_stricter rd_indel_combined rd_cancer_by_type"

charts() {
  for s in snv indel cancer; do
    uv run -q --with matplotlib python "$HERE/charts_pdf.py" --study "$s"
  done
  for n in $CHARTS; do
    cp "$STUDY/figures/out/$n.pdf" "$OUT/$n.pdf"
    cp "$STUDY/figures/out/$n.png" "$STUDY/figures/$n.png"
    echo "wrote $OUT/$n.pdf"
  done
}

diagrams() {
  local pconf="${MMDC_PUPPETEER:-$HERE/diagrams/puppeteer.json}"
  for mmd in "$HERE"/diagrams/*.mmd; do
    for ext in pdf png; do
      out="$OUT/diagram_$(basename "$mmd" .mmd).$ext"
      if [[ $ext == png ]]; then
        mmdc -q -i "$mmd" -o "$out" -c "$HERE/diagrams/config.json" -p "$pconf" -b white -s 2
      else
        mmdc -q -i "$mmd" -o "$out" -c "$HERE/diagrams/config.json" -p "$pconf" -f -b white
      fi
      echo "wrote $out"
    done
  done
}

case "$WHAT" in
  all) charts; diagrams ;;
  charts) charts ;;
  diagrams) diagrams ;;
  *) echo "usage: make_figures.sh [all|charts|diagrams]" >&2; exit 2 ;;
esac
