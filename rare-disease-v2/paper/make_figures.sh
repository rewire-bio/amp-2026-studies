#!/usr/bin/env bash
# Make the paper's figures in figures/: the charts from ../figures/rd_charts.py as PDF, and the Mermaid diagrams in
# diagrams/*.mmd as PDF.
#
#   bash make_figures.sh            # charts and diagrams
#   bash make_figures.sh charts     # charts only
#   bash make_figures.sh diagrams   # diagrams only
#
# Charts: rd_charts.py writes SVG and PNG to figures/charts/ (it is not edited here); each SVG is converted to
# figures/<name>.pdf with rsvg-convert (cairosvg and Inkscape are not installed). The SVGs keep text as text
# (svg.fonttype none) and name DejaVu Sans first, which is not installed here, so rsvg-convert falls back to the next
# family in the list (Lucida Grande or Verdana on macOS).
# evo2_gains is drawn from whatever Evo 2 results exist: rerun this after the Evo 2 run before using it.
#
# Diagrams: mmdc with diagrams/config.json (neutral theme) and -f (page fitted to the diagram); PDF output is white.
# mmdc's bundled Chrome is not installed, so diagrams/puppeteer.json points it at Google Chrome; set
# MMDC_PUPPETEER to another puppeteer config to override.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
STUDY="$(dirname "$HERE")"
OUT="$HERE/figures"
WHAT="${1:-all}"
mkdir -p "$OUT/charts"

charts() {
  (cd "$STUDY" && uv run -q --python 3.12 --with-requirements requirements.txt \
     python figures/rd_charts.py --out "$OUT/charts")
  for svg in "$OUT"/charts/*.svg; do
    pdf="$OUT/$(basename "$svg" .svg).pdf"
    rsvg-convert -f pdf -o "$pdf" "$svg"
    echo "wrote $pdf"
  done
}

diagrams() {
  local pconf="${MMDC_PUPPETEER:-$HERE/diagrams/puppeteer.json}"
  for mmd in "$HERE"/diagrams/*.mmd; do
    pdf="$OUT/diagram_$(basename "$mmd" .mmd).pdf"
    mmdc -q -i "$mmd" -o "$pdf" -c "$HERE/diagrams/config.json" -p "$pconf" -f -b white
    echo "wrote $pdf"
  done
}

case "$WHAT" in
  all) charts; diagrams ;;
  charts) charts ;;
  diagrams) diagrams ;;
  *) echo "usage: make_figures.sh [all|charts|diagrams]" >&2; exit 2 ;;
esac
