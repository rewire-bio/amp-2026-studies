#!/usr/bin/env bash
# Make the paper's figures: the charts from charts.py (PDF in figures/, PNG in ../figures/ for the companion page)
# and the Mermaid diagrams in diagrams/*.mmd (PDF in figures/, PNG in ../figures/).
#
#   bash make_figures.sh            # charts and diagrams
#   bash make_figures.sh charts     # charts only
#   bash make_figures.sh diagrams   # diagrams only
#
# Diagrams: mmdc with diagrams/config.json (neutral theme) and -f (page fitted to the diagram). mmdc's bundled Chrome
# is not installed here, so diagrams/puppeteer.json points it at Google Chrome; set MMDC_PUPPETEER to override.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
STUDY="$(dirname "$HERE")"
WHAT="${1:-all}"
mkdir -p "$HERE/figures" "$STUDY/figures"

charts() {
  uv run -q --with matplotlib python "$HERE/charts.py"
}

diagrams() {
  local pconf="${MMDC_PUPPETEER:-$HERE/diagrams/puppeteer.json}"
  for mmd in "$HERE"/diagrams/*.mmd; do
    local name="diagram_$(basename "$mmd" .mmd)"
    mmdc -q -i "$mmd" -o "$HERE/figures/$name.pdf" -c "$HERE/diagrams/config.json" -p "$pconf" -f -b white
    mmdc -q -i "$mmd" -o "$STUDY/figures/$name.png" -c "$HERE/diagrams/config.json" -p "$pconf" -b white -s 2
    echo "wrote $name"
  done
}

case "$WHAT" in
  all) charts; diagrams ;;
  charts) charts ;;
  diagrams) diagrams ;;
  *) echo "usage: make_figures.sh [all|charts|diagrams]" >&2; exit 2 ;;
esac
