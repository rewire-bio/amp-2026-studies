#!/usr/bin/env bash
# Make the paper's figures: charts (PDF in figures/, PNG in ../figures/ for the companion) and Mermaid diagrams (PDF in
# figures/; the companion draws the same diagrams from Mermaid blocks).
#
#   bash gpu-drift/paper/make_figures.sh            # charts and diagrams
#   bash gpu-drift/paper/make_figures.sh charts     # charts only
#   bash gpu-drift/paper/make_figures.sh diagrams   # diagrams only
#
# Charts need a checkout of rewire-bio/evo2-cloud-operations (commit f8ae11c9) next to this repository.
# Diagrams use mmdc with diagrams/config.json; diagrams/puppeteer.json points mmdc at Google Chrome (set
# MMDC_PUPPETEER to another puppeteer config to override).
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
STUDY="$(dirname "$HERE")"
WHAT="${1:-all}"
mkdir -p "$HERE/figures" "$STUDY/figures"

charts() {
  uv run -q --with matplotlib --with pandas python "$HERE/make_figures.py"
}

diagrams() {
  local pconf="${MMDC_PUPPETEER:-$HERE/diagrams/puppeteer.json}"
  for mmd in "$HERE"/diagrams/*.mmd; do
    local name="diagram_$(basename "$mmd" .mmd)"
    mmdc -q -i "$mmd" -o "$HERE/figures/$name.pdf" -c "$HERE/diagrams/config.json" -p "$pconf" -f -b white
    echo "wrote $name"
  done
}

case "$WHAT" in
  all) charts; diagrams ;;
  charts) charts ;;
  diagrams) diagrams ;;
  *) echo "usage: make_figures.sh [all|charts|diagrams]" >&2; exit 2 ;;
esac
