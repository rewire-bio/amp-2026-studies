#!/usr/bin/env bash
# Make the paper's figures in figures/ and the companion's PNGs in ../figures/: the two charts from
# ../evee_charts.py, converted from SVG, and the Mermaid diagrams in diagrams/*.mmd. Then rerun stats.py, which reads
# the chart SVGs back as a check.
#
#   bash make_figures.sh            # charts, stats and diagrams
#   bash make_figures.sh charts     # charts and stats only
#   bash make_figures.sh diagrams   # diagrams only
#
# Needs the EVEE figure data: goodfire-ai/evee-manuscript cloned next to this repository (EVEE=... to override).
# evee_charts.py writes its SVGs to ../out/ (not committed); rsvg-convert makes the PDF and PNG copies.
# Diagrams: mmdc with diagrams/config.json; diagrams/puppeteer.json points mmdc at Google Chrome (MMDC_PUPPETEER to
# override).
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
STUDY="$(dirname "$HERE")"
REPO="$(dirname "$STUDY")"
EVEE="${EVEE:-$(dirname "$REPO")/evee-manuscript/artifacts}"
WHAT="${1:-all}"
mkdir -p "$HERE/figures" "$STUDY/figures"

charts() {
  (cd "$REPO" && uv run -q --python 3.11 --with matplotlib==3.11.2 --with pandas==3.1.0 --with pyarrow==26.0.0 \
     python evee/evee_charts.py --evee "$EVEE")
  for name in evee_confounding evee_results; do
    rsvg-convert -f pdf -o "$HERE/figures/$name.pdf" "$STUDY/out/$name.svg"
    rsvg-convert -f png -w 2468 -b '#FDFCFA' -o "$STUDY/figures/$name.png" "$STUDY/out/$name.svg"
    echo "wrote figures/$name.pdf and ../figures/$name.png"
  done
  (cd "$REPO" && uv run -q --python 3.11 --with pandas==3.1.0 --with pyarrow==26.0.0 --with numpy==2.4.6 \
     python evee/paper/stats.py --evee "$EVEE")
}

diagrams() {
  local pconf="${MMDC_PUPPETEER:-$HERE/diagrams/puppeteer.json}"
  for mmd in "$HERE"/diagrams/*.mmd; do
    local name; name="diagram_$(basename "$mmd" .mmd)"
    mmdc -q -i "$mmd" -o "$HERE/figures/$name.pdf" -c "$HERE/diagrams/config.json" -p "$pconf" -f -b white
    mmdc -q -i "$mmd" -o "$HERE/figures/$name.png" -c "$HERE/diagrams/config.json" -p "$pconf" -b white -s 2
    echo "wrote figures/$name.pdf and .png"
  done
}

case "$WHAT" in
  all) charts; diagrams ;;
  charts) charts ;;
  diagrams) diagrams ;;
  *) echo "usage: make_figures.sh [all|charts|diagrams]" >&2; exit 2 ;;
esac
