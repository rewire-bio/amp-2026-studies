#!/usr/bin/env bash
# Make the paper's figures in figures/ (PDF) and the companion's PNG copies in ../figures/.
#
#   bash make_figures.sh            # charts and diagrams
#   bash make_figures.sh charts     # charts only
#   bash make_figures.sh diagrams   # diagrams only
#
# Charts: make_charts.py reads run/ (written by record_run.sh) and writes SVG and PNG to figures/charts/; each SVG
# becomes figures/<name>.pdf with rsvg-convert. Diagrams: mmdc with diagrams/config.json (neutral theme), PDF and
# PNG. diagrams/puppeteer.json points mmdc at Google Chrome; set MMDC_PUPPETEER to override.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
OUT="$HERE/figures"
SITE="$(dirname "$HERE")/figures"
WHAT="${1:-all}"
mkdir -p "$OUT/charts" "$SITE"

charts() {
  (cd "$HERE" && uv run -q --python 3.12 --with matplotlib==3.11.2 python make_charts.py --out "$OUT/charts")
  for svg in "$OUT"/charts/*.svg; do
    name="$(basename "$svg" .svg)"
    rsvg-convert -f pdf -o "$OUT/$name.pdf" "$svg"
    cp "$OUT/charts/$name.png" "$SITE/$name.png"
    echo "wrote $OUT/$name.pdf and $SITE/$name.png"
  done
}

diagrams() {
  local pconf="${MMDC_PUPPETEER:-$HERE/diagrams/puppeteer.json}"
  for mmd in "$HERE"/diagrams/*.mmd; do
    name="diagram_$(basename "$mmd" .mmd)"
    mmdc -q -i "$mmd" -o "$OUT/$name.pdf" -c "$HERE/diagrams/config.json" -p "$pconf" -f -b white
    mmdc -q -i "$mmd" -o "$OUT/$name.png" -c "$HERE/diagrams/config.json" -p "$pconf" -b white -s 2
    echo "wrote $OUT/$name.pdf and .png"
  done
}

case "$WHAT" in
  all) charts; diagrams ;;
  charts) charts ;;
  diagrams) diagrams ;;
  *) echo "usage: make_figures.sh [all|charts|diagrams]" >&2; exit 2 ;;
esac
