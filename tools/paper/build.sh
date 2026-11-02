#!/usr/bin/env bash
# Build a study's paper: <study>/paper/paper.md to paper.pdf (default) or paper.tex (--tex), next to the source.
#
#   bash tools/paper/build.sh rare-disease-v2/paper           # rare-disease-v2/paper/paper.pdf
#   bash tools/paper/build.sh rare-disease-v2/paper --tex     # paper.tex (standalone LaTeX, same options)
#   bash tools/paper/build.sh rare-disease-v2/paper/other.md  # another Markdown file (output next to it, same stem)
#
# Needs pandoc 3 and xelatex (TeX Live basic scheme is enough). Make the study's figures first.
#
# Shared files in this folder: header.tex (LaTeX preamble), nature.csl (citation style), nodash.lua (dash check) and
# references.bib (the shared bibliography). A references.bib in the paper folder is added to the shared one; if a key
# is in both, the study's entry is used and the build says so.
#
# Notes for paper.md:
# - YAML header: title, author, date, abstract. The affiliation goes in the author field as raw LaTeX (authblk is
#   not installed and the article template has no institute field). Either form works:
#     author: '`Tim Richardson \\ Rewire`{=latex}'
#     author: |
#       Tim Richardson`\\[4pt]\small`{=latex} Rewire
# - Citations as [@key], numbered in Nature style. End the file with "# References {-}" if you want a heading; the
#   list is placed there or at the end.
# - Figures: ![Caption.](figures/name.pdf){width=100%}, paths relative to the paper folder. pandoc-crossref is not
#   installed, so number references to figures and tables by hand ("Figure 1"); LaTeX numbers the captions.
# - Tables: pipe tables become booktabs longtables. Keep to plain cells; row spans would need multirow, which is
#   not installed.
# - nodash.lua turns en and em dashes into hyphens. In the reference list and citations that is silent; a dash
#   anywhere else fails the build (rewrite the sentence). ALLOW_DASH=1 lets it through as a hyphen with a warning.
set -euo pipefail
TOOLS="$(cd "$(dirname "$0")" && pwd)"

usage() { echo "usage: build.sh <study>/paper [--tex]   (or a .md file in place of the folder)" >&2; exit 2; }
[[ $# -ge 1 ]] || usage
TARGET="$1"; shift
FORMAT=pdf
case "${1:-}" in
  --tex) FORMAT=tex ;;
  "") ;;
  *) usage ;;
esac

if [[ -d "$TARGET" ]]; then
  SRC="$TARGET/paper.md"
else
  SRC="$TARGET"
fi
[[ -f "$SRC" ]] || { echo "build.sh: no such file: $SRC" >&2; exit 1; }
SRC_DIR="$(cd "$(dirname "$SRC")" && pwd)"
SRC="$SRC_DIR/$(basename "$SRC")"
STEM="$SRC_DIR/$(basename "$SRC" .md)"

# Bibliography: the shared file, then the study's own if there is one. Pandoc keeps the last entry for a repeated key.
BIBS=(--bibliography "$TOOLS/references.bib")
STUDY_BIB="$SRC_DIR/references.bib"
if [[ -f "$STUDY_BIB" ]]; then
  BIBS+=(--bibliography "$STUDY_BIB")
  keys() { sed -nE 's/^[[:space:]]*@[A-Za-z]+[[:space:]]*\{[[:space:]]*([^,[:space:]]+)[[:space:]]*,.*/\1/p' "$1" | sort -u; }
  DUP="$(comm -12 <(keys "$TOOLS/references.bib") <(keys "$STUDY_BIB") | tr '\n' ' ')"
  if [[ -n "$DUP" ]]; then
    echo "build.sh: keys in both the shared and the study references.bib, using the study's entry: $DUP" >&2
  fi
fi

ARGS=(
  "$SRC"
  --from markdown
  --standalone
  --citeproc
  "${BIBS[@]}"
  --csl "$TOOLS/nature.csl"
  --lua-filter "$TOOLS/nodash.lua"
  --number-sections
  --resource-path "$SRC_DIR"
  --include-in-header "$TOOLS/header.tex"
  --pdf-engine xelatex
  -V documentclass=article
  -V papersize=a4
  -V fontsize=11pt
  -V geometry:margin=2.5cm
  -V colorlinks=true
  -V linkcolor=black
  -V urlcolor=blue!50!black
  -V citecolor=black
  -M link-citations=true
  -M reference-section-title=References
)

cd "$SRC_DIR"   # so relative figure paths in the .tex resolve when it is compiled by hand
OUT="$STEM.$FORMAT"
LOG="$(mktemp)"
trap 'rm -f "$LOG"' EXIT
if ! pandoc "${ARGS[@]}" -o "$OUT" 2> >(tee "$LOG" >&2); then exit 1; fi
wait
if grep -q '^nodash.lua:' "$LOG" && [[ "${ALLOW_DASH:-0}" != 1 ]]; then
  echo "build.sh: en or em dash in the text (see above); rewrite it, or set ALLOW_DASH=1" >&2
  rm -f "$OUT"
  exit 1
fi
echo "wrote $OUT"
