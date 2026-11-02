#!/usr/bin/env bash
# Build the paper with the shared build in tools/paper (see the notes there): paper.md to paper.pdf, or paper.tex.
#
#   bash build.sh                   # paper.pdf
#   bash build.sh --tex             # paper.tex
#   PAPER=other.md bash build.sh    # build another Markdown file (output next to it, same stem)
#
# The style, citation format, dash check and shared references.bib live in tools/paper. Entries only this paper
# cites can go in a references.bib here; the build adds it to the shared one.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
exec bash "$HERE/../../tools/paper/build.sh" "${PAPER:-$HERE}" "$@"
