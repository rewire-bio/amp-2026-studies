#!/usr/bin/env bash
# Genotypes of the 50 replication people in the panel regions (PLAN.md phase 1, PROTOCOL.md section 5).
# Runs the first run's fetch_kgp.sh for each chromosome, 8 at a time, retrying a chromosome up to 6 times because
# the 1000 Genomes FTP drops connections. Output: data/rep/kgp/chrN.bcf.
#
#   cd rare-disease-v2 && bash fetch_rep_kgp.sh
set -euo pipefail
cd "$(dirname "$0")"
HERE=$(pwd)
mkdir -p data/rep/kgp logs
one() {
  for i in 1 2 3 4 5 6; do
    out=$(cd ../rare-disease/data && ../fetch_kgp.sh "$1" "$HERE/data/rep/samples50.txt" "$HERE/data/rep/kgp" 2>/dev/null) \
      && [ "$out" = "$1 done" ] && { echo "$1 done"; return 0; }
    echo "$1 attempt $i failed"; sleep 10
  done
  return 1
}
export -f one
export HERE
printf 'chr%s\n' {1..22} X | xargs -P 8 -I{} bash -c 'one {}' | tee logs/fetch_kgp_rep.log
