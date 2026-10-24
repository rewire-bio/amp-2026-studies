#!/usr/bin/env bash
# Fetch 1000 Genomes high-coverage genotypes for the panel regions on one chromosome, 50 samples, sites present in them.
# Run from rare-disease/data:  printf "chr%s\n" {1..22} X | xargs -P 8 -n 1 ../fetch_kgp.sh
# Optional sample file and output directory (default samples50.txt, kgp), e.g. for the replication backgrounds:
#   printf "chr%s\n" {1..22} X | xargs -P 8 -I{} ../fetch_kgp.sh {} replication/samples50.txt replication/kgp
set -o pipefail
c=$1
samples=${2:-samples50.txt}
out=${3:-kgp}
mkdir -p "$out"
awk -v c="$c" '$1==c' regions.tsv > "$out/regions_$c.tsv"
[ -s "$out/regions_$c.tsv" ] || exit 0
bcftools view -R "$out/regions_$c.tsv" -S "$samples" --force-samples "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.$c.filtered.SNV_INDEL_SV_phased_panel$([ "$c" = chrX ] && echo .v2).vcf.gz" -Ou 2>/dev/null \
  | bcftools view -c 1 -Ob -o "$out/$c.bcf" && echo "$c done"
