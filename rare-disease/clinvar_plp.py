"""Keep ClinVar pathogenic or likely pathogenic variants with at least two review stars, one row per allele.

Input is the bcftools query output from prepare_data.sh; output adds the most severe consequence across
transcripts. Indels longer than 50 bp are dropped.

  python3 clinvar_plp.py clinvar_snv_panel.tsv clinvar_plp2.tsv
"""
import csv
import sys

STARS2 = {"criteria_provided,_multiple_submitters,_no_conflicts", "reviewed_by_expert_panel", "practice_guideline"}
PLP = {"Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"}
MAX_LEN = 51  # 50 bp change plus the VCF anchor base
# MC lists one SO term per transcript ("SO:0001583|missense_variant,SO:0001627|intron_variant"); keep the most
# severe, in VEP's order.
ORDER = ["splice_acceptor_variant", "splice_donor_variant", "nonsense", "frameshift_variant", "stop_lost",
         "initiator_codon_variant", "inframe_insertion", "inframe_deletion", "inframe_indel", "missense_variant",
         "synonymous_variant", "5_prime_UTR_variant", "3_prime_UTR_variant", "non-coding_transcript_variant",
         "intron_variant", "genic_upstream_transcript_variant", "genic_downstream_transcript_variant",
         "no_sequence_alteration"]


def most_severe(mc: str) -> str:
    terms = [t.split("|")[-1] for t in mc.split(",") if "|" in t]
    ranked = [t for t in ORDER if t in terms]
    return ranked[0] if ranked else "."


def main(src: str, dst: str) -> None:
    seen = set()
    with open(src) as h, open(dst, "w") as out:
        out.write("chrom\tpos\tclinvar_id\tref\talt\tclnsig\trevstat\tconsequence\tgene\tdisease\n")
        for r in csv.reader(h, delimiter="\t"):
            key = (r[0], r[1], r[3], r[4])
            if r[5] not in PLP or r[6] not in STARS2 or key in seen or "," in r[4]:
                continue
            if max(len(r[3]), len(r[4])) > MAX_LEN or not set(r[3] + r[4]) <= set("ACGT"):
                continue
            seen.add(key)
            out.write("\t".join(["chr" + r[0], r[1], r[2], r[3], r[4], r[5], r[6], most_severe(r[7]),
                                 r[8].split(":")[0], r[9]]) + "\n")
    print(f"{dst}: {len(seen)} P/LP variants with 2+ stars")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
