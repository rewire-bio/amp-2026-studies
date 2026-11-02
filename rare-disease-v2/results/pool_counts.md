# Diagnostic pool counts

Made by `pool.py` (PROTOCOL.md sections 3 to 5). Counts are of eligible variants (full gnomAD v3.1.2 genome grpmax below 0.1%) and of their distinct positions, since a diagnostic position is used once.
SNVs with no CADD row: 0; with no MANE row in CADD (class from all rows): 0.

## Feasibility

- SNV missense: feasible. 8543 positions for 210 together (spare 8333); 8331 not first-run for the 500's 150 (spare 8181); under the cap the 500 can draw 4083 (spare 3933). After a simulated 200, the spare for the 500 is 3905 at worst and 3917 at the median. Feasible even if the 200 took only non-first-run positions. The cap bites: 192 genes with more than 10 positions (PAH 338, SCN1A 271, NF1 222).
- SNV canonical_splice: feasible. 3660 positions for 140 together (spare 3520); 3521 not first-run for the 500's 100 (spare 3421); under the cap the 500 can draw 2412 (spare 2312). After a simulated 200, the spare for the 500 is 2285 at worst and 2295 at the median. Feasible even if the 200 took only non-first-run positions. The cap bites: 89 genes with more than 10 positions (NF1 189, DMD 111, TSC2 77).
- SNV nonsense: feasible. 8488 positions for 105 together (spare 8383); 8382 not first-run for the 500's 75 (spare 8307); under the cap the 500 can draw 4538 (spare 4463). After a simulated 200, the spare for the 500 is 4448 at worst and 4454 at the median. Feasible even if the 200 took only non-first-run positions. The cap bites: 205 genes with more than 10 positions (NF1 419, DMD 337, CEP290 170).
- SNV noncoding: feasible. 515 positions for 245 together (spare 270); 276 not first-run for the 500's 175 (spare 101); under the cap the 500 can draw 246 (spare 71). After a simulated 200, the spare for the 500 is 29 at worst and 40 at the median. Feasible even if the 200 took only non-first-run positions. The cap bites: 5 genes with more than 10 positions (NF1 75, PAH 31, TSC2 16).
- Indel frameshift: feasible. 7561 positions, 3792 under the cap, for 40 (spare 3752).
- Indel inframe: feasible. 187 positions, 187 under the cap, for 30 (spare 157).
- Indel canonical_splice: feasible. 337 positions, 323 under the cap, for 30 (spare 293).
- Indel noncoding: feasible. 36 positions, 36 under the cap, for 20 (spare 16).

## Notes

- CADD class against ClinVar's consequence on the whole SNV pool: 3 of 9190 missense_variant, 9 of 8645 nonsense, 4 of 1859 splice_acceptor_variant, 4 of 2427 splice_donor_variant differ (table below). PROTOCOL.md section 4 reports 0 on its smaller check.
- Inheritance mode unknown (PanelApp gives no mode): GNAS, 60 eligible variants. Section 3 defines no diagnostic genotype for these.

## SNV pool

22968 ClinVar SNVs; 22878 eligible (full gnomAD grpmax below 0.1%). First-run diagnostic positions: 694.

The 200 may use any eligible position. The 500 may not use first-run diagnostic positions or this run's 200. Cap: at most 10 diagnostics per gene per group in each study. Capacity under the cap is the sum over genes of min(positions, 10).

| Group | Eligible variants | Positions | Not first-run | Genes | Genes over 10 | Cap capacity, 200 | Cap capacity, 500 (not first-run) | Quota 200 | Quota 500 | Together |
|---|---|---|---|---|---|---|---|---|---|---|
| missense | 9405 | 8543 | 8331 | 898 | 192 | 4132 | 4083 | 60 | 150 | 210 |
| canonical_splice | 4270 | 3660 | 3521 | 609 | 89 | 2463 | 2412 | 40 | 100 | 140 |
| nonsense | 8613 | 8488 | 8382 | 858 | 205 | 4569 | 4538 | 30 | 75 | 105 |
| noncoding | 550 | 515 | 276 | 212 | 5 | 417 | 246 | 70 | 175 | 245 |

Simulated draws (200 seeds): the 200 drawn at random under the cap, then the positions left for the 500 (not first-run, not the 200) under the cap, minus the 500's quota. Negative means a shortfall.

| Group | Spare for the 500: min | median | max | Largest genes (positions) |
|---|---|---|---|---|
| missense | 3905 | 3917 | 3927 | PAH 338, SCN1A 271, NF1 222 |
| canonical_splice | 2285 | 2295 | 2303 | NF1 189, DMD 111, TSC2 77 |
| nonsense | 4448 | 4454 | 4461 | NF1 419, DMD 337, CEP290 170 |
| noncoding | 29 | 40 | 51 | NF1 75, PAH 31, TSC2 16 |

## Indel pool

8490 ClinVar indels up to 50 bp, not counting 159 variants with REF and ALT both longer than 1 bp, which are left out (by class: frameshift 78, inframe 60, canonical_splice 21); 8458 eligible.

| Group | Eligible variants | Positions | Genes | Genes over 10 | Cap capacity | Quota 120 | Spare |
|---|---|---|---|---|---|---|---|
| frameshift | 7887 | 7561 | 781 | 175 | 3792 | 40 | 3752 |
| inframe | 191 | 187 | 122 | 0 | 187 | 30 | 157 |
| canonical_splice | 343 | 337 | 183 | 2 | 323 | 30 | 293 |
| noncoding | 37 | 36 | 28 | 0 | 36 | 20 | 16 |

## SNV non-coding group by CADD class

| CADD class | Eligible variants | Positions |
|---|---|---|
| splice_region | 179 | 176 |
| synonymous | 33 | 33 |
| utr | 22 | 22 |
| intronic | 316 | 284 |

## SNV groups by inheritance mode (eligible positions)

| Group | mono | biallelic | xl_mono | xl_bi | unknown | none |
|---|---|---|---|---|---|---|
| missense | 3836 | 3929 | 358 | 390 | 30 | 0 |
| canonical_splice | 998 | 2373 | 87 | 197 | 5 | 0 |
| nonsense | 2972 | 4653 | 293 | 558 | 12 | 0 |
| noncoding | 186 | 267 | 23 | 38 | 1 | 0 |

## Indel groups by inheritance mode (eligible positions)

| Group | mono | biallelic | xl_mono | xl_bi | unknown | none |
|---|---|---|---|---|---|---|
| frameshift | 2311 | 4766 | 205 | 273 | 6 | 0 |
| inframe | 76 | 96 | 1 | 12 | 2 | 0 |
| canonical_splice | 89 | 220 | 8 | 19 | 1 | 0 |
| noncoding | 16 | 16 | 3 | 1 | 0 | 0 |

## ClinVar consequence against the SNV class (all SNVs)

| ClinVar consequence | nonsense | missense | canonical_splice | splice_region | synonymous | utr | intronic | other |
|---|---|---|---|---|---|---|---|---|
| missense_variant | 0 | 9187 | 0 | 0 | 0 | 0 | 3 | 0 |
| nonsense | 8636 | 0 | 0 | 1 | 0 | 7 | 1 | 0 |
| splice_donor_variant | 0 | 1 | 2423 | 0 | 0 | 1 | 2 | 0 |
| splice_acceptor_variant | 2 | 1 | 1855 | 0 | 0 | 1 | 0 | 0 |
| intron_variant | 0 | 0 | 0 | 97 | 0 | 0 | 310 | 1 |
| initiator_codon_variant | 0 | 257 | 0 | 0 | 0 | 2 | 0 | 0 |
| synonymous_variant | 0 | 0 | 0 | 81 | 33 | 0 | 1 | 0 |
| stop_lost | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 21 |
| non-coding_transcript_variant | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 13 |
| 5_prime_UTR_variant | 0 | 0 | 0 | 0 | 0 | 7 | 1 | 5 |
| 3_prime_UTR_variant | 0 | 0 | 0 | 0 | 0 | 9 | 0 | 0 |
| . | 3 | 3 | 0 | 0 | 0 | 0 | 0 | 0 |

## ClinVar consequence against the indel class (all indels)

| ClinVar consequence | frameshift | inframe | canonical_splice | utr | intronic | other |
|---|---|---|---|---|---|---|
| frameshift_variant | 7388 | 0 | 10 | 0 | 2 | 0 |
| nonsense | 449 | 14 | 2 | 0 | 0 | 0 |
| splice_donor_variant | 46 | 0 | 198 | 1 | 2 | 0 |
| splice_acceptor_variant | 17 | 1 | 113 | 0 | 2 | 0 |
| inframe_deletion | 0 | 121 | 0 | 0 | 0 | 0 |
| inframe_insertion | 0 | 32 | 0 | 0 | 0 | 0 |
| . | 10 | 4 | 16 | 0 | 0 | 0 |
| intron_variant | 1 | 1 | 4 | 0 | 21 | 0 |
| inframe_indel | 0 | 11 | 0 | 0 | 0 | 0 |
| non-coding_transcript_variant | 0 | 2 | 0 | 7 | 0 | 0 |
| initiator_codon_variant | 0 | 6 | 0 | 0 | 0 | 0 |
| 5_prime_UTR_variant | 0 | 0 | 0 | 3 | 1 | 0 |
| genic_downstream_transcript_variant | 1 | 0 | 3 | 0 | 0 | 0 |
| genic_upstream_transcript_variant | 1 | 0 | 0 | 0 | 0 | 0 |
