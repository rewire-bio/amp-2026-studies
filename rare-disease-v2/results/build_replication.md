# build_cases.py --study replication

500 cases, 50 people (rep set), 80870 case rows, 8143 distinct variants. Seed 20261109.

Cases per person: 50 people with 10.

Background variants removed as ClinVar P/LP (any review status): 110. Background variants at the diagnostic's position: 0; overlapping the diagnostic's reference bases (VCF POS to the last base touched): 0.

Background per case, SNVs and indels counted apart. Both are kept for the inheritance filter; only the study's type (SNVs) is ranked. The inheritance filter is applied at analysis time and is not counted here.

| Type | Build filter (grpmax below 1%) | Main filter (grpmax below 0.1%) | Total rows |
|---|---|---|---|
| SNVs | median 127 (55 to 275) | median 33.5 (4 to 73) | 70330 |
| Indels | median 20 (8 to 41) | median 4 (0 to 13) | 10040 |
| Both | median 145.5 (71 to 310) | median 36 (5 to 81) | 80370 |

## Quotas

Excluded 200 diagnostic positions of rare-disease-v2/data/snv/cases.tsv.
Excluded 200 diagnostic positions of rare-disease/data/cases.tsv.
Excluded 494 diagnostic positions of rare-disease/data/replication/cases.tsv.
Excluded positions in all: 852.
Every quota filled.

Available: positions the group could still supply under the gene cap when it was drawn (after exclusions and positions taken by earlier groups).

Largest gene: the most diagnostics from one gene in the group (the cap is 10 per group).

| Group | Quota | Available | Drawn | Genes | Largest gene | mono | biallelic | xl_mono | xl_bi | unknown | het | hom | hemi |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| missense | 150 | 4068 | 150 | 110 | 6 | 74 | 67 | 6 | 2 | 1 | 78 | 68 | 4 |
| canonical_splice | 100 | 2398 | 100 | 76 | 8 | 31 | 63 | 0 | 5 | 1 | 32 | 67 | 1 |
| nonsense | 75 | 4529 | 75 | 60 | 5 | 24 | 44 | 1 | 6 | 0 | 24 | 48 | 3 |
| noncoding | 175 | 214 | 175 | 117 | 10 | 48 | 103 | 9 | 15 | 0 | 53 | 110 | 12 |
| all | 500 | | 500 | 256 | | 177 | 277 | 16 | 28 | 2 | 187 | 293 | 20 |

Diagnostic classes: missense 150, canonical_splice 100, intronic 84, nonsense 75, splice_region 63, synonymous 15, utr 13.
Diagnostics in more than one gene: 0.
Background zygosity: het 78470, hemi 1240, hom 660.
Background variants with no gene: 0.

## Checks

- every diagnostic grpmax below 0.1%: yes
- every background grpmax below 1%: yes
- every diagnostic's pool grpmax equals full gnomAD grpmax: yes
- no gene over 10 in a group: yes
- one diagnostic per position: yes
- every hemi row is male, chrX, outside the PARs: yes
- every background row joins to a frequency: yes
- every gene has a mode of inheritance: yes
- no background variant has REF and ALT both longer than 1 bp: yes
- every diagnostic is an SNV: yes
