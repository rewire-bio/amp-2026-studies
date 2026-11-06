# build_cases.py --study snv

200 cases, 50 people (main set), 99488 case rows, 19294 distinct variants. Seed 20261102.

Cases per person: 50 people with 4.

Background variants removed as ClinVar P/LP (any review status): 52. Background variants at the diagnostic's position: 0; overlapping the diagnostic's reference bases (VCF POS to the last base touched): 0.

Background per case, SNVs and indels counted apart. Both are kept for the inheritance filter; only the study's type (SNVs) is ranked. The inheritance filter is applied at analysis time and is not counted here.

| Type | Build filter (grpmax below 5%) | Main filter (grpmax below 0.1%) | Total rows |
|---|---|---|---|
| SNVs | median 299 (207 to 903) | median 51 (16 to 98) | 83268 |
| Indels | median 61 (42 to 167) | median 6 (1 to 14) | 16020 |
| Both | median 360.5 (250 to 1060) | median 58 (17 to 107) | 99288 |

## Quotas

Every quota filled.

Available: positions the group could still supply under the gene cap when it was drawn (after exclusions and positions taken by earlier groups).

Largest gene: the most diagnostics from one gene in the group (the cap is 10 per group).

| Group | Quota | Available | Drawn | Genes | Largest gene | mono | biallelic | xl_mono | xl_bi | unknown | het | hom | hemi |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| missense | 60 | 4132 | 60 | 50 | 4 | 22 | 26 | 4 | 8 | 0 | 24 | 29 | 7 |
| canonical_splice | 40 | 2463 | 40 | 33 | 4 | 12 | 24 | 2 | 2 | 0 | 13 | 26 | 1 |
| nonsense | 30 | 4569 | 30 | 26 | 3 | 13 | 10 | 4 | 3 | 0 | 15 | 11 | 4 |
| noncoding | 70 | 417 | 70 | 46 | 10 | 24 | 36 | 4 | 6 | 0 | 25 | 39 | 6 |
| all | 200 | | 200 | 133 | | 71 | 96 | 14 | 19 | 0 | 77 | 105 | 18 |

Diagnostic classes: missense 60, canonical_splice 40, intronic 38, nonsense 30, splice_region 23, synonymous 6, utr 3.
Diagnostics in more than one gene: 0.
Background zygosity: het 96584, hemi 1464, hom 1240.
Background variants with no gene: 0.

## Checks

- every diagnostic grpmax below 0.1%: yes
- every background grpmax below 5%: yes
- every diagnostic's pool grpmax equals full gnomAD grpmax: yes
- no gene over 10 in a group: yes
- one diagnostic per position: yes
- every hemi row is male, chrX, outside the PARs: yes
- every background row joins to a frequency: yes
- every gene has a mode of inheritance: yes
- no background variant has REF and ALT both longer than 1 bp: yes
- every diagnostic is an SNV: yes
