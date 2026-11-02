# gnomAD frequency checks

Made by `gnomad_af.py checks` (PROTOCOL.md section 1). Full counts: MyVariant.info `gnomad_genome` (gnomAD v3.1.2 genomes). Subset rows: `gnomad.genomes.v3.1.2.hgdp_tgp.chrN.vcf.bgz` at the positions of the 100 people's variants.

## 1. Adj counts against the subset file's own AC and AN

The subset file's INFO AC and AN count the samples that are `high_quality` and not `related` in the sample metadata (3,380 of 4,151; 2,510 of them are 1000 Genomes samples, which is not the same set as the 2,500 released ones, so the released samples' counts have no INFO field to compare with directly). Recounting those samples with the adj filter used here checks the filter and the haploid rule. Rows that were multiallelic before splitting are left out.

| Rows | Checked | AC and AN both equal (adj as used here, haploid DP 10) | Same, with haploid DP 5 |
|---|---|---|---|
| autosomes | 30,842 | 30,842 of 30,842 (100.00%) | 30,842 of 30,842 (100.00%) |
| chrX outside the PARs | 1,592 | 1,592 of 1,592 (100.00%) | 4 of 1,592 (0.25%) |
| chrX PARs | 11 | 11 of 11 (100.00%) | 11 of 11 (100.00%) |
| all | 32,445 | 32,445 of 32,445 (100.00%) | 30,857 of 32,445 (95.11%) |

On chrX outside the PARs the subset file's own counts hold haploid calls to DP 10, not the protocol's 5, so the frequencies here use DP 10 for haploid calls (a deviation). The rules differ only for a haploid call with GQ 20 or more and DP 5 to 9, which changes a family member's AN by 1 allele.

Rows that came from a multiallelic row: 0 of 32,445. After splitting, a 1/2 heterozygote's allele balance uses the reference and that one alt allele only.

## 2. Full counts: the subset file's gnomad_AC and gnomad_AN against MyVariant.info

- Rows with gnomad_ fields and a MyVariant.info record: 25,535.
- Overall AC and AN equal: 25,535 of 25,535 (100.00%).
- Overall and all five groups' AC and AN equal: 25,535 of 25,535 (100.00%).
- Subset rows without gnomad_ fields (sites not in the full release): 107.

## 3. Released 1000 Genomes adj counts within gnomAD's group counts

Rows where the 2,500 released samples' count exceeds the full gnomAD count in a group: none (of 25,535 rows).

## 4. Group sizes for the all-1000-Genomes fallback

| Group | gnomAD v3.1.2 genomes (N) | Largest AN / 2 on autosomes in gnomad_full.tsv | Released 1000 Genomes |
|---|---|---|---|
| afr | 20,744 | 20,743 | 632 |
| amr | 7,647 | 7,647 | 340 |
| eas | 2,604 | 2,604 | 501 |
| nfe | 34,029 | 34,029 | 423 |
| sas | 2,419 | 2,419 | 507 |

## 5. Matching, absent rows and clamping

- main: 17,584 person variants; in the subset file 17,559 of 17,584 (99.86%); not found in gnomAD (all frequencies 0) 37; in the subset file but not found by MyVariant.info 14; negative AC clamped at 0 (group or overall counts): family left out 0, all 1000 Genomes left out 0.
- rep: 16,028 person variants; in the subset file 16,006 of 16,028 (99.86%); not found in gnomAD (all frequencies 0) 19; in the subset file but not found by MyVariant.info 0; negative AC clamped at 0 (group or overall counts): family left out 0, all 1000 Genomes left out 0.

## 6. Hand-worked examples

Twenty person variants: 3 on chrX outside the PARs in XY people, 2 in XX people, 1 in a PAR (17 person variants lie in a PAR), 2 with no subset row (fallback AN), the rest autosomal with a subset row. For each: the full counts from MyVariant.info and from the gnomAD v3.1.2 sites VCF (what the browser shows), each released family member's adj call (alt alleles / alleles), and the subtraction in the person's group and overall.

Full counts agree with the sites VCF on 20 of 20 (overall and five groups).

**chrX:40073696 G>A**, NA19258 (male, gnomAD group afr, released false)
- Full: AC 894 / AN 113514; afr 99/31286, amr 65/10888, eas 0/3617, nfe 623/53434, sas 1/2834; grpmax 0.0117. Sites VCF agrees.
- Family: NA19256 (afr, male) 0/1; NA19257 (afr, female) 1/2.
- Family left out: overall (894 - 1) / (113514 - 3) = 0.007867; afr (99 - 1) / (31286 - 3) = 0.003133; grpmax_fam 0.01166. All 1000 Genomes left out: af 0.008044, grpmax 0.01167.

**chrX:24208272 A>ATGC**, NA20867 (male, gnomAD group sas, released true)
- Full: AC 28 / AN 113333; afr 0/31261, amr 0/10762, eas 0/3612, nfe 7/53449, sas 15/2825; grpmax 0.00531. Sites VCF agrees.
- Family: NA20867 (sas, male) 1/1.
- Family left out: overall (28 - 1) / (113333 - 1) = 0.0002382; sas (15 - 1) / (2825 - 1) = 0.004958; grpmax_fam 0.004958. All 1000 Genomes left out: af 0.000219, grpmax 0.005276.

**chrX:124100892 A>C**, NA19098 (male, gnomAD group afr, released true)
- Full: AC 655 / AN 108709; afr 618/30035, amr 25/10073, eas 0/3548, nfe 2/52273, sas 1/2586; grpmax 0.0206. Sites VCF agrees.
- Family: NA19098 (afr, male) 1/1.
- Family left out: overall (655 - 1) / (108709 - 1) = 0.006016; afr (618 - 1) / (30035 - 1) = 0.02054; grpmax_fam 0.02054. All 1000 Genomes left out: af 0.005933, grpmax 0.02017.

**chrX:148997737 G>C**, HG02010 (female, gnomAD group afr, released true)
- Full: AC 355 / AN 111925; afr 336/30747, amr 14/10572, eas 0/3578, nfe 0/53208, sas 0/2692; grpmax 0.0109. Sites VCF agrees.
- Family: HG02010 (afr, female) 1/2.
- Family left out: overall (355 - 1) / (111925 - 2) = 0.003163; afr (336 - 1) / (30747 - 2) = 0.0109; grpmax_fam 0.0109. All 1000 Genomes left out: af 0.003142, grpmax 0.01081.

**chrX:129590704 C>T**, NA20768 (female, gnomAD group nfe, released true)
- Full: AC 1959 / AN 111451; afr 100/30633, amr 207/10519, eas 1/3530, nfe 1367/53063, sas 10/2632; grpmax 0.0258. Sites VCF agrees.
- Family: NA20768 (nfe, female) 1/2.
- Family left out: overall (1959 - 1) / (111451 - 2) = 0.01757; nfe (1367 - 1) / (53063 - 2) = 0.02574; grpmax_fam 0.02574. All 1000 Genomes left out: af 0.01791, grpmax 0.02573.

**chrX:2220713 G>C**, HG01508 (male, gnomAD group nfe, released false), PAR
- Full: AC 1814 / AN 152166; afr 106/41442, amr 116/15274, eas 3/5194, nfe 1221/68034, sas 133/4830; grpmax 0.0275. Sites VCF agrees.
- Family: HG01506 (nfe, male) 0/2; HG01507 (nfe, female) 1/2.
- Family left out: overall (1814 - 1) / (152166 - 4) = 0.01191; nfe (1221 - 1) / (68034 - 4) = 0.01793; grpmax_fam 0.02754. All 1000 Genomes left out: af 0.01199, grpmax 0.02725.

**chr11:125958569 G>A**, NA19098 (male, gnomAD group afr, released true)
- Full: AC 6 / AN 116114; afr 0/34970, amr 1/11360, eas 0/3818, nfe 2/50518, sas 0/3318; grpmax 8.8e-05. Sites VCF agrees.
- Family: NA19098 (afr, male) no row, 0/2.
- Family left out: overall (6 - 0) / (116114 - 2) = 5.167e-05; afr (0 - 0) / (34970 - 2) = 0; grpmax_fam 8.803e-05. All 1000 Genomes left out: af 5.343e-05, grpmax 9.212e-05.

**chr11:125958569 G>A**, HG04194 (male, gnomAD group sas, released true)
- Full: AC 6 / AN 116114; afr 0/34970, amr 1/11360, eas 0/3818, nfe 2/50518, sas 0/3318; grpmax 8.8e-05. Sites VCF agrees.
- Family: HG04194 (sas, male) no row, 0/2.
- Family left out: overall (6 - 0) / (116114 - 2) = 5.167e-05; sas (0 - 0) / (3318 - 2) = 0; grpmax_fam 8.803e-05. All 1000 Genomes left out: af 5.343e-05, grpmax 9.212e-05.

**chr9:107012057 T>C**, HG00380 (female, gnomAD group nfe, released true)
- Full: AC 757 / AN 151428; afr 51/41260, amr 102/15152, eas 0/5166, nfe 505/67940, sas 18/4788; grpmax 0.00743. Sites VCF agrees.
- Family: HG00380 (nfe, female) 1/2.
- Family left out: overall (757 - 1) / (151428 - 2) = 0.004993; nfe (505 - 1) / (67940 - 2) = 0.007419; grpmax_fam 0.007419. All 1000 Genomes left out: af 0.00504, grpmax 0.007393.

**chr2:27234145 G>A**, NA19900 (male, gnomAD group afr, released true)
- Full: AC 509 / AN 152182; afr 497/41432, amr 8/15282, eas 0/5202, nfe 2/68028, sas 0/4834; grpmax 0.012. Sites VCF agrees.
- Family: NA19900 (afr, male) 1/2.
- Family left out: overall (509 - 1) / (152182 - 2) = 0.003338; afr (497 - 1) / (41432 - 2) = 0.01197; grpmax_fam 0.01197. All 1000 Genomes left out: af 0.003384, grpmax 0.0121.

**chr7:40127519 C>A**, HG03197 (male, gnomAD group afr, released false)
- Full: AC 5 / AN 151850; afr 4/41304, amr 0/15242, eas 0/5166, nfe 1/67964, sas 0/4818; grpmax 9.68e-05. Sites VCF agrees.
- Family: HG03195 (afr, female) 0/2; HG03196 (afr, male) 1/2.
- Family left out: overall (5 - 1) / (151850 - 4) = 2.634e-05; afr (4 - 1) / (41304 - 4) = 7.264e-05; grpmax_fam 7.264e-05. All 1000 Genomes left out: af 2.724e-05, grpmax 7.493e-05.

**chr8:23019157 T>TGGCCCTG**, NA19781 (female, gnomAD group amr, released false)
- Full: AC 90 / AN 152206; afr 79/41454, amr 8/15286, eas 0/5188, nfe 0/68032, sas 0/4834; grpmax 0.00191. Sites VCF agrees.
- Family: NA19779 (amr, female) 1/2; NA19780 (amr, male) 0/2.
- Family left out: overall (90 - 1) / (152206 - 4) = 0.0005847; amr (8 - 1) / (15286 - 4) = 0.0004581; grpmax_fam 0.001906. All 1000 Genomes left out: af 0.0005978, grpmax 0.001941.

**chr15:52321396 G>C**, NA19735 (male, gnomAD group amr, released true)
- Full: AC 780 / AN 152150; afr 54/41436, amr 355/15274, eas 0/5186, nfe 306/68030, sas 26/4828; grpmax 0.0232. Sites VCF agrees.
- Family: NA19735 (amr, male) 1/2.
- Family left out: overall (780 - 1) / (152150 - 2) = 0.00512; amr (355 - 1) / (15274 - 2) = 0.02318; grpmax_fam 0.02318. All 1000 Genomes left out: af 0.00511, grpmax 0.02295.

**chr4:46244524 C>T**, NA12485 (male, gnomAD group nfe, released false)
- Full: AC 2320 / AN 151392; afr 154/41350, amr 266/15116, eas 1/5116, nfe 1443/67594, sas 63/4824; grpmax 0.0213. Sites VCF agrees.
- Family: NA12413 (nfe, male) 1/2; NA12414 (nfe, female) 0/2.
- Family left out: overall (2320 - 1) / (151392 - 4) = 0.01532; nfe (1443 - 1) / (67594 - 4) = 0.02133; grpmax_fam 0.02133. All 1000 Genomes left out: af 0.01561, grpmax 0.02124.

**chr2:209976120 G>A**, NA19055 (male, gnomAD group eas, released true)
- Full: AC 493 / AN 152154; afr 462/41428, amr 15/15278, eas 4/5190, nfe 5/68030, sas 0/4822; grpmax 0.0112. Sites VCF agrees.
- Family: NA19055 (eas, male) 1/2.
- Family left out: overall (493 - 1) / (152154 - 2) = 0.003234; eas (4 - 1) / (5190 - 2) = 0.0005783; grpmax_fam 0.01115. All 1000 Genomes left out: af 0.003241, grpmax 0.01118.

**chr17:42565760 C>T**, HG03872 (male, gnomAD group sas, released true)
- Full: AC 794 / AN 152242; afr 66/41452, amr 61/15292, eas 1/5204, nfe 518/68042, sas 119/4836; grpmax 0.0246. Sites VCF agrees.
- Family: HG03872 (sas, male) 1/2.
- Family left out: overall (794 - 1) / (152242 - 2) = 0.005209; sas (119 - 1) / (4836 - 2) = 0.02441; grpmax_fam 0.02441. All 1000 Genomes left out: af 0.005073, grpmax 0.0225.

**chr12:13562044 C>T**, NA19454 (male, gnomAD group afr, released true)
- Full: AC 111 / AN 152196; afr 108/41462, amr 2/15286, eas 0/5172, nfe 1/68032, sas 0/4832; grpmax 0.0026. Sites VCF agrees.
- Family: NA19454 (afr, male) 1/2.
- Family left out: overall (111 - 1) / (152196 - 2) = 0.0007228; afr (108 - 1) / (41462 - 2) = 0.002581; grpmax_fam 0.002581. All 1000 Genomes left out: af 0.0007065, grpmax 0.002513.

**chr9:4586790 C>A**, NA19137 (female, gnomAD group afr, released true)
- Full: AC 214 / AN 152126; afr 193/41430, amr 19/15270, eas 0/5182, nfe 1/68024, sas 0/4820; grpmax 0.00466. Sites VCF agrees.
- Family: NA19137 (afr, female) 1/2.
- Family left out: overall (214 - 1) / (152126 - 2) = 0.0014; afr (193 - 1) / (41430 - 2) = 0.004635; grpmax_fam 0.004635. All 1000 Genomes left out: af 0.00138, grpmax 0.004556.

**chr1:156868532 C>A**, HG03131 (male, gnomAD group afr, released false)
- Full: AC 61 / AN 152196; afr 58/41448, amr 2/15290, eas 0/5192, nfe 0/68022, sas 0/4824; grpmax 0.0014. Sites VCF agrees.
- Family: HG03129 (afr, female) 1/2; HG03130 (afr, male) 0/2.
- Family left out: overall (61 - 1) / (152196 - 4) = 0.0003942; afr (58 - 1) / (41448 - 4) = 0.001375; grpmax_fam 0.001375. All 1000 Genomes left out: af 0.0004076, grpmax 0.001418.

**chr5:151267815 C>CT**, HG01389 (male, gnomAD group amr, released true)
- Full: AC 739 / AN 148152; afr 627/40462, amr 50/14828, eas 2/5102, nfe 47/66750, sas 1/4724; grpmax 0.0155. Sites VCF agrees.
- Family: HG01389 (amr, male) 1/2.
- Family left out: overall (739 - 1) / (148152 - 2) = 0.004981; amr (50 - 1) / (14828 - 2) = 0.003305; grpmax_fam 0.0155. All 1000 Genomes left out: af 0.00491, grpmax 0.01515.

