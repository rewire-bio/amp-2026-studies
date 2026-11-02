# How these results compare with existing benchmarks

Checked on 2026-11-02 against Rewire's benchmark repositories (`rewire-benchmark-data`, `rewire-benchmarks`) and the
published literature. Every number below was checked by a second pass against its source or, where noted, against
Rewire's stored transcription. Numbers that could not be
checked are listed at the end and are not used in the argument.

## Short answer

No published benchmark has our design: variant-level top 1, no phenotype, one PanelApp panel (1,510 genes), a grpmax
filter below 0.1%, an inheritance filter, ClinVar 2-star P/LP diagnostics and a median list of 19 variants. So our
absolute figures (lookups combined 90.6%, rule 87.2%, AVI alone 83.6%) cannot be checked against outside numbers. Our
top 1 is higher than any published no-phenotype or real-case ranking result, and the reason is the design, not a
better method.

What can be checked is the direction of each finding. Most agree with the literature:

- Adding a newer score to a strong pipeline gives small gains elsewhere too. Cooperstein 2025 found +2.7% at top 10
  from adding AlphaMissense and SpliceAI to REVEL and MVP in Exomiser. Our AVI gain is +1.8 (-0.4 to +4.0).
- Evo 2 and GPN-MSA sit at or below CADD and AlphaMissense on SNVs, even in their own papers. This fits our finding
  that GPN adds nothing and Evo 2 adds +0.8 on SNVs.
- AVI's advantage is in non-coding variants in its own paper (AlphaGenome Atlas) and in the one independent
  calibration (Babazade 2026). Our only AVI gain was in non-coding diagnostics (+5.1).
- EVEE comes out about level with CADD and AlphaMissense in its own paper (0.977 against 0.963 and 0.972 AUROC). Ours
  agrees with ClinVar labels as often as the lookups do (87.0 against 87.0).
- Frequency and filter settings move ranking results by tens of points. Bone 2016 reported the same circularity we
  found in the first run: frequencies from the same people as the backgrounds make the frequency filter work too well.

Where we differ, the difference is explained by the metric. Our single-tool order (CADD 61.0, REVEL 28.6,
AlphaMissense 24.2 at top 1) is the reverse of missense benchmarks (Livesey 2025, ProteinGym), because REVEL and
AlphaMissense score only missense variants and our lists hold every consequence class.

The like-for-like published figures are the no-phenotype or variant-only results, and all are far lower than ours
(1% to 77% top 1). The gap comes from list length. The 100,000 Genomes pilot found a median of 1 variant (IQR 0 to 2)
in PanelApp virtual panels after filtering, against 221 per proband genome-wide. A gene panel, not the score, does
most of the work.

## Summary table

| Benchmark | Design | List size and filters | Metric | What it found | How it compares with ours |
|---|---|---|---|---|---|
| This study | 1000 Genomes backgrounds, ClinVar 2-star P/LP spiked in, 500 cases | PanelApp 285 panel, grpmax below 0.1% (person and family left out), inheritance filter; median 19 SNVs, 4 indels | Variant top 1, no phenotype | Lookups 90.6%, rule 87.2%, AVI alone 83.6%; AVI added +1.8 (-0.4 to +4.0) | Reference |
| Robinson 2014 (Exomiser PHIVE) | 28,516 HGMD SNVs spiked into 1,092 1000 Genomes exomes | Exome-wide; 1% MAF in some runs; mean genes 8,388, then 400, then 379 (AD) or 37 (AR) | Gene top 1 | Variant only: 28% AD, 77% AR. With phenotype: 66% AD, 83% AR | Same design family. Lists far longer; 1000 Genomes people in the frequency source |
| Bone 2016 (Exomiser hiPHIVE) | 10,000 simulated exomes, HGMD variants in 1000 Genomes exomes | Exome-wide, 1% MAF; ESP plus 1000 Genomes, or ESP only | Variant top 1 | Variant only: 20 to 77% with 1000 Genomes frequencies, 1 to 10% with ESP only. Full phenotype: 96 to 97% | Direct precedent for our first-run circularity |
| Smedley 2016 (Genomiser) | 453 non-coding Mendelian variants in 1,092 1000 Genomes genomes | Genome-wide, MAF above 1% removed (frequency source includes the backgrounds) | Variant top 1 | 77% full phenotype, 23% no phenotype; CADD alone 0.2% top 10 | Only spike-in benchmark for non-coding diagnostics. Not comparable in list size |
| Jacobsen 2022, Hum Mutat (100,000 Genomes) | 4,877 solved probands | Genome-wide, MAF 0.1% (dominant) or 2% (recessive), segregation, coding only | Variant top 1 | 82.6% with phenotype; variant score only 33%; 25 variants per case at variant score 0.8 or above | Variant-only 33% is the closest real-patient figure. Their list at a variant-score threshold is similar to our 19 |
| 100,000 Genomes pilot, NEJM 2021 | 490 SNV and indel diagnoses | PanelApp virtual panels plus frequency, consequence, segregation | Variant top 1 | Exomiser 77% top 1. Median 1 variant (IQR 0 to 2) in panels, 221 genome-wide | Shows list size sets the ceiling |
| Cipriani 2020 (retinal exomes) | 134 solved patients | Exome-wide, mean 474 variants after filtering | Variant top 1 | 73.9% with phenotype; 3.0% variant only. CADD 72.4, REVEL 71.6 with phenotype | Variant-only 3% on lists 25 times longer than ours |
| Talos, Nat Med 2026 (ACG, RGP) | Real exomes and genomes of ill children | Rule-based, PanelApp-driven unranked short list | Recall within list | 77% to 90% recall at 0.8 to 2.2 candidates per proband | Closest analogue to our rule (87.2%). Not a top-1 metric |
| Talos paper, Exomiser v14 (ACG trios) | 194 families (190 in Extended Data Fig. 1) | Genome-wide, default settings | Variant top k | Top 1 65%, top 5 81%, top 10 85%, any 89% | Real cases, phenotype likely used. Lower level expected |
| Jacobsen 2022 reply, Brief Bioinform | 161 100,000 Genomes singletons, Exomiser 12.1.0 | Settings changed one at a time | Gene top 1 | 34%, then 44% (MOI frequency), 62% (non-PASS removed), 72% (REVEL and MVP) | Settings move top 1 by 38 points on the same cases |
| Yuan 2022 and 2024 (DDD 305, plus in-house) | Real solved exomes, HPO | Genome-wide | Gene top 1 | Exomiser 15.1% (2022 defaults) to 61.6% (2024 Protocol D, trio, REVEL and MVP) on DDD | Same cases, fourfold range from settings |
| Cooperstein 2025 (UDN) | 231 GS and 125 ES coding probands | Genome-wide, HPO, ClinVar whitelist off | Variant top 10 | Tuning 49.7% to 85.5% (GS). AlphaMissense plus SpliceAI added to REVEL plus MVP: +2.7% | Only published added-score result. Same size as ours |
| AlphaGenome Atlas, GREGoR (developer) | 112 P/LP variants in solved genomes | Whole genome; arm with gnomAD AF below 0.001 (105 variants) | Recall at top K | AF-filtered top 10: AVI 36.2%, CADD 27.6% | Closest model-paper design. Baseline is CADD alone, lists far longer |
| CAGI6 RGP (Stenton 2024) | 14 solved test families, blinded | Genome-wide, HPO, diagnostics not already in ClinVar | Variant rank | Best 9 of 14 at rank 1; Exomiser 6 of 14 | Harder diagnostics. Not comparable |
| Kafkas 2025 | Synthetic gene sets of 5 to 100 genes, ClinVar cases, HPO | Set size 25 near our 19 | Gene Hits@1 | Exomiser 84% at 25 genes | Similar level, but gene-level, phenotype-driven, and the authors state a circularity |
| Evo 2 paper, SD1 (developer) | ClinVar Oct 2024, 2+ star, P/LP vs B/LB | No frequency filter, no case lists | AUROC | Coding SNV: AlphaMissense 0.958, CADD 0.927, Evo 2 7B 0.830 | Agrees: Evo 2 below the lookups on SNVs |
| EVEE v4 (developer) | Deconfounded ClinVar, 158,616 SNVs | Consequence-balanced, not frequency-matched | AUROC | EVEE 0.977, AlphaMissense 0.972, CADD 0.963, AlphaGenome 0.574 | Agrees: EVEE level with the lookups |
| GPN-Star, Nature 2026 (developer) | ClinVar missense; OMIM non-coding vs gnomAD common | No case lists | AUPRC | GPN-Star (M) best on OMIM non-coding (bars read from image only) | Not comparable; controls defined by frequency |
| TraitGym Mendelian | Non-coding causal variants, 9 matched controls | No case lists | AUPRC | Zero-shot CADD 0.712, GPN-MSA 0.694, Evo 2 7B 0.385 | Agrees: CADD at or above GPN-MSA and Evo 2 on non-coding |
| Babazade 2026 (independent, abstract) | ClinGen-style calibration of AVI and SpliceAI2 | ClinVar likely benign or population controls | Evidence strength | AVI no better than AlphaMissense on missense; strong in exon-adjacent non-coding classes | Agrees with our AVI pattern |
| Livesey and Marsh 2025 | 819 proteins, 83 VEPs, ClinVar vs gnomAD missense | Missense only | Rank score from per-protein AUROC | AlphaMissense 77.5, REVEL 76.0, CADD 1.7 52.1 | Reverse of our single-tool order; explained by coverage |

## 1. Rare-disease ranking benchmarks

Almost every published ranking benchmark uses phenotype (HPO terms) and ranks across the whole exome or genome. Their
headline top-1 rates with phenotype:

- Exomiser on real 100,000 Genomes cases: 77% (NEJM 2021) and 82.6% (Jacobsen 2022, 4,877 cases).
- Exomiser on retinal exomes: 73.9% (Cipriani 2020).
- Exomiser v14 on Australian Genomics trios: 65% (Talos paper).
- LIRICAL and Exomiser on 384 case reports spiked into one exome: 80.7% and 77.3% at gene rank 1 (Robinson 2020).
- AMELIE on 215 DDD singletons: 66% at gene rank 1, Exomiser 38%, with a median of 127 candidate genes (Birgmeier 2020).
- Simulated exomes with full phenotype: 96 to 97% (Bone 2016).

These are not comparable with ours because phenotype carries much of the ranking. Jacobsen 2022 split Exomiser's
top hit by score component: combined 82%, phenotype score only 55%, variant score only 33%.

The no-phenotype results are the fair comparison, and all are lower than ours:

| Source | List | Variant-only or no-phenotype result |
|---|---|---|
| Robinson 2014 | Exome-wide, 379 genes (AD) or 37 (AR) after filters | 28% AD, 77% AR top 1 |
| Bone 2016 | Exome-wide, 1% MAF | 20 to 77% top 1 with 1000 Genomes frequencies; 1 to 10% with ESP only |
| Genomiser 2016 | Genome-wide, non-coding diagnostics | 23% top 1 |
| Jacobsen 2022 | Genome-wide, 100,000 Genomes | 33% top 1 |
| Cipriani 2020 | Mean 474 variants | 3.0% top 1, 55% top 10 |
| Cooperstein 2025 | Genome-wide, filtered VCF | 40.9% top 10 |
| Ruzicka 2025 (abstract) | Median shortlist of 9, nephrology | 42% top 1 without HPO |
| AlphaGenome Atlas, GREGoR | Whole genome, AF below 0.001 | AVI 36.2%, CADD 27.6% at top 10 |
| This study | Median 19, one panel | Lookups 90.6%, AVI alone 83.6% top 1 |

List length explains most of this. The NEJM 2021 pilot found a median of 1 variant (IQR 0 to 2) left in PanelApp
virtual panels after filtering, against 221 per proband genome-wide. Our median of 19 sits between these because our
cases are singletons with no segregation filter. Ruzicka 2025 is the one exception that does not fit simply: 42% on
lists of about 9. Its diagnostics are real nephrology cases, not ClinVar 2-star variants, which are easier.

Settings also move results a great deal on identical cases. On 161 100,000 Genomes singletons, Exomiser 12.1.0 went
from 34% top 1 to 44% with inheritance-specific frequency cutoffs (0.1% dominant or homozygous recessive, 2% compound
heterozygous), to 62% after removing non-PASS (low quality) variants, and to 72% with REVEL and MVP in place of
PolyPhen-2, MutationTaster and SIFT (Jacobsen reply, 2022). On 305 DDD cases, Exomiser top 1 ranged from 15.1% (Yuan
2022, defaults) to 31.8% (Yuan 2024, Protocol A) and 61.6% (Protocol D, trio mode with REVEL and MVP). The REVEL
and MVP step adds about 10 points in the Jacobsen reply and 21.3 points in Yuan 2024 (Protocol A to B). This agrees
with our finding that a combined set of established lookups beats any single score (90.6% against CADD 61.0%).

Two results test adding newer scores. Cooperstein 2025 (UDN, Exomiser 14) found that adding AlphaMissense and
SpliceAI together to REVEL and MVP raised the share of diagnostics in the top 10 "by 2.7%" over default settings. The
paper does not say whether that is points or a relative change, and the pair is not tested separately. Cipriani 2020
found that swapping in CADD or REVEL moved top 1 by about 1.5 to 2.3 points from the default, with phenotype in use,
and the differences were not significant. Both are the same size as our AVI gain. Neither gives an interval.

The rule-based Talos is the closest analogue to our ClinGen-points rule. It returns an unranked list of 1.3
candidates per proband and recovers 90% of known diagnoses in ACG trios, 82% in ACG singletons at 1.8 candidates, and
86% to 87% in RGP. Our rule ranks the diagnostic first in 87.2% of cases. Talos measures recall in a short list, not
rank, so the numbers are not directly comparable, but rule-based filtering on real cases reaches a similar level.

The blinded CAGI6 Rare Genomes Project challenge is harder by design. Its 14 solved test families had diagnostics
that were not already P/LP in ClinVar or DM in HGMD. The best teams placed the causal variant at rank 1 in 9 of 14,
and Exomiser in 6 of 14, with phenotype and segregation available. Our ClinVar 2-star diagnostics are the opposite
case: variants with enough evidence to be classified.

## 2. What the model papers claim, and why ClinVar AUROC does not carry over

Every model paper except AlphaGenome Atlas reports classification: pathogenic against benign (or against gnomAD
common variants), scored by AUROC or AUPRC over large label sets, with no case structure, no phenotype and no frequency
filter on candidates. The claims, as checked:

- **Evo 2** (Brixi 2026, SD1, ClinVar October 2024, 2+ star). Coding SNV AUROC: AlphaMissense 0.958, CADD 0.927,
  ESM-1b 0.923, GPN-MSA 0.896, Evo 2 40B 0.841, Evo 2 7B 0.830. Non-coding SNV: SpliceAI 0.992, CADD 0.991, Evo 2 40B
  0.987, GPN-MSA 0.981, Evo 2 7B 0.978. Coding indel: Evo 2 7B 0.906, ESM-1b 0.809; CADD and AlphaMissense were not
  scored for indels. On non-coding non-SNVs only the 40B model leads (0.969); 7B (0.930) is just below PhyloP241
  (0.936). Intergenic SNVs: CADD 0.913, Evo 2 7B 0.744.
- **EVEE** (Pearce 2026, v4, deconfounded ClinVar, 158,616 SNVs). EVEE 0.977, AlphaMissense 0.972 (missense only),
  Evo 2 mean probe 0.967, CADD 1.7 0.963, GPN-MSA 0.952, Evo 2 zero-shot 0.929, AlphaGenome 0.574. Indels (73,961,
  SNV probe applied zero-shot): EVEE 0.991, CADD 1.7 0.980. The authors state that the probe "is supervised on ClinVar
  labels", which gives it a home-domain advantage on ClinVar benchmarks.
- **GPN-MSA** (Benegas, bioRxiv v2). Against gnomAD common controls (MAF above 5%) it beats CADD. With ClinVar
  benign controls, the setting closest to ours, the authors say it "performs slightly behind CADD and ESM-1b".
  Ablation ClinVar AUROC 0.969.
- **GPN-Star** (Ye 2026). GPN-Star (V) is best on ClinVar missense, level with ESM-1b, and GPN-Star (M) is
  best on OMIM and HGMD non-coding variants against 9 matched gnomAD common controls. The figure values are not
  printed; from the image, OMIM non-coding AUPRC is about 0.76 for GPN-Star (M) and about 0.70 for CADD. The authors
  note that AlphaMissense and PrimateAI-3D use population frequencies, which also define benign variants.
- **AlphaGenome Atlas, AVI** (Cheng 2026, developer). ClinVar SNV AUPRC by consequence, AVI against the best other
  method: intronic 0.76 against 0.44, synonymous 0.57 against 0.35, 3' UTR 0.50 against 0.18, protein-altering 0.90
  against 0.86, 5' UTR 0.26 against 0.27. Retraining without consequence features lowers ClinVar AUPRC by 13.71%;
  without AlphaGenome features by 2.56%; without AlphaMissense by 2.51%. On GREGoR (the only ranking analysis among the
  model papers) AVI beats CADD 1.7 at top 50, 29.5% against 12.5% on whole genomes, and 74.3% against 61.0% after an
  AF filter below 0.001.
- **AlphaMissense** (Cheng 2023). ClinVar test auROC 0.940 against EVE 0.911; DDD 0.809 against PrimateAI 0.797, a
  difference that was not significant.
- **CADD 1.7** (Schubach 2024). ClinVar 2+ star AUROC 0.982. The authors recommend ranking by CADD rather than using a
  cutoff.

The direction of these claims agrees with ours. Evo 2 and GPN-MSA trail CADD and AlphaMissense on SNVs. AVI's margin
is large in intronic, synonymous and 3' UTR variants and small on protein-altering variants. AlphaGenome alone is
weak on coding-dominated ClinVar.

The size of the claims does not carry over, for three reasons.

1. **The benign class is partly defined by frequency, and our filter removes that signal.** GPN-MSA, GPN-Star and
   TraitGym use gnomAD common variants as negatives. AlphaMissense defines benign as frequent in human and primate
   populations. AVI's training label is gnomAD v4.1 grpmax FAF: 0.001 or above is benign, below 0.001 is impactful.
   ClinVar benign labels also often rest on frequency (BA1, BS1). ProteinGym gives a within-tool example: BayesDel
   scores 0.972 AUC with its allele frequency feature and 0.925 without. Our grpmax filter below 0.1% removes every
   common variant before any score sees the list, so a score's ability to recognise common variants earns nothing.
2. **Consequence class carries much of the ClinVar signal, and our lookups already include it.** EVEE's aggregate
   SNV AUROC is 0.997, which its authors call confounded by consequence; within a consequence class it ranges from
   0.900 (nonsense) to 0.984 (intronic). AVI loses most when consequence features are removed (13.71%). Our combined
   lookups have consequence-class terms, so a new score must add information within a class.
3. **There is little room at the top.** The lookups rank the diagnostic first in 90.6% of cases and in every
   canonical splice and nonsense case. An added score can only gain in the remaining tenth, mostly non-coding
   diagnostics, which is where AVI gained.

CADD shows the gap between the two tasks clearly. Its AUROC on ClinVar 2+ star variants is 0.982, yet alone it ranks
the diagnostic first in 61.0% of our cases.

The Evo 2 indel results cannot be tested against the literature. The Evo 2 paper reports no CADD or AlphaMissense
scores for indels, so its indel lead (0.906 coding) is over conservation and protein models only. No outside source
tests Evo 2 against CADD on indels. Our indel result (+1.7, -5.9 to +9.9) and the first run's +15.8 have no external
comparator.

## 3. Independent benchmarks and calibration

Most comparator values above were reported by the developers of the model being tested (Evo 2, EVEE, GPN-MSA,
GPN-Star, AlphaGenome Atlas, SpliceAI2), and TraitGym was built by the GPN-MSA group. Independent results:

- **Babazade 2026** (bioRxiv, abstract checked) calibrated AVI and SpliceAI2 as PP3/BP4 evidence. On missense, AVI
  reached moderate to strong but did not outperform its component AlphaMissense. In non-coding classes it reached
  strong evidence in the three classes next to exons, with ClinVar controls. SpliceAI2 discriminated better than AVI
  in two classes and in a splicing assay. The number of ClinVar likely benign variants changed about 60-fold at 20
  nucleotides from the exon without a measured change in splicing, so thresholds depend on the choice of controls.
  This is the only independent clinical evaluation of AVI found, and it agrees with ours.
- **Livesey and Marsh 2025** ranked 97 VEPs on 36 deep mutational scans and 83 VEPs on ClinVar missense in 819
  proteins. Clinical rank scores: AlphaMissense 77.5 (4th), REVEL 76.0 (8th), CADD 1.7 52.1 (44th). On the mutational
  scans REVEL (66.1) and CADD (40.3) rank lower than on ClinVar, and AlphaMissense (90.7) and EVE (77.3) higher. The
  authors attribute the clinical lead of ClinVar-trained tools "almost certainly" to circularity. They class
  AlphaMissense as population-tuned.
- **Livesey and Marsh 2023** (26 mutational scans, 55 VEPs) found that mutational-scan and clinical rankings agree
  closely for unsupervised VEPs (Spearman 0.946). Supervised tools such as MetaRNN, ClinPred, BayesDel, VEST4 and
  MutPred do much better on clinical data than on mutational scans. The clinical positives were ClinVar plus HGMD.
- **ProteinGym clinical** (commit 144fe22). Supervised: ClinPred 0.981, MetaRNN 0.977, REVEL 0.928, CADD 0.905.
  Zero-shot: PoET and TranceptEVE L 0.920, EVE 0.914. AlphaMissense, SpliceAI, Evo 2 and GPN are not included.
- **Grimm 2015** defined two types of circularity: overlap of training and test variants, and genes whose variants
  are all labelled one way. Tools that seemed better than CADD and SIFT lost their lead on independent data.
- **Splicing.** On Genomics England rare-disease blood RNA-seq, AlphaGenome does not beat SpliceAI (auPRC: SpliceAI
  0.48, AlphaGenome splice-site 0.45, composite 0.36; SpliceAI2 0.59; Jaganathan 2026, developer-reported). In
  Alavi 2026, AlphaGenome's splice-junction score reached 96.2% accuracy (v2 rerun, after a coordinate error in v1
  was fixed) against SpliceAI's 94.5% (v1 table) on 183 held-out variants. Both fit a small AVI gain over a pipeline that already has
  SpliceAI.
- **Mathur 2026** is an independent Evo 2 evaluation on mitochondrial variants only (AUROC 0.896) and does not bear
  on our nuclear panel.
- No independent clinical evaluation of GPN-MSA or GPN-Star was found.

Calibration papers are inputs to our rule, not comparisons. Pejaver 2022 set the REVEL levels (PP3 supporting 0.644,
moderate 0.773, strong 0.932) on ClinVar 2019 missense variants. Bergquist 2025 added AlphaMissense (+1 at 0.792 up
to +4 at 0.990) and notes that AlphaMissense used gnomAD variants as a proxy for benign. Walker 2023 set SpliceAI
PP3 at 0.2 and BP4 at 0.1 from assays in 8 genes (BRCA1, BRCA2, MLH1, MSH2, MSH6, PMS2, NF1, POU1F1). So the rule's
missense levels are calibrated on ClinVar, which favours the rule on ClinVar diagnostics, and its splice level is
calibrated on 8 genes that are not in an intellectual disability panel.

## 4. Circularity and frequency leakage

Our first run put the background people in gnomAD. 1000 Genomes is part of gnomAD, so every background variant had
been seen there, while most inserted diagnostics had not. That made "not in gnomAD" a strong signal for the
diagnostic. The same problem has been seen before:

- **Bone 2016** stated that using 1000 Genomes frequencies with 1000 Genomes backgrounds "provides unrealistically
  effective frequency filtration". Their variant-only top 1 fell from 20 to 77% to 1 to 10% when those frequencies
  were removed. This is the closest published precedent for our first-run finding.
- **Robinson 2014** used 1000 Genomes frequencies on 1000 Genomes backgrounds.
- **Genomiser 2016** did the same and does not discuss it.
- **LIRICAL** (one Genome in a Bottle exome against gnomAD), **PhEval** (one unnamed template exome), **Cipriani
  2020** (patients in the local UCLex frequency source), **Jacobsen 2022** (local 100,000 Genomes frequencies) and the
  **GREGoR** analysis in AlphaGenome Atlas (gnomAD v4.1) do not say whether the background or patient was removed
  from the frequency source.
- **Xrare** kept separate background sets for features, training and evaluation, and used internal genomes, so it
  avoided the problem by design.

Our fix leaves each case's person and family out of gnomAD. It does not remove presence as a signal. In the first
run's 200 cases, 76% of the spiked-in diagnostics are absent from gnomAD (AF 0). Of the background variants, 0% are
absent at full frequency, 2% once the family is left out, and 4% once all of 1000 Genomes is left out, counted before
the frequency filter (`results/gnomad_overlap.md`). In this run, 386 of the 500 diagnostics (77.2%) are absent. A
variant absent from gnomAD always passes the grpmax filter, so absent background variants make up a larger share of
the filtered lists than 2%, but they remain a small minority while diagnostics are mostly absent.

This matters most for scores trained on frequency. AVI's label threshold is grpmax FAF 0.001, the same cut as our
filter. After the filter, nearly every candidate is on the "impactful" side of AVI's label (the label uses gnomAD v4.1
FAF95 grpmax; our filter uses v3.1.2 grpmax with the family left out), so AVI cannot use frequency
classes to separate them. What is left is the difference between very rare and absent, which AVI may still have
learned. AlphaMissense is in the same position: its benign class is defined by being seen in populations. The paper's
split check did not flag the primary contrast, but the gain is entirely in cases whose diagnostic is absent from
gnomAD, and the check cannot exclude a presence effect. No published benchmark has measured this effect for a
frequency-trained score.

The first run's larger gains fit this picture: AVI +6.2 (0.0 to +13.0) and Evo 2 on indels +15.8, against +1.8 and
+1.7 now. The AlphaGenome Atlas GREGoR analysis shows the same pattern from the other side: the AVI gap over CADD at
top 50 falls from 17.0 points to 13.3 once an AF filter is applied, and CADD's recall rises from 12.5% to 61.0%.

ClinVar labels add a second circularity:

- Many 2-star P/LP diagnostics used PP3 evidence from CADD, REVEL or SpliceAI, so they favour the lookups and the rule.
- REVEL's training set overlaps ClinVar pathogenic missense, and the rule's thresholds were calibrated on ClinVar.
- EVEE's released scores come from a probe trained on ClinVar labels with an 80/20 gene split. If the released
  scores come from those models, about 80% of genes would be in-sample; the paper does not say, so this is our
  inference. Our EVEE result (87.0 against 87.0) is therefore a measure of agreement with ClinVar labels, not of
  diagnostic value.
- Evo 2, GPN-MSA and GPN-Star use no human variant labels in training, so they carry no label leak. AlphaMissense
  uses ClinVar for early stopping and calibration. CADD includes PolyPhen.

## 5. Rewire's own benchmark data

`rewire-benchmark-data` holds no record that matches our design. The relevant stores are:

- `use-case-coverage-20260930/clinical`: Talos and the Exomiser v14 ACG comparison.
- `use-case-coverage-rare-ranking-20261009`: Yuan 2022 (SM Table 3), Yuan 2024 (Table 1) and Kafkas 2025.
- The 2026-10-08 dossier `docs/reviews/use-cases/rare-disease-candidate-ranking-2026-10-08.md`: Cooperstein 2025 and
  the Genomics England Exomiser performance page (v13.1.0 top 3 91%, 95% CI 90 to 93%, on 1,869 diagnostic variants;
  v12.0.0 80%). These two are in the dossier, not stored as records.
- Model and VEP stores: Evo 2 SD1, EVEE v4, GPN-MSA, TraitGym, Livesey 2025, ProteinGym, Pejaver 2022, Bergquist
  2025, AlphaGenome Nature Supplementary Table 4, AlphaGenome Atlas (release de2862b), Alavi 2026 v1 and v2,
  SpliceAI2 and Feng 2025.
- GPN-Star has only a predictor record (status "discovered") and no benchmark values. Our GPN-Star results have no
  stored external comparison.

`rewire-benchmarks` holds one relevant benchmark: MFASS v2 (Chong 2018 minigene splicing labels, 8,324 held-out
variants). AUROC: k-mer baseline 0.778, SpliceAI 0.806, Pangolin 0.876, DNABERT-2 frozen pair with a logistic head
0.550. SpliceAI minus baseline was +0.028 (-0.004 to +0.064). It uses the same reporting style as this study (paired
intervals by group resampling, a reading stated in advance). The one genomic language model protocol tested there
fell below a simple baseline, in line with our GPN and Evo 2 results. Feng 2025 (stored) agrees in direction:
frozen DNA model embeddings separate pathogenic from common SNPs weakly (best AUROC 0.732, NT-v2).

Stored values that are relevant but less direct:

- AlphaGenome Supplementary Table 4 (ClinVar 2025-03-23, held-out chromosomes, AUPRC): deep intronic and synonymous
  0.66 against Pangolin 0.641; splice-site region 0.57 against 0.554. Small margins in non-coding, as in our study.
- Exomiser top 1 across the stored records ranges from 15.1% (Yuan 2022 defaults) to 68.4% (Yuan 2024, Protocol D,
  in-house trios). The range is set by configuration and cohort.

## 6. What this means and what this study adds

What the comparison supports:

- Our absolute top 1 is a property of the design. It should not be compared with Exomiser's 65 to 83% on real cases
  or used to claim that lookups beat phenotype-driven tools. The fair statement is that within one panel, with rare
  variants and an inheritance filter, the lists are short enough that established scores rank a ClinVar 2-star
  diagnostic first about 9 times in 10.
- Our main result, that AVI adds little to the established lookups, agrees with every source that tests something
  similar: Cooperstein (+2.7% at top 10), Cipriani (CADD or REVEL change top 1 by about 2 points), Babazade (AVI adds
  nothing over AlphaMissense on missense) and the AlphaGenome Atlas's own AF-filtered arm.
- Our non-coding AVI gain matches where AVI leads in its own paper and in Babazade 2026.
- Our GPN and Evo 2 SNV results match the order in the Evo 2, EVEE and TraitGym papers.
- Our first-run inflation matches the circularity Bone 2016 described.

What this study adds that we did not find elsewhere:

- A paired test of adding a genomic foundation model (AVI, GPN-MSA, GPN-Star, Evo 2) to a combined lookup baseline,
  inside case lists, with intervals. Published model papers compare single scores on classification; Cooperstein
  tested AlphaMissense and SpliceAI together without intervals.
- A measured frequency fix: leaving each person and family out of gnomAD, and a check on how much presence signal
  remains (76% of diagnostics absent against 2% of background).
- A rule-based ClinGen-points baseline (87.2%) compared directly with a fitted model (90.6%), on the same cases.
- No-phenotype results by consequence group, which show where an added score can still gain.

What it cannot show: performance on patients' genomes, on longer lists, on recessive diagnoses made of two variants,
or on diagnostics not already classified in ClinVar. CAGI6 and the 100,000 Genomes results show those cases are harder.

## Could not confirm

These were not checked against the source and are not used above as evidence:

- Ruzicka 2025: REVEL 12.9% and AlphaMissense and CADD 0% top 1 (search snippet only; full text returned 403).
- Babazade 2026: full text, thresholds and counts (abstract only).
- GPN-Star Figure 2 bar values (read from the image; no source data).
- Mathur 2026 and the Evo 2 SD1 non-coding cells cited in one claim were checked only against Rewire's transcription.
- Talos: the 194 against 190 denominator is unresolved.
- Cooperstein 2025: top-1 values (Figure 3) and frequency settings (Additional file 4).
- PhEval top-k values, AI-MARRVEL and Xrare comparator figures (figures only).
- CAGI6 Table 2 was truncated in the fetch.
- Moon vendor claims (94% top 3 on 150 exomes) from press coverage.
- Pejaver 2022: the text and Table 2 disagree on whether CADD 20 is BP4 moderate or supporting.

## Sources

Rare-disease ranking

- Robinson PN et al. Genome Res 2014;24:340. PMC3912424, Results "Overall, our method ranked", Figure 4.
- Smedley D et al. Nat Protoc 2015;10:2004. PMC5467691, hiPHIVE subsection.
- Bone WP et al. Genet Med 2016;18(6):608. doi 10.1038/gim.2015.137, PMC4916229, Table 1, Figure 2a.
- Smedley D et al. (Genomiser). Am J Hum Genet 2016;99(3):595. PMC5011059, Performance Evaluation, Figures 3 and 4.
- Cipriani V et al. Genes 2020;11(4):460. PMC7230372, Tables 3 and 4.
- Robinson PN et al. (LIRICAL). Am J Hum Genet 2020;107(3):403. PMC7477017, Table 1, Figures 2 and 4.
- Birgmeier J et al. (AMELIE). Sci Transl Med 2020;12:eaau9113. PMC9366928.
- Smedley D et al. (100,000 Genomes pilot). N Engl J Med 2021;385:1868. PMC7613219, Table 3, Figure 2D.
- Jacobsen JOB et al. Hum Mutat 2022;43(8):1071. doi 10.1002/humu.24380, PMC9288531, Figures 1 to 3.
- Jacobsen JOB et al. (reply to Yuan). Brief Bioinform 2022;23(5):bbac188. PMC9487604, Figure 1.
- Yuan X et al. Brief Bioinform 2022;23(2):bbac019. PMC8921623, SM Table 3; stored in
  `rewire-benchmark-data/data/omics/use-case-coverage-rare-ranking-20261009/claims.csv`.
- Yuan X et al. Sci Rep 2024;14:2845. PMC10838329, Table 1; same store.
- Kafkas S et al. Sci Rep 2025. PMC12041562, Table 3; same store.
- Li Q et al. (Xrare). Genet Med 2019;21:2126. PMC6752318.
- Mao D et al. (AI-MARRVEL). NEJM AI 2024;1(5). PMC11221788.
- Cooperstein IB et al. Genome Med 2025;17:127. doi 10.1186/s13073-025-01546-1, PMC12539062, Figures 2 and 3.
- Bridges Y et al. (PhEval). BMC Bioinformatics 2025;26:87. PMC11929307.
- Stenton SL et al. (CAGI6 RGP). Hum Genomics 2024;18:44. PMC11057178, Table 2.
- Welland MJ et al. (Talos). Nat Med 2026. doi 10.1038/s41591-026-04477-5, Table 1 and "Comparison with other tools";
  stored in `rewire-benchmark-data/data/omics/use-case-coverage-20260930/clinical/claims.csv`.
- O'Brien TD et al. Genet Med 2022;24(1):192 (abstract).
- Ruzicka J et al. medRxiv 2025, doi 10.1101/2025.02.04.25321641 (abstract).
- Genomics England. Exomiser performance page,
  https://pipeline-rd-help.genomicsengland.co.uk/Mira/variant-prioritisation-approaches/exomiser/exomiser-performance/

Model papers

- Brixi G et al. (Evo 2). Nature 2026;652:1349. doi 10.1038/s41586-026-10176-5, Supplementary Data 1 Sheet1, SI A.3.13;
  stored in `rewire-benchmark-data/data/omics/gfm-guides-20261014-evo2-brixi2026/claims.csv`.
- Pearce et al. (EVEE). bioRxiv v4, doi 10.64898/2026.04.10.717844, Figure 1c and 1d, lines 405 to 407; stored in
  `rewire-benchmark-data/data/omics/gfm-guide-sources-20261014/batch.jsonl`.
- Benegas G et al. (GPN-MSA). bioRxiv v2, doi 10.1101/2023.10.10.561776, Supplementary Table 4, Supplementary Figure 4.
- Ye et al. (GPN-Star). Nature 2026, doi 10.1038/s41586-026-11005-5, Figure 2, SI Table 3, SI Figure 1.
- Cheng J, Taylor KR et al. (AlphaGenome Atlas). medRxiv v1, doi 10.64898/2026.09.16.26363192, Results p.5 and p.7,
  Figure 3A, figure S9A, Table S2, Methods p.20 to 22 and p.29 to 30.
- AlphaGenome. Nature 2026, Supplementary Table 4 cells L2:M4 and L8:M8; stored in
  `rewire-benchmark-data/data/evidence/results`.
- Cheng J et al. (AlphaMissense). Science 2023;381:eadg7492.
- Schubach M et al. (CADD 1.7). Nucleic Acids Res 2024;52:D1143. PMC10767851.
- Benegas G et al. (TraitGym). Hugging Face revision 1fde19555fe8, `mendelian_traits_matched_9`; stored in
  `rewire-benchmark-data/data/omics/lit-genomics-20261011-traitgym/claims.csv`.
- Jaganathan K et al. (SpliceAI2). bioRxiv 2026, Figure 4g; stored in
  `rewire-benchmark-data/data/omics/use-case-coverage-rna-splicing-20261013-spliceai2/claims.csv`.
- Alavi 2026. bioRxiv doi 10.64898/2026.10.02.756289 v2, Table IV; stored in
  `rewire-benchmark-data/data/omics/alavi-v2-20261101/`.

Independent benchmarks and calibration

- Babazade et al. bioRxiv 2026, doi 10.64898/2026.10.21.760715 (abstract).
- Livesey BJ, Marsh JA. Genome Biol 2025;26:104. PMC12016141, Tables S3 and S8; stored in
  `rewire-benchmark-data/data/omics/vep-literature-20261014/claims.csv`.
- Livesey BJ, Marsh JA. Mol Syst Biol 2023;19:e11474. PMC10407742, Figure 5, Figure EV6.
- Notin P et al. (ProteinGym). Clinical summaries at commit 144fe22; stored in
  `rewire-benchmark-data/data/omics/vep-literature-20261014/artifacts/`.
- Grimm DG et al. Hum Mutat 2015;36(5):513. PMC4409520, Tables 3 and 4.
- Mathur, Sachidanandam. bioRxiv 2026, doi 10.64898/2026.03.10.710786, Table 2 (store transcription).
- Feng et al. 2025, Table 5; stored as `amp-feng-20261007-*` in `rewire-benchmark-data/data/evidence/`.
- Pejaver V et al. Am J Hum Genet 2022;109:2163. PMC9748256, Tables 1 to 3.
- Bergquist T et al. Genet Med 2025;27:101402, Table 1.
- Walker LC et al. Am J Hum Genet 2023;110(7):1046. PMC10357475, Table 1.

Rewire

- `rewire-benchmark-data/docs/reviews/use-cases/rare-disease-candidate-ranking-2026-10-08.md`, sections 1 to 3.
- `rewire-benchmark-data/data/omics/vep-use-cases-20261014/ledger.jsonl` (GPN-Star, no benchmark values).
- `rewire-benchmarks/benchmarks/mfass/README.md`, "MFASS-v2 results".
- This study: `results/gnomad_overlap.md`, `paper/paper.md` sections 2.9, 2.10, 3.1 and Limits.
