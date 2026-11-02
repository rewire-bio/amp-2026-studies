# Draft for the Discussion: Comparison with other benchmarks

No published benchmark ranks variants without phenotype inside one gene panel after a grpmax and an inheritance
filter, so our top 1 cannot be set against outside figures. With phenotype, Exomiser ranks the diagnostic first in 65
to 83% of real cases genome-wide in Australian Genomics and 100,000 Genomes cohorts [@welland2026talos; @smedley2021pilot; @jacobsen2022phenotype], and variant-only
scoring reaches 1 to 77% on exome or genome lists [@robinson2014improved; @bone2016computational;
@cipriani2020improved; @jacobsen2022phenotype]. Our higher rate reflects short lists, not a better method: after
filtering, the 100,000 Genomes pilot left a median of 1 variant in PanelApp panels against 221 genome-wide
[@smedley2021pilot]. The direction of our results agrees with the literature. Adding AlphaMissense and SpliceAI to
an Exomiser pipeline raised top 10 by a reported 2.7% [@cooperstein2025optimized], close to our AVI gain of +1.8. AVI leads
mainly on non-coding variants in its own paper [@cheng2026alphagenome] and adds nothing over AlphaMissense on
missense in an independent calibration [@babazade2026calibration]. Evo 2 trails CADD and AlphaMissense on coding SNVs
[@brixi2026genome], and GPN-MSA trails CADD with ClinVar benign controls [@benegas2025dna]. High AUROC on ClinVar does not carry over to our lists,
because many benchmarks define the benign class by population frequency, which our filter has already removed, and
because consequence class, which our lookups include, carries much of the signal [@pearce2026interpretable]. Bone et al. saw
the circularity of our first run: frequencies from the same 1000 Genomes people as the backgrounds gave
"unrealistically effective frequency filtration" [@bone2016computational].

Word count: about 220.

## Citation keys

Already in `tools/paper/references.bib`:

- `cheng2026alphagenome` (AlphaGenome Atlas, medRxiv 10.64898/2026.09.16.26363192)
- `brixi2026genome` (Evo 2, Nature)
- `benegas2025dna` (GPN-MSA, Nat Biotechnol). The "slightly behind CADD and ESM-1b" statement was checked in bioRxiv v2
  (Supplementary Fig. 4); confirm it is in the Nat Biotechnol version before citing this key for it.
- `pearce2026interpretable` (EVEE paper, bioRxiv v4; `pearce2026evee` is the dataset entry)

New entries needed:

- `welland2026talos`: Welland MJ et al. Automated reanalysis of genomic data for rare disease diagnostics at scale.
  Nat Med 2026. doi 10.1038/s41591-026-04477-5
- `smedley2021pilot`: Smedley D et al. 100,000 Genomes Pilot on Rare-Disease Diagnosis in Health Care. N Engl J Med
  2021;385:1868-1880. doi 10.1056/NEJMoa2035790
- `jacobsen2022phenotype`: Jacobsen JOB et al. Phenotype-driven approaches to enhance variant prioritization and
  diagnosis of rare disease. Hum Mutat 2022;43(8):1071-1081. doi 10.1002/humu.24380
- `robinson2014improved`: Robinson PN et al. Improved exome prioritization of disease genes through cross-species
  phenotype comparison. Genome Res 2014;24:340-348. doi 10.1101/gr.160325.113
- `bone2016computational`: Bone WP et al. Computational evaluation of exome sequence data using human and model
  organism phenotypes improves diagnostic efficiency. Genet Med 2016;18(6):608-617. doi 10.1038/gim.2015.137
- `cipriani2020improved`: Cipriani V et al. An improved phenotype-driven tool for rare Mendelian variant
  prioritization: benchmarking Exomiser on real patient whole-exome data. Genes 2020;11(4):460.
  doi 10.3390/genes11040460
- `cooperstein2025optimized`: Cooperstein IB et al. An optimized variant prioritization process for rare disease
  diagnostics: recommendations for Exomiser and Genomiser. Genome Med 2025;17:127. doi 10.1186/s13073-025-01546-1
- `babazade2026calibration`: Babazade et al. Calibration of AlphaGenome Variant Impact and SpliceAI2 scores as
  PP3/BP4 evidence for missense, non-coding, and synonymous variants. bioRxiv 2026.
  doi 10.64898/2026.10.21.760715

The DOIs for Robinson 2014 and Smedley 2021 were not among the checked claims; confirm them before adding.
