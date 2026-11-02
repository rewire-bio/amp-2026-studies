# Protocol: rare-disease case tests, second run

Status: revised 2 November 2026 after a design review (below), not yet frozen. The freeze (section 13) comes
before the 500 replication cases are built. After it, any change goes in the deviation log (section 14) with its
date and reason.

Why this run exists and what it asks: [`PROPOSAL.md`](PROPOSAL.md). Work plan: [`PLAN.md`](PLAN.md).

## 1. Frequencies a patient's pipeline would see

A real patient is not in gnomAD, nor usually are their relatives. Everyone else is. So for each case, the gnomAD
frequency of every background variant is recomputed **without that person and their released relatives**
("leave the family out"). This is the main frequency. Removing every 1000 Genomes sample (the literal fix in issue
[#1](https://github.com/rewire-bio/amp-2026-studies/issues/1)) is a sensitivity check (section 9): it also removes
other people of the same population, who would be in gnomAD for a real patient, and about a fifth of gnomAD's eas
and sas genomes.

Source: gnomAD v3.1.2 genomes (76,156 genomes). Here grpmax is the highest AF over afr, amr, eas, nfe and sas, the
groups gnomAD v4 uses for genome grpmax (v3.1.2 calls the related field `AF_popmax` and the leftover group `oth`).

- **Full counts.** AC and AN, overall and per group, from MyVariant.info's `gnomad_genome` field (gnomAD v3.1.2
  genomes; the first run used the same source for AF).
- **Family counts.** gnomAD's own genotype calls for the 1000 Genomes samples, from the HGDP + 1KG subset release
  (`gnomad.genomes.v3.1.2.hgdp_tgp.chrN.vcf.bgz`, per-sample GT, GQ, DP, AD), read by position. Only samples with
  `gnomad_release = true` in `gnomad.genomes.v3.1.2.hgdp_1kg_subset_sample_meta.tsv.bgz` are counted, each in its
  `gnomad_population_inference.pop` group. A genotype counts only if it passes gnomAD's adj filter: GQ at least 20,
  DP at least 10 (5 for haploid calls), and for heterozygotes an alt allele balance of at least 0.2. Sample IDs are
  matched after stripping a trailing A or B (NA12546B, NA12830A, NA18874A).
- **Family.** The person, plus their parents, children and full siblings in the 1000 Genomes pedigree file
  (`1kGP.3202_samples.pedigree_info.txt`). Undocumented relatedness is not handled (a stated limit).
- **Formula.** For group g: AF = (AC_g minus AC_fam,g) / (AN_g minus AN_fam,g), with AC clamped at 0. AN_fam is
  the adj-called alleles of the family members in that group (2 per diploid call, 1 per haploid call). On chrX
  outside the pseudoautosomal regions (PAR1 chrX:10,001-2,781,479; PAR2 chrX:155,701,383-156,030,895), haploid
  calls count 1 allele. Where the HGDP + 1KG file has no row, no family member carries the allele, so AC_fam is 0
  and AN_fam is taken as 2 (or 1) per released family member in the group.
- **Absent variants.** A variant not in gnomAD, or with no allele left, has AF 0 and is kept, as a pipeline keeps a
  variant it cannot find. A gnomAD site that is not PASS is treated as present with its counts.
- **Diagnostic variants** are not carried by the person, so their frequency is the full gnomAD frequency.

Checks before use: on 20 variants (five on chrX outside the PARs, one in a PAR), the full counts match the gnomAD
browser and the subtraction is worked by hand. The adj counts of all 2,500 released 1000 Genomes samples must
reproduce the subset file's own counts on at least 99% of biallelic rows checked.

The first run's people include 16 (of 50) who are not in the gnomAD release; for them only released relatives are
subtracted. 702 of the 3,202 1000 Genomes high-coverage samples are not in gnomAD, mostly relatives removed by
gnomAD's relatedness pruning.

## 2. Frequency filter

- **Main filter (dominant mode).** Keep a variant only if its grpmax AF (section 1) is below 0.1%. This copies the
  raw per-group maximum AF filter that Exomiser and the Genomics England tiering use, on purpose. The filtering
  allele frequency (faf95) and a disease-specific maximum credible AF (Whiffin et al. 2017) were considered and not
  used.
- **Looser background filter.** The same cases with grpmax AF below 1%. This tests how much a denser background
  lowers top 1. It is not a test of recessive diagnoses. Its locked models are refitted on the 200 under the same
  1% filter.
- Cases are built once with background variants below 1% grpmax AF, and every case variant's frequency columns are
  stored in `case_variants.tsv`. The main filter is applied at analysis time by one function shared by every score.
  All diagnostics pass the main filter, so no case loses its diagnostic under either filter and case counts are
  fixed before any score is looked at.
- Diagnostic candidates must have full gnomAD grpmax AF below 0.1%.

## 3. Mode of inheritance

In the first run, 94 of 200 and 231 of 500 diagnostics were single heterozygous variants in genes where PanelApp
gives only biallelic inheritance: a carrier state, not a diagnosis. This run makes every case a diagnosis and
filters the background as a pipeline would.

Each green gene's mode of inheritance comes from PanelApp 285 (`mode_of_inheritance`), matched to the MANE gene of
the region the variant falls in (by Ensembl gene ID, then symbol). Each person's sex comes from the 1000 Genomes
pedigree file.

- **Diagnostic genotype.** Heterozygous in MONOALLELIC and BOTH genes, and in X-linked genes marked "monoallelic
  mutations in females may cause disease". Hemizygous for XY people in any X-linked gene. Homozygous in BIALLELIC
  genes, and for XX people in X-linked genes marked "biallelic mutations in females", as in a consanguineous family.
  Scores do not use zygosity; it only decides which background variants a pipeline would keep.
- **Background filter (main analysis).** A background variant is kept if its gene's mode allows its genotype:
  any genotype in MONOALLELIC, BOTH and Unknown genes; any genotype for XY people in X-linked genes; in BIALLELIC
  genes (and X-linked "biallelic in females" genes for XX people), only if homozygous or if the person has at least
  two variants in that gene after the frequency filter (a possible compound heterozygote; phase is ignored). The
  diagnostic counts toward that total. Variants outside any green gene's region are kept.
- **Sensitivity.** The same analysis without the background inheritance filter (section 9).

## 4. One consequence annotation

The same annotation chooses the diagnostic variants, sets the model's classes and defines the reported groups.

- **SNVs:** CADD v1.7's most severe consequence over its rows for MANE v1.5 Select transcripts (CADD `FeatureID`,
  version stripped), using the severity order in `lookup_scores.py`. A variant with no MANE row takes the most
  severe consequence over all rows, and the count is logged. CADD is looked up for the whole ClinVar candidate pool
  before cases are built.
- **Indels:** `annotate_indels.py` (MANE Select), run on the ClinVar indel pool and on the case variants.

| Class | SNV (CADD) | Indel (`annotate_indels.py`) |
|---|---|---|
| nonsense | `STOP_GAINED` | |
| frameshift | | frameshift |
| inframe | | inframe |
| missense | `NON_SYNONYMOUS` | |
| canonical splice | `CANONICAL_SPLICE` | splice |
| splice region | `SPLICE_SITE` | |
| synonymous | `SYNONYMOUS` | |
| UTR | `5PRIME_UTR`, `3PRIME_UTR` | utr |
| intronic | `INTRONIC` | intronic |
| other | `STOP_LOST` and anything else | other |

For this pool, CADD and ClinVar's MC field agree on every canonical splice, nonsense and missense variant checked
(0 of 367, 398 and 415 differ). The change of annotation moves variants within the non-coding group (about a third
become splice region) and moves ClinVar initiator codon variants to missense.

**Model classes.** For the fitted models only, synonymous, UTR and other are merged into one class, and any class
with fewer than 10 diagnostics in the 200 is merged into it too. The merge list is fixed from the 200 and recorded at
the freeze. Reported groups keep the full list.

## 5. Cases

| Study | Cases | People | Diagnostic quotas |
|---|---|---|---|
| SNV | 200, 4 per person | the first run's 50 (`../rare-disease/data/samples50.txt`) | missense 60, canonical splice 40, nonsense 30, non-coding 70 (synonymous, UTR, intronic, splice region) |
| Indel (up to 50 bp) | 120, ranked among the person's background indels | the same 50 | frameshift 40, inframe 30, canonical splice 30, UTR or intronic 20 |
| Replication | 500, 10 per person | 50 new people (below) | missense 150, canonical splice 100, nonsense 75, non-coding 175 |

- **Replication people.** 50 people drawn from the 2,504 unrelated 1000 Genomes samples, leaving out both of the
  first run's sets of 50 and their families (seed 20261110). They share no one with the first run, so the 500 test
  new people and new diagnostic draws.
- **Replication diagnostics.** No position used as a diagnostic in this run's 200 or in the first run's 200 or 500.
  If a class has too few positions, its quota is reduced to what is available, and the shortfall is logged before
  the build.
- **Gene cap.** No gene supplies more than 10 diagnostics in one class of one study (NF1 is 17% of the non-coding
  pool).
- Inputs fixed to the first run's: ClinVar fileDate 2026-10-11, P/LP with 2+ review stars; PanelApp 285 v11.45 green
  genes; MANE v1.5 Select exons plus 20 bp; 1000 Genomes high-coverage phased panel (20220422).
- Background: the person's variants with 1000 Genomes AF below 1% in the regions, then grpmax AF below 1%
  (section 2). Background variants that are ClinVar P/LP are removed, so each case has one right answer.
- Seeds: SNV 20261102, indel 20261103, replication 20261109.
- The pool holds enough: after the 0.1% filter, non-coding SNVs have about 540 positions for the 245 needed by the
  200 and 500 together, and UTR or intronic indels 30 to 37 for 20. These counts are rechecked with this run's
  annotation and inheritance rules before building.

## 6. Scores

| Score | Role | Covers |
|---|---|---|
| CADD v1.7 PHRED | lookup | all SNVs; indels seen in gnomAD (not used for indels, below) |
| SpliceAI max delta | lookup | SNVs from CADD's annotation; indels run ourselves (`run_spliceai.sh`) |
| REVEL v1.3 | lookup | missense |
| AlphaMissense | lookup | missense |
| AlphaGenome AVI | primary added score | SNVs |
| GPN-MSA, GPN-Star V, GPN-Star M | secondary added scores | SNVs |
| EVEE | descriptive, its own analysis (section 10) | ClinVar variants only |

- A missing score ranks last, as a failed lookup would.
- CADD and precomputed SpliceAI are unavailable for indels in every model and in the rule-based baseline: CADD's
  indel table has a row only for indels seen in gnomAD, which would give the answer away.

## 7. Baselines

**Lookups combined.** Logistic regression on the four lookups, model class (section 4) and score-by-class terms, as
in the first run.

**Rule-based.** Each variant gets ACMG points from the rule or tool for its class, and variants are sorted by
points. This is what a lab applying ClinGen guidance would do with the same tools.

| Variants | Points |
|---|---|
| Nonsense SNVs, frameshift indels, canonical splice SNVs and indels | +8 (PVS1 tier) |
| Missense | the higher of the REVEL level and the SpliceAI level |
| Inframe indels | the higher of +2 (PM4) and the SpliceAI level |
| Everything else | the SpliceAI level |

- **REVEL level** (Bergquist et al. 2025, Genet Med 27:101402, Table 1): +4 at 0.932 or more, +3 at 0.879, +2 at
  0.773, +1 at 0.644; -1 at 0.290 or less, -2 at 0.183, -3 at 0.052, -4 at 0.016; otherwise 0. These match Pejaver
  et al. 2022 (AJHG 109:2163, Table 2), which also gives a BP4 very strong level at 0.003 or less, merged into -4
  here. Where REVEL is missing, AlphaMissense from the same table: +4 at 0.990, +3 at 0.972, +2 at 0.906, +1 at
  0.792; -1 at 0.169 or less, -2 at 0.099, -3 at 0.070; otherwise 0. ClinGen asks for one tool chosen in advance, so
  REVEL is the tool and AlphaMissense only fills gaps.
- **SpliceAI level** (Walker et al. 2023, AJHG 110:1046, Table 1): +1 at a max delta of 0.2 or more, -1 at 0.1 or
  less, 0 between. Walker calibrates these to moderate but recommends applying them at supporting weight, as here.
- CADD has no calibration outside missense (Pejaver's thresholds are missense only), so it sets no level.
- **Ties.** Within a points level, SNVs are sorted by the score that set the level (REVEL, or AlphaMissense where
  REVEL is missing, for missense; SpliceAI otherwise; none for nonsense), then by CADD PHRED. Indels are sorted by
  our SpliceAI score. A variant with no score for its class ranks below every scored variant at the same level.
  Remaining ties count against the score, as everywhere, and the share of diagnostics in a tie is reported.
- **Rule + AVI.** The same points, with ties within a level broken by AVI first.

The rule does not check whether loss of function is the disease mechanism for a gene, or whether a nonsense or
frameshift variant escapes nonsense-mediated decay.

## 8. Models, test sets and outcomes

- Cross-validated: 5-fold, grouped by person, on the 200 SNV cases and on the 120 indel cases.
- Locked: every model fitted once on all 200 SNV cases (main filter, inheritance filter) and applied unchanged to
  the 500.
- Combinations: lookups combined; lookups plus one added score at a time. Adding a score never changes another
  combination.

**Primary.** Top-1 of lookups plus AVI minus top-1 of lookups combined, both locked, on the 500, with the main
filter and the inheritance filter, with a paired 95% interval. The two models differ only by AVI.

Reading, fixed in advance (smallest gain of interest 5 points):
- "AVI adds a gain that matters": the interval lies above 0 and the estimate is at least 5 points, and the key
  secondary points the same way.
- "No gain that matters": the interval lies inside -5 to +5.
- Otherwise "inconclusive".

**Key secondary.** Rule + AVI minus the rule-based baseline, on the same cases: does AVI help the stronger
per-type baseline, or only make up for the combined model's weakness (issue #2, problem 3)?

**Secondary, adjusted.** GPN-MSA, GPN-Star V and GPN-Star M added to lookups combined: a family of three, top 1,
locked, on the 500, each with a Bonferroni-adjusted 98.33% paired interval. A secondary score is said to add
something only if that interval excludes zero.

**Descriptive** (unadjusted 95% intervals, labelled descriptive): top 5 and top 10; each diagnostic group with at
least 20 cases, and the non-coding group by CADD class; lookups combined minus the rule-based baseline, and lookups
plus AVI minus the rule-based baseline; the 200 cross-validated and the 500 cross-validated; the indel cases;
section 9's checks; each single tool's top 1 by class.

## 9. Checks

- **AVI training chromosomes.** AVI was trained on gnomAD frequency classes on chromosomes 1, 4, 7, 8, 10, 13 and
  15, from a gnomAD that includes the 1000 Genomes people; the new frequencies do not change that. Cases are split
  by whether the diagnostic is on a training chromosome, and each case is ranked only among its background variants
  on chromosomes of the same split, with the locked models unchanged. In each split: the primary contrast and the
  key secondary, with paired intervals, and the difference between splits. If the gain is larger on training
  chromosomes and the interval of the difference excludes zero, the primary result is reported as possibly
  inflated by training overlap.
- **Looser background filter** (grpmax below 1%, section 2).
- **No background inheritance filter** (section 3).
- **All of 1000 Genomes left out** (issue #1's literal fix): every released 1000 Genomes sample is subtracted instead
  of the family, using the adj counts of all 2,500 samples; where the subset file has no row, AN is reduced in
  proportion to each group's share of 1000 Genomes samples.
- **Ties** in the rule-based baseline (section 7).

## 10. EVEE: descriptive only

EVEE's released table (Zenodo 19701997, CC BY 4.0) has one row per ClinVar variant; its labels date from about
2026-03-21 (latest `last_evaluated` 2026-03-20; the release date is not stated). Measured on the first run's cases:

- 199 of 200 SNV diagnostics and 71 of 120 indel diagnostics have a row, against 45% of background SNVs.
- With the 0.1% filter the median case keeps 12 background variants with a row; with the 1% filter, 57.
- 197 of the 199 SNV diagnostics were already P/LP in EVEE's labels, and most covered background variants are B/LB,
  EVEE's training negatives. Whether the released scores are held-out predictions is not documented (the
  manuscript's five-fold scheme covers its benchmark only); their AUROC on ClinVar SNVs (0.997) equals the held-out
  benchmark, which rules out heavy overfitting but not training-label exposure.

So EVEE cannot be tested fairly here. It is reported once, descriptively: on the 200 and 500 SNV cases with the
looser 1% filter, each list restricted to variants with an EVEE row (the "ClinVar-covered" lists), EVEE and every
single score ranked on those same lists. The result is labelled as agreement with ClinVar labels on in-sample
variants, not as diagnostic performance, and is not compared with the main analysis. EVEE rows are joined at
`pos` equal to VCF POS minus 1 (EVEE's positions are 0-based), checked against ClinVar variation IDs.

## 11. Intervals

- 1,000 bootstrap resamples of people (each with all their cases).
- Each stratum's draws come from `np.random.default_rng(int.from_bytes(sha256(f"{study}/{analysis}/{stratum}").digest()[:8], "big"))`,
  with people sorted by ID before drawing. The draws are made before any score is evaluated and reused by every
  score and every paired difference in that stratum, so adding a score never moves another interval.
- Every gain has a paired interval. Charts are dot-and-interval.
- `case_benchmark.py` runs with the package versions pinned in `requirements.txt`.

## 12. Adding Evo 2 later

- Case files are frozen with SHA-256 checksums in `results/case_checksums.txt`: the 200 SNV and 120 indel files at
  tag `rd-v2-protocol`, the 500 replication files at tag `rd-v2-cases`, before they are scored.
- `evo2/make_windows.py` takes a required `--fasta` (hg38 from `../evo2-size-by-use-case/data/raw/hg38.fa.gz`) and an
  `--indels-only` option.
- Evo 2 on the indel cases is an exploratory analysis, run only if Tim asks for it (about 2,100 indels, roughly
  USD 3 on one H100). With 120 cases a top-1 gain has an interval of about plus or minus 7 points, and there is no
  held-out indel set, so it cannot confirm the first run's indel gain.
- Check after adding it: every row of every result file that does not involve Evo 2 is byte-for-byte the same as a
  rerun without it in the same environment.

## 13. Freeze

Before the 500 replication cases are built: this protocol, the analysis code, `requirements.txt` and the 200 SNV and
120 indel case files with their checksums are committed and tagged `rd-v2-protocol`. The freeze guards against
tuning on this run's 500; it does not undo what is already known from the first run (for example AVI's +0.2 under
grpmax, which helped make grpmax the main filter).

Decisions closed at the freeze and recorded below: the frequency method (section 1, done), the rule-based baseline
(section 7, done), the model class merges (section 4, from the 200), any quota reductions (section 5).

## 14. Deviation log

| Date | Section | Change | Reason |
|---|---|---|---|
| 2026-11-02 | 1 | Frequencies leave out the person and their family, not every 1000 Genomes sample; all of 1000 Genomes left out becomes a check | Neither gnomAD release carries subtractable 1000 Genomes counts; per-sample genotypes are available, and leaving out the family is what a real patient's situation looks like (design review) |
| 2026-11-02 | 3 | Inheritance mode added | About half of the first run's diagnostics were carrier states (design review) |
| 2026-11-02 | 5 | 50 new replication people | The first run's replication people and their results were already seen (design review) |
| 2026-11-02 | 7, 8 | Rule-based baseline uses ACMG points; primary contrast is lookups plus AVI against lookups combined | The first draft's primary mixed the model-versus-rule difference with AVI's gain; "CADD for the rest" had no calibration (design review) |
| 2026-11-02 | scope | Cancer panel dropped, including the PP3 test that issue #2 keeps | Tim's decision: the rare-disease panel only |
| 2026-11-02 | 6, 10 | GPN-MSA, GPN-Star V and M kept as secondary scores, EVEE added as descriptive, although issue #2 drops GPN-Star | Tim's decision |
