# Guide to the studies

A visual way into this repository. What each study asks, why it was worth asking, how it works, what it found, and
what to read first. The READMEs in each folder have the full method and every number; this page is the map.

Every claim has a source. Published work is cited by first author and year, with the full reference and DOI in the
[reading list](#5-reading-list-in-order). Our own numbers link to the results file they come from.

Contents:

1. [The question in one picture](#1-the-question-in-one-picture)
2. [Primer: the ideas the studies rest on](#2-primer-the-ideas-the-studies-rest-on)
3. [How the studies led to each other](#3-how-the-studies-led-to-each-other)
4. [The studies, one card each](#4-the-studies-one-card-each)
5. [Reading list, in order](#5-reading-list-in-order)
6. [Glossary](#6-glossary)

## 1. The question in one picture

A genome differs from the reference at 4 to 5 million sites (1000 Genomes 2015). A diagnostic pipeline filters
them by gene panel, frequency and consequence to a short list of candidates, then a scientist reviews them in order
(Smedley 2015; Whiffin 2017). Scores decide that order. In our test cases the list has a median of 235 rare variants
([`rare-disease/results/results.json`](rare-disease/results/results.json)). Genomic foundation models promise better
scores. Do they put the diagnosis nearer the top?

```mermaid
flowchart LR
    g["Patient genome<br/>4 to 5 million variants"] --> f["Pipeline filters<br/>panel genes, rare,<br/>protein or splice effect"]
    f --> c["Candidate list<br/>a few hundred variants"]
    c --> s["Sort by scores"]
    s --> r["Scientist reviews<br/>from the top"]
    r --> d["Diagnosis<br/>(if it is near the top)"]

    old["Established scores<br/>CADD, SpliceAI,<br/>REVEL, AlphaMissense"] -.-> s
    new["Foundation-model scores<br/>GPN-MSA, GPN-Star,<br/>AlphaGenome AVI, Evo 2"] -.-> s

    style new fill:#fde2ef,stroke:#df007d
    style old fill:#e3eefa,stroke:#005eb8
```

Every study here is a version of one test: **does adding a foundation-model score move the diagnostic variant
up the list, compared with the scores labs already use?**

## 2. Primer: the ideas the studies rest on

### 2.1 A variant score

A variant effect predictor gives each variant a number: higher means more likely to damage a gene. Labs use these
numbers as one line of evidence in ACMG/AMP classification, the PP3 and BP4 codes (Richards 2015; Pejaver 2022),
and to sort candidates (Smedley 2015).

```mermaid
flowchart LR
    v["A variant, for example<br/>a BRCA1 SNV"] --> p["Predictor"]
    p --> n["A number:<br/>higher = more damaging"]
    n --> u1["Sorts the<br/>candidate list"]
    n --> u2["Evidence for<br/>classification (PP3)"]
```

### 2.2 Two kinds of score, and two ways to get one

```mermaid
flowchart TB
    subgraph est["Established scores"]
        cadd["CADD<br/>many annotations,<br/>simulated vs observed variants"]
        splice["SpliceAI<br/>predicts splicing<br/>from sequence"]
        revel["REVEL<br/>ensemble of missense tools"]
        am["AlphaMissense<br/>built on AlphaFold,<br/>missense only"]
    end
    subgraph fm["Genomic foundation models"]
        gpn["GPN-MSA, GPN-Star<br/>DNA language models<br/>over many species"]
        avi["AlphaGenome AVI<br/>AlphaGenome effects plus<br/>AlphaMissense, conservation,<br/>consequence: one score"]
        evo["Evo 2 7B<br/>DNA language model,<br/>2.4 trillion bases"]
    end
    est --> look["Lookup: precomputed table,<br/>milliseconds per variant"]
    gpn --> look
    avi --> look
    evo --> run["Run it yourself:<br/>GPU, about 0.7 s per variant"]
    style fm fill:#fde2ef,stroke:#df007d
    style est fill:#e3eefa,stroke:#005eb8
```

Sources: CADD (Kircher 2014; Schubach 2024), SpliceAI (Jaganathan 2019), REVEL (Ioannidis 2016), AlphaMissense
(Cheng 2023), GPN-MSA (Benegas 2025), GPN-Star (Ye 2026), AlphaGenome and AVI (Avsec 2026; Cheng 2026), Evo 2
(Brixi 2026). The 0.7 s per variant is our own H100 run: 2.9 windows per second, two windows per variant
([`rare-disease/evo2/README.md`](rare-disease/evo2/README.md)).

A **foundation model** is trained on DNA sequence, and for AlphaGenome on measured molecular data, rather than on
disease labels, and is then used for many tasks. Some publish genome-wide **score tables** (GPN-MSA, GPN-Star, AVI),
so a lab can look them up like CADD. AVI is not a pure foundation-model score: it combines AlphaGenome's predicted
effects with AlphaMissense, conservation and consequence, using a model trained on gnomAD rarity (Cheng 2026). Evo 2
has no genome-wide table. EVEE publishes scores built on Evo 2, but only for the about 4.25 million variants already
in ClinVar (Pearce 2026), so it cannot score most of a patient's rare variants. For those, a lab has to run Evo 2
itself, which raises cost and reproducibility questions.

### 2.3 Why rank cases, not variants

Most papers report how well a score separates pathogenic from benign variants across a whole dataset (AUROC). A lab
asks something different: within one patient, is the diagnosis near the top? A score can look good on AUROC by
learning that nonsense variants tend to be pathogenic, yet not help inside a case where the candidates already
passed a consequence filter. In ClinVar, 98% of nonsense SNVs are pathogenic and 0.1% of synonymous ones
(Pearce 2026). Benchmarks can also be inflated by training-data overlap (Grimm 2015) and by differences between
genes (Lu 2026).

```mermaid
flowchart LR
    subgraph v["Variant-level test"]
        direction TB
        a["Thousands of variants<br/>pathogenic vs benign"] --> au["AUROC"]
    end
    subgraph c["Case-level test"]
        direction TB
        b["One patient's ~235 rare variants<br/>1 is the diagnosis"] --> rk["Rank of the diagnosis"] --> tk["Top 1 / top 5 / top 10<br/>across many cases"]
    end
    v -.->|"can mislead because<br/>consequence alone<br/>predicts labels"| c
```

**Top k** is the share of cases in which the diagnostic variant is in the first k positions. Top 1 asks whether the
score finds the answer outright; top 5 asks whether the answer lands on the first screen.

### 2.4 Combining scores, and "what does it add?"

No lab uses one score alone. The studies fit a simple model (logistic regression) on the established scores plus
the variant's consequence class: "lookups combined". Then they add one foundation-model score and measure the
change. That change, not the score on its own, is the result that matters to a lab.

```mermaid
flowchart LR
    l["Lookups combined<br/>CADD + SpliceAI + REVEL<br/>+ AlphaMissense + consequence"] --> b["Baseline<br/>top 1 = 59.5%"]
    l2["Lookups combined<br/>+ one new score"] --> n["New<br/>top 1 = 67.0%"]
    b & n --> g["Gain = +7.5 points<br/>with a 95% interval"]
```

The model is fitted with cross-validation grouped by person, so no case is scored by a model that saw it. Intervals
come from resampling the background people (a bootstrap). Numbers: 200 SNV cases,
[`rare-disease/results/results.json`](rare-disease/results/results.json); method in
[`rare-disease/README.md`](rare-disease/README.md).

### 2.5 Four ways a benchmark can fool you, and the check for each

```mermaid
flowchart TB
    t1["Circularity:<br/>the score was trained on the<br/>same labels it is tested on"] --> c1["Name it as a limit;<br/>test AVI on chromosomes<br/>it was not trained on"]
    t2["Easy backgrounds:<br/>common variants a pipeline<br/>would have removed"] --> c2["Rerun with a gnomAD<br/>frequency filter"]
    t3["Winner's curse:<br/>the best of several scores<br/>on 200 cases looks too good"] --> c3["Replicate on 500 new cases<br/>with the model locked"]
    t4["Weak baseline:<br/>a new score 'helps' a<br/>combination that was poor"] --> c4["Compare with all<br/>established scores<br/>combined, fitted the same way"]
```

Circularity: Grimm 2015. AVI's training chromosomes: Cheng 2026; our held-out check in
[`rare-disease/REPLICATION.md`](rare-disease/REPLICATION.md). Frequency filtering: Whiffin 2017.

## 3. How the studies led to each other

```mermaid
flowchart TB
    q["Question: do genomic foundation models<br/>help rank a case's candidates?"]
    q --> b["BRCA1 case ranking<br/>one gene, assay labels<br/>Evo 2 adds little beyond CADD + SpliceAI"]
    q --> p["Evo 2 pipeline<br/>what running Evo 2 in a lab<br/>pipeline involves"]
    b --> rd["Rare-disease benchmark<br/>real backgrounds, ClinVar diagnoses<br/>SNVs, indels, inherited cancer"]
    rd --> rep["Replication<br/>500 new cases<br/>GPN-MSA gain does not hold"]
    rep --> avi["AlphaGenome AVI added<br/>gain holds on new cases"]
    avi --> held["Held-out chromosomes<br/>AVI lead unchanged"]
    avi --> gn["gnomAD frequency filter<br/>AVI gain holds; indel gain halves"]
    ev["EVEE charts<br/>why variant-level AUROC misleads"] -.-> rd
    style rd fill:#fde2ef,stroke:#df007d
```

## 4. The studies, one card each

### 4.1 Rare-disease case benchmark (the main study)

> **Question.** Put one known disease-causing variant into a real person's rare variants. How often does each score
> rank it first, or in the top 5? Does adding a foundation-model score to the established lookups help?
>
> **Why.** This is the closest public approximation of a lab's candidate list. It uses real backgrounds, real
> ClinVar diagnoses and the scores labs already have.
>
> **Folder.** [`rare-disease/`](rare-disease/), with the follow-up checks in
> [`rare-disease/REPLICATION.md`](rare-disease/REPLICATION.md).

```mermaid
flowchart LR
    person["1000 Genomes person"] --> rare["Rare variants in<br/>panel genes (~235)"]
    cv["ClinVar P/LP, 2+ stars"] --> dx["Diagnostic variant"]
    rare & dx --> case["Case"]
    case --> scores["9 scores and<br/>combined models"]
    scores --> rank["Rank of the diagnosis"]
```

Four studies share this design: 200 SNV cases, 120 indel cases, 140 inherited cancer cases, and a 500-case
replication with 50 new people ([`rare-disease/README.md`](rare-disease/README.md)). Backgrounds come from 1000
Genomes (Byrska-Bishop 2022), genes from PanelApp (Martin 2019) on MANE transcripts (Morales 2022), and diagnoses
from ClinVar (Landrum 2018).

Each single score on its own, 200 SNV cases, share of cases with the diagnosis in the top 5:

```mermaid
xychart-beta
    title "Top 5, single scores, 200 SNV cases (%)"
    x-axis ["REVEL", "AM", "Evo 2", "GPN-MSA", "CADD", "SpliceAI", "GPN-Star", "AVI", "Lookups"]
    y-axis "Cases with diagnosis in top 5 (%)" 0 --> 100
    bar [27, 26.5, 43, 48.5, 51.5, 54.5, 57.5, 87, 92]
```

AM is AlphaMissense, GPN-Star the vertebrate model, and Lookups the established scores combined. REVEL and
AlphaMissense score only missense variants, so they miss every other diagnosis. AVI is the best single
score; the established lookups combined are better still. Source:
[`rare-disease/results/results.json`](rare-disease/results/results.json), `strata.all.scores`.

What adding AlphaGenome AVI to the lookups gains at top 1, as the checks got stricter:

```mermaid
xychart-beta
    title "AVI added to the lookups: gain at top 1 (points)"
    x-axis ["A", "B", "C", "D", "E"]
    y-axis "Gain (points)" 0 --> 10
    bar [7.5, 6.2, 8.0, 7.4, 0.8]
```

A: the first 200 cases. B: 500 new cases, model locked. C: 200 cases with a gnomAD AF filter below 0.1%. D: 500 new
cases, locked, with that filter. E: 500 new cases, locked, with the strict grpmax filter below 0.1%.

Sources: A, [`results.json`](rare-disease/results/results.json); B,
[`results_replication.json`](rare-disease/results/results_replication.json); C,
[`results_gnomad0.001.json`](rare-disease/results/results_gnomad0.001.json); D,
[`results_replication_gnomad0.001.json`](rare-disease/results/results_replication_gnomad0.001.json); E,
[`results_replication_gnomad0.001_grpmax.json`](rare-disease/results/results_replication_gnomad0.001_grpmax.json).

The 95% intervals run 5 to 8 points either side of each bar, so read the bars together. The gain holds on new cases and with a pipeline's
frequency filter. Under a strict filter, lists shrink to a median of 34 to 49 variants. The lookups alone then rank
the diagnosis first in 70 to 79% of cases, which leaves little room for any score to add.

**Findings in one line each** (all from [`rare-disease/results/`](rare-disease/results/)):

- The established lookups, combined, do most of the work.
- AVI is the one foundation-model table whose gain holds on new cases.
- GPN-MSA's +6 point gain on the first 200 cases did not replicate.
- On indels, Evo 2 adds +16 points at top 1, or +8 after a frequency filter.
- On a realistic short list, little room is left for any new score to help.

**Read first:** Richards 2015 (how labs classify variants), Pejaver 2022 (score thresholds as evidence), Martin 2019
(PanelApp), Byrska-Bishop 2022 (the 1000 Genomes backgrounds), Cheng 2026 (the AlphaGenome Atlas preprint and its
own rare-disease test), Grimm 2015 (circularity).

### 4.2 Evo 2 on the case variants

> **Question.** What does it take to score every case variant with a model that has no genome-wide table, and does its score
> help?
>
> **Why.** Evo 2 is the largest openly released DNA foundation model (Brixi 2026), and a lab could only use it on
> a patient's variants by running it.
>
> **Folder.** [`rare-disease/evo2/`](rare-disease/evo2/)

```mermaid
flowchart LR
    v["16,357 case variants"] --> w["Two 8,192 bp windows each<br/>reference and variant"]
    w --> h["One rented H100<br/>pinned container"]
    h --> d["Score = change in<br/>mean log-likelihood"]
    d --> b["Into the case benchmark<br/>as 'Evo 2 7B'"]
    h --> rc["Receipts: GPU, image digest,<br/>timings, cost"]
```

Two runs took about 3.5 hours of pod time and cost about USD 14. Throughput was 2.9 windows per second
([`rare-disease/evo2/README.md`](rare-disease/evo2/README.md)). On SNVs, Evo 2 alone was among the weaker scores
(43% top 5, [`results.json`](rare-disease/results/results.json)). On indels, where few tables exist, it was the
strongest single score (92.5% top 5) and added the largest gain to the lookups
([`results_indel.json`](rare-disease/results/results_indel.json)).

**Read first:** Brixi 2026 (Evo 2), Arc Institute's
[BRCA1 zero-shot notebook](https://github.com/ArcInstitute/evo2/blob/main/notebooks/brca1/brca1_zero_shot_vep.ipynb)
(the scoring method followed here).

### 4.3 BRCA1 case ranking

> **Question.** On one gene with lab-measured function, does Evo 2 put the loss-of-function variant near the top of
> a simulated case?
>
> **Why.** It was the first case-level test. Saturation genome editing gives a functional label for 96.5% of
> possible SNVs in 13 BRCA1 exons (Findlay 2018), so it does not depend on ClinVar.
>
> **Folder.** [`brca1-case-ranking/`](brca1-case-ranking/)

```mermaid
flowchart LR
    s["3,893 BRCA1 SNVs<br/>with assay labels"] --> l["823 loss of function"]
    s --> f["3,070 functional<br/>or intermediate"]
    l --> c["Case: 1 LOF<br/>+ 100 functional"]
    f --> c
    c --> r["Rank: CADD, SpliceAI,<br/>REVEL, Evo 2, combinations"]
```

```mermaid
xychart-beta
    title "BRCA1 cases (K = 100): loss-of-function variant in the top 5 (%)"
    x-axis ["REVEL", "SpliceAI", "Evo 2", "CADD", "C+S", "C+S+Evo 2"]
    y-axis "Cases (%)" 0 --> 100
    bar [30.4, 37.2, 48.1, 54.9, 62.0, 65.0]
```

C+S is the type-aware CADD + SpliceAI combination. Evo 2 adds about 3 points to it, but against the best existing score for each variant type no gain is detected. The study
also showed that how scores are combined matters as much as which scores are combined. That lesson shaped the
rare-disease design. Source: [`brca1-case-ranking/README.md`](brca1-case-ranking/README.md) and
[`results.json`](brca1-case-ranking/results.json).

**Read first:** Findlay 2018 (BRCA1 saturation genome editing).

### 4.4 Evo 2 pipeline

> **Question.** What does adding an Evo 2 score to a lab pipeline involve in practice?
>
> **Why.** A score is only useful if it can run inside a validated, reproducible pipeline.
>
> **Folder.** [`evo2-pipeline/`](evo2-pipeline/)

```mermaid
flowchart LR
    vcf["Case VCF"] --> n["Normalise"] --> a["Add CADD,<br/>SpliceAI"] --> w["Windows"] --> s["Evo 2<br/>(GPU or replay)"] --> r["Rank the case"]
```

It is a Nextflow pipeline with two profiles. `replay` runs on a laptop with recorded H100 scores and has been run
end to end. `gpu` uses the same scoring code but has not been run end to end through Nextflow. The GPU, cost and
reproducibility findings shown in the talk, including scores that drift between GPU types, come from the separate
[evo2-cloud-operations](https://github.com/rewire-bio/evo2-cloud-operations) study.

Source: [`evo2-pipeline/README.md`](evo2-pipeline/README.md).

**Read first:** Di Tommaso 2017 (Nextflow), NVIDIA's
[cuBLAS results reproducibility notes](https://docs.nvidia.com/cuda/cublas/index.html#results-reproducibility) (why
GPU results can differ between GPU types).

### 4.5 EVEE charts

> **Question.** None of our own. This folder redraws two figures from Goodfire's EVEE preprint for the slides.
>
> **Why they matter here.** They show the trap from section 2.3. In ClinVar, a variant's consequence alone predicts
> much of its label, so variant-level benchmarks reward scores that learn consequence. Within each consequence,
> AUROCs are lower (for CADD v1.7, 0.99 overall against about 0.79 to 0.95 within a consequence) and the gaps between
> methods change (Pearce 2026; data from
> [goodfire-ai/evee-manuscript](https://github.com/goodfire-ai/evee-manuscript)).
>
> **Folder.** [`evee/`](evee/)

```mermaid
flowchart LR
    cv["ClinVar SNVs<br/>about 834,000"] --> con["By consequence:<br/>nonsense 98% pathogenic,<br/>synonymous 0.1%"]
    con --> trap["A score can win on AUROC<br/>by learning consequence"]
    trap --> fix["EVEE: compare within<br/>each consequence"]
```

**Read first:** Pearce 2026 (EVEE).

## 5. Reading list, in order

A path from the clinical setting to the models, then to the evaluation problems. The bold items are the ones to read
if time is short.

```mermaid
flowchart LR
    A["1. How labs classify<br/>Richards 2015<br/>Pejaver 2022"] --> B["2. The inputs<br/>PanelApp, ClinVar,<br/>1000 Genomes, gnomAD"]
    B --> C["3. Established scores<br/>CADD, SpliceAI,<br/>REVEL, AlphaMissense"]
    C --> D["4. Foundation models<br/>GPN-MSA, GPN-Star,<br/>AlphaGenome, Evo 2"]
    D --> E["5. Evaluation traps<br/>Grimm 2015, EVEE,<br/>Lu 2026"]
```

**1. How labs classify variants**

- **Richards S et al. (2015).** Standards and guidelines for the interpretation of sequence variants (ACMG/AMP).
  *Genetics in Medicine* 17:405-424. [doi:10.1038/gim.2015.30](https://doi.org/10.1038/gim.2015.30)
- **Pejaver V et al. (2022).** Calibration of computational tools for missense variant pathogenicity
  classification and ClinGen recommendations for PP3/BP4 criteria. *AJHG* 109:2163-2177.
  [doi:10.1016/j.ajhg.2022.10.013](https://doi.org/10.1016/j.ajhg.2022.10.013)
- Walker LC et al. (2023). Using the ACMG/AMP framework to capture evidence related to predicted and observed impact
  on splicing: recommendations from the ClinGen SVI Splicing Subgroup. *AJHG* 110:1046-1067. [doi:10.1016/j.ajhg.2023.06.002](https://doi.org/10.1016/j.ajhg.2023.06.002)
- Bergquist T et al. (2025). Calibration of additional computational tools expands ClinGen recommendation options
  for variant classification with PP3/BP4 criteria. *Genetics in Medicine* 27:101402.
  [doi:10.1016/j.gim.2025.101402](https://doi.org/10.1016/j.gim.2025.101402)
- Smedley D et al. (2015). Next-generation diagnostics and disease-gene discovery with the Exomiser. *Nature
  Protocols* 10:2004-2015. [doi:10.1038/nprot.2015.124](https://doi.org/10.1038/nprot.2015.124). How a real
  pipeline ranks candidates.

**2. The inputs the cases are built from**

- **Martin AR et al. (2019).** PanelApp crowdsources expert knowledge to establish consensus diagnostic gene panels.
  *Nature Genetics* 51:1560-1565. [doi:10.1038/s41588-019-0528-2](https://doi.org/10.1038/s41588-019-0528-2)
- 1000 Genomes Project Consortium (2015). A global reference for human genetic variation. *Nature* 526:68-74.
  [doi:10.1038/nature15393](https://doi.org/10.1038/nature15393). Source of "4 to 5 million variants per genome".
- Byrska-Bishop M et al. (2022). High-coverage whole-genome sequencing of the expanded 1000 Genomes Project cohort
  including 602 trios. *Cell* 185:3426-3440.e19. [doi:10.1016/j.cell.2022.08.004](https://doi.org/10.1016/j.cell.2022.08.004)
- Landrum MJ et al. (2018). ClinVar: improving access to variant interpretations and supporting evidence. *Nucleic
  Acids Research* 46:D1062-D1067. [doi:10.1093/nar/gkx1153](https://doi.org/10.1093/nar/gkx1153). Review stars:
  [ClinVar review status](https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/).
- Chen S et al. (2024). A genomic mutational constraint map using variation in 76,156 human genomes (gnomAD).
  *Nature* 625:92-100. [doi:10.1038/s41586-023-06045-0](https://doi.org/10.1038/s41586-023-06045-0). The grpmax
  definition: [gnomAD help](https://gnomad.broadinstitute.org/help/grpmax).
- Whiffin N et al. (2017). Using high-resolution variant frequencies to empower clinical genome interpretation.
  *Genetics in Medicine* 19:1151-1158. [doi:10.1038/gim.2017.26](https://doi.org/10.1038/gim.2017.26). Why
  pipelines filter on frequency first.
- Findlay GM et al. (2018). Accurate classification of BRCA1 variants with saturation genome editing. *Nature*
  562:217-222. [doi:10.1038/s41586-018-0461-z](https://doi.org/10.1038/s41586-018-0461-z). The BRCA1 assay labels.
- Morales J et al. (2022). A joint NCBI and EMBL-EBI transcript set for clinical genomics and research (MANE).
  *Nature* 604:310-315. [doi:10.1038/s41586-022-04558-8](https://doi.org/10.1038/s41586-022-04558-8)

**3. Established scores**

- Kircher M et al. (2014). A general framework for estimating the relative pathogenicity of human genetic variants
  (CADD). *Nature Genetics* 46:310-315. [doi:10.1038/ng.2892](https://doi.org/10.1038/ng.2892)
- **Schubach M et al. (2024).** CADD v1.7: using protein language models, regulatory CNNs and other
  nucleotide-level scores to improve genome-wide variant predictions. *Nucleic Acids Research* 52:D1143-D1154.
  [doi:10.1093/nar/gkad989](https://doi.org/10.1093/nar/gkad989)
- **Jaganathan K et al. (2019).** Predicting splicing from primary sequence with deep learning (SpliceAI). *Cell*
  176:535-548.e24. [doi:10.1016/j.cell.2018.12.015](https://doi.org/10.1016/j.cell.2018.12.015)
- Ioannidis NM et al. (2016). REVEL: an ensemble method for predicting the pathogenicity of rare missense variants.
  *AJHG* 99:877-885. [doi:10.1016/j.ajhg.2016.08.016](https://doi.org/10.1016/j.ajhg.2016.08.016)
- **Cheng J et al. (2023).** Accurate proteome-wide missense variant effect prediction with AlphaMissense. *Science*
  381:eadg7492. [doi:10.1126/science.adg7492](https://doi.org/10.1126/science.adg7492)

**4. Genomic foundation models**

- **Benegas G et al. (2025).** A DNA language model based on multispecies alignment predicts the effects of
  genome-wide variants (GPN-MSA). *Nature Biotechnology* 43:1960-1965.
  [doi:10.1038/s41587-024-02511-w](https://doi.org/10.1038/s41587-024-02511-w)
- Ye C et al. (2026). Predicting genome-wide functional constraints with GPN-Star. *Nature*.
  [doi:10.1038/s41586-026-11005-5](https://doi.org/10.1038/s41586-026-11005-5)
- **Avsec Ž et al. (2026).** Advancing regulatory variant effect prediction with AlphaGenome. *Nature*
  649:1206-1218. [doi:10.1038/s41586-025-10014-0](https://doi.org/10.1038/s41586-025-10014-0)
- **Cheng J et al. (2026).** AlphaGenome Atlas: in silico mutagenesis of the entire human genome improves
  prioritization and interpretation of non-coding variants. *medRxiv*.
  [doi:10.64898/2026.09.16.26363192](https://doi.org/10.64898/2026.09.16.26363192). The source of the AVI table,
  its training chromosomes and the 0.1% frequency filter used here.
- **Brixi G et al. (2026).** Genome modelling and design across all domains of life with Evo 2. *Nature*
  652:1349-1361. [doi:10.1038/s41586-026-10176-5](https://doi.org/10.1038/s41586-026-10176-5)

**5. Evaluation traps**

- **Grimm DG et al. (2015).** The evaluation of tools used to predict the impact of missense variants is hindered
  by two types of circularity. *Human Mutation* 36:513-523.
  [doi:10.1002/humu.22768](https://doi.org/10.1002/humu.22768)
- **Pearce MT et al. (2026).** Interpretable variant effect prediction from genomic foundation model
  representations (EVEE). *bioRxiv*. [doi:10.64898/2026.04.10.717844](https://doi.org/10.64898/2026.04.10.717844)
- Lu B et al. (2026). Genomic heterogeneity inflates the performance of variant pathogenicity predictions.
  *Nature Communications*, online 30 October 2026.
  [doi:10.1038/s41467-026-78623-5](https://doi.org/10.1038/s41467-026-78623-5). Preprint:
  [doi:10.1101/2025.09.05.674459](https://doi.org/10.1101/2025.09.05.674459)
- Kong SW et al. (2025). Discordance between a deep learning model and clinical-grade variant pathogenicity
  classification in a rare disease cohort. *npj Genomic Medicine* 10:17.
  [doi:10.1038/s41525-025-00480-w](https://doi.org/10.1038/s41525-025-00480-w)
- Curtis D (2025). Assessment of ability of a DNA language model to predict pathogenicity of rare coding variants.
  *Journal of Human Genetics* 70:603-607.
  [doi:10.1038/s10038-025-01385-3](https://doi.org/10.1038/s10038-025-01385-3)
- Di Tommaso P et al. (2017). Nextflow enables reproducible computational workflows. *Nature Biotechnology*
  35:316-319. [doi:10.1038/nbt.3820](https://doi.org/10.1038/nbt.3820)

## 6. Glossary

| Term | Meaning |
|---|---|
| ACMG/AMP | The standard framework for classifying variants into five classes: pathogenic, likely pathogenic, uncertain significance, likely benign, benign (Richards 2015) |
| PP3, BP4 | ACMG/AMP evidence codes for computational predictions supporting pathogenic or benign |
| P/LP | Pathogenic or likely pathogenic |
| ClinVar review stars | How well supported a ClinVar classification is. 2 stars: two or more submitters with criteria and no conflict; 3: expert panel; 4: practice guideline |
| Background | The rare variants of one 1000 Genomes person, into which a diagnostic variant is inserted |
| Top k | Share of cases with the diagnostic variant in the first k positions |
| Lookups combined | Logistic regression on CADD, SpliceAI, REVEL, AlphaMissense and consequence class |
| Locked | A model fitted on the first 200 cases and applied unchanged to the 500 new ones |
| Bootstrap interval | A 95% range from resampling the background people 1,000 times |
| AF, grpmax AF | Allele frequency in gnomAD. gnomAD's grpmax is the highest frequency among its large, non-bottlenecked ancestry groups. Our filter takes the highest of afr, amr, eas, mid, nfe and sas, which for genomes adds mid to gnomAD's set |
| Zero-shot | Using a model's own likelihoods as a score, without training on labels |
| Log-likelihood delta | How much less likely the model finds the variant sequence than the reference |
| SGE | Saturation genome editing: nearly every possible SNV in a region made in cells and its effect measured (Findlay 2018) |
| AUROC | How well a score separates two groups across all thresholds; 1 is perfect, 0.5 is chance |
