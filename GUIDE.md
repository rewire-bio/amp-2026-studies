# Guide to the studies

A visual way into this repository. What each study asks, why it was worth asking, how it works, what it found, and
what to read first. The READMEs in each folder have the full method and every number; this page is the map.

Contents:

1. [The question in one picture](#1-the-question-in-one-picture)
2. [Primer: the ideas the studies rest on](#2-primer-the-ideas-the-studies-rest-on)
3. [How the studies led to each other](#3-how-the-studies-led-to-each-other)
4. [The studies, one card each](#4-the-studies-one-card-each)
5. [Reading list, in order](#5-reading-list-in-order)
6. [Glossary](#6-glossary)

## 1. The question in one picture

A rare disease genome has millions of variants. A diagnostic pipeline filters them to a few hundred candidates, then
a scientist reviews them in order. Scores decide that order. Genomic foundation models promise better scores. Do
they put the diagnosis nearer the top?

```mermaid
flowchart LR
    g["Patient genome<br/>~5 million variants"] --> f["Pipeline filters<br/>panel genes, rare,<br/>protein or splice effect"]
    f --> c["Candidate list<br/>~100 to 300 variants"]
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
numbers as one line of evidence in ACMG/AMP classification (PP3 and BP4) and to sort candidates.

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
        am["AlphaMissense<br/>protein structure model,<br/>missense only"]
    end
    subgraph fm["Genomic foundation models"]
        gpn["GPN-MSA, GPN-Star<br/>DNA language models<br/>over many species"]
        avi["AlphaGenome AVI<br/>predicted molecular effects,<br/>combined into one score"]
        evo["Evo 2<br/>DNA language model,<br/>9.3 trillion bases"]
    end
    est --> look["Lookup: precomputed table,<br/>milliseconds per variant"]
    gpn --> look
    avi --> look
    evo --> run["Run it yourself:<br/>GPU, about 0.7 s per variant"]
    style fm fill:#fde2ef,stroke:#df007d
    style est fill:#e3eefa,stroke:#005eb8
```

A **foundation model** learns from raw DNA sequence without disease labels and is then used for many tasks. Some
publish genome-wide **score tables** (GPN-MSA, GPN-Star, AVI), so a lab can look them up like CADD. Evo 2 has no
table, so a lab has to run it, which raises cost and reproducibility questions.

### 2.3 Why rank cases, not variants

Most papers report how well a score separates pathogenic from benign variants across a whole dataset (AUROC). A lab
asks something different: within one patient, is the diagnosis near the top? A score can look good on AUROC by
learning that nonsense variants tend to be pathogenic, yet not help inside a case where the candidates already
passed a consequence filter.

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
come from resampling the background people (a bootstrap).

### 2.5 Four ways a benchmark can fool you, and the check for each

```mermaid
flowchart TB
    t1["Circularity:<br/>the score was trained on the<br/>same labels it is tested on"] --> c1["Name it as a limit;<br/>test AVI on chromosomes<br/>it was not trained on"]
    t2["Easy backgrounds:<br/>common variants a pipeline<br/>would have removed"] --> c2["Rerun with a gnomAD<br/>frequency filter"]
    t3["Winner's curse:<br/>the best of several scores<br/>on 200 cases looks too good"] --> c3["Replicate on 500 new cases<br/>with the model locked"]
    t4["Weak baseline:<br/>a new score 'helps' a<br/>combination that was poor"] --> c4["Compare with the best<br/>existing score per type"]
```

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
replication with 50 new people.

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
score; the established lookups combined are better still.

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

The 95% intervals are about ±7 points, so read the bars together. The gain holds on new cases and with a pipeline's
frequency filter. Under a strict filter, lists shrink to a median of 34 to 49 variants. The lookups alone then rank
the diagnosis first in 70 to 79% of cases, which leaves little room for any score to add.

**Findings in one line each:**

- The established lookups, combined, do most of the work.
- AVI is the one foundation-model table whose gain holds on new cases.
- GPN-MSA's +6 point gain on the first 200 cases did not replicate.
- On indels, Evo 2 adds +16 points at top 1, or +8 after a frequency filter.
- On a realistic short list, little room is left for any new score to help.

**Read first:** Richards 2015 (how labs classify variants), Pejaver 2022 (score thresholds as evidence), Martin 2019
(PanelApp), Byrska-Bishop 2022 (the 1000 Genomes backgrounds), Cheng 2026 (the AlphaGenome Atlas preprint and its
own rare-disease test), Grimm 2015 (circularity).

### 4.2 Evo 2 on the case variants

> **Question.** What does it take to score every case variant with a model that has no table, and does its score
> help?
>
> **Why.** Evo 2 is the best-known DNA foundation model, and a lab could only use it by running it.
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

Two runs took about 3.5 hours of pod time and cost about USD 14. Throughput was 2.9 windows per second. On SNVs,
Evo 2 alone was among the weaker scores (43% top 5). On indels, where few tables exist, it was the strongest single
score (92.5% top 5) and added the largest gain to the lookups.

**Read first:** Brixi 2026 (Evo 2), Arc Institute's BRCA1 zero-shot notebook (the scoring method followed here).

### 4.3 BRCA1 case ranking

> **Question.** On one gene with lab-measured function, does Evo 2 put the loss-of-function variant near the top of
> a simulated case?
>
> **Why.** It was the first case-level test. Saturation genome editing gives a functional label for almost every
> possible SNV in the tested BRCA1 exons, so it does not depend on ClinVar.
>
> **Folder.** [`brca1-case-ranking/`](brca1-case-ranking/)

```mermaid
flowchart LR
    s["3,893 BRCA1 SNVs<br/>with assay labels"] --> l["823 loss of function"]
    s --> f["3,070 functional"]
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
rare-disease design.

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

**Read first:** Di Tommaso 2017 (Nextflow), the cuBLAS reproducibility notes (why GPU results can differ).

### 4.5 EVEE charts

> **Question.** None of our own. This folder redraws two figures from Goodfire's EVEE preprint for the slides.
>
> **Why they matter here.** They show the trap from section 2.3. In ClinVar, a variant's consequence alone predicts
> much of its label, so variant-level benchmarks reward scores that learn consequence. Within each consequence the
> scores separate much less well.
>
> **Folder.** [`evee/`](evee/)

```mermaid
flowchart LR
    cv["ClinVar SNVs<br/>834,000"] --> con["By consequence:<br/>nonsense 98% pathogenic,<br/>synonymous 0.1%"]
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
    D --> E["5. Evaluation traps<br/>Grimm 2015, EVEE,<br/>Lu 2025"]
```

**1. How labs classify variants**

- **Richards S et al. (2015).** Standards and guidelines for the interpretation of sequence variants (ACMG/AMP).
  *Genetics in Medicine* 17:405-424. [doi:10.1038/gim.2015.30](https://doi.org/10.1038/gim.2015.30)
- **Pejaver V et al. (2022).** Calibration of computational tools for missense variant pathogenicity
  classification and ClinGen recommendations for PP3/BP4 criteria. *AJHG* 109:2163-2177.
  [doi:10.1016/j.ajhg.2022.10.013](https://doi.org/10.1016/j.ajhg.2022.10.013)
- Walker LC et al. (2023). Using the ACMG/AMP framework to capture evidence related to predicted and observed impact
  on splicing. *AJHG* 110:1046-1067. [doi:10.1016/j.ajhg.2023.06.002](https://doi.org/10.1016/j.ajhg.2023.06.002)
- Bergquist T et al. (2025). Calibration of additional computational tools expands ClinGen recommendation options
  for variant classification. *Genetics in Medicine* 27:101402.
  [doi:10.1016/j.gim.2025.101402](https://doi.org/10.1016/j.gim.2025.101402)
- Smedley D et al. (2015). Next-generation diagnostics and disease-gene discovery with the Exomiser. *Nature
  Protocols* 10:2004-2015. [doi:10.1038/nprot.2015.124](https://doi.org/10.1038/nprot.2015.124). How a real
  pipeline ranks candidates.

**2. The inputs the cases are built from**

- **Martin AR et al. (2019).** PanelApp crowdsources expert knowledge to establish consensus diagnostic gene panels.
  *Nature Genetics* 51:1560-1565. [doi:10.1038/s41588-019-0528-2](https://doi.org/10.1038/s41588-019-0528-2)
- Byrska-Bishop M et al. (2022). High-coverage whole-genome sequencing of the expanded 1000 Genomes Project cohort.
  *Cell* 185:3426-3440. [doi:10.1016/j.cell.2022.08.004](https://doi.org/10.1016/j.cell.2022.08.004)
- Landrum MJ et al. (2018). ClinVar: improving access to variant interpretations and supporting evidence. *Nucleic
  Acids Research* 46:D1062-D1067. [doi:10.1093/nar/gkx1153](https://doi.org/10.1093/nar/gkx1153)
- Chen S et al. (2024). A genomic mutational constraint map using variation in 76,156 human genomes (gnomAD).
  *Nature* 625:92-100. [doi:10.1038/s41586-023-06045-0](https://doi.org/10.1038/s41586-023-06045-0)
- Whiffin N et al. (2017). Using high-resolution variant frequencies to empower clinical genome interpretation.
  *Genetics in Medicine* 19:1151-1158. [doi:10.1038/gim.2017.26](https://doi.org/10.1038/gim.2017.26). Why
  pipelines filter on frequency first.
- Morales J et al. (2022). A joint NCBI and EMBL-EBI transcript set for clinical genomics and research (MANE).
  *Nature* 604:310-315. [doi:10.1038/s41586-022-04558-8](https://doi.org/10.1038/s41586-022-04558-8)

**3. Established scores**

- Kircher M et al. (2014). A general framework for estimating the relative pathogenicity of human genetic variants
  (CADD). *Nature Genetics* 46:310-315. [doi:10.1038/ng.2892](https://doi.org/10.1038/ng.2892)
- **Schubach M et al. (2024).** CADD v1.7. *Nucleic Acids Research* 52:D1143-D1154.
  [doi:10.1093/nar/gkad989](https://doi.org/10.1093/nar/gkad989)
- **Jaganathan K et al. (2019).** Predicting splicing from primary sequence with deep learning (SpliceAI). *Cell*
  176:535-548. [doi:10.1016/j.cell.2018.12.015](https://doi.org/10.1016/j.cell.2018.12.015)
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
- **Cheng J et al. (2026).** AlphaGenome Atlas: in silico mutagenesis of the entire human genome. *medRxiv*.
  [doi:10.64898/2026.09.16.26363192](https://doi.org/10.64898/2026.09.16.26363192). The source of the AVI table,
  its training chromosomes and the 0.1% frequency filter used here.
- **Brixi G et al. (2026).** Genome modelling and design across all domains of life with Evo 2. *Nature*
  652:1349-1361. [doi:10.1038/s41586-026-10176-5](https://doi.org/10.1038/s41586-026-10176-5)
- Findlay GM et al. (2018). Accurate classification of BRCA1 variants with saturation genome editing. *Nature*
  562:217-222. [doi:10.1038/s41586-018-0461-z](https://doi.org/10.1038/s41586-018-0461-z)

**5. Evaluation traps**

- **Grimm DG et al. (2015).** The evaluation of tools used to predict the impact of missense variants is hindered
  by two types of circularity. *Human Mutation* 36:513-523.
  [doi:10.1002/humu.22768](https://doi.org/10.1002/humu.22768)
- **Pearce MT et al. (2026).** Interpretable variant effect prediction from genomic foundation model
  representations (EVEE). *bioRxiv*. [doi:10.64898/2026.04.10.717844](https://doi.org/10.64898/2026.04.10.717844)
- Lu B et al. (2025). Genomic heterogeneity inflates the performance of variant pathogenicity predictions.
  *bioRxiv*. [doi:10.1101/2025.09.05.674459](https://doi.org/10.1101/2025.09.05.674459)
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
| ACMG/AMP | The standard framework for classifying variants as pathogenic, benign or uncertain |
| PP3, BP4 | ACMG/AMP evidence codes for computational predictions supporting pathogenic or benign |
| P/LP | Pathogenic or likely pathogenic |
| ClinVar review stars | How well supported a ClinVar classification is; 2+ stars means multiple submitters agree or an expert panel |
| Background | The rare variants of one 1000 Genomes person, into which a diagnostic variant is inserted |
| Top k | Share of cases with the diagnostic variant in the first k positions |
| Lookups combined | Logistic regression on CADD, SpliceAI, REVEL, AlphaMissense and consequence class |
| Locked | A model fitted on the first 200 cases and applied unchanged to the 500 new ones |
| Bootstrap interval | A 95% range from resampling the background people 1,000 times |
| AF, grpmax AF | Allele frequency in gnomAD; grpmax is the highest frequency in any one ancestry group |
| Zero-shot | Using a model's own likelihoods as a score, without training on labels |
| Log-likelihood delta | How much less likely the model finds the variant sequence than the reference |
| SGE | Saturation genome editing: every possible SNV in a region made in cells and its effect measured |
| AUROC | How well a score separates two groups across all thresholds; 1 is perfect, 0.5 is chance |
