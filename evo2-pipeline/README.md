# Evo 2 case pipeline

A Nextflow DSL2 pipeline that adds an Evo 2 7B zero-shot score to a case's variants next to CADD and SpliceAI,
and ranks the case. Shown on the talk's `pipeline-overview` slide and the backup `code-*` slides.

```sh
cd evo2-pipeline
python3 scripts/prepare_resources.py --study ../../evo2-cloud-operations   # builds resources/ and test/
JAVA_HOME=/opt/homebrew/opt/openjdk@21 nextflow run main.nf -profile replay  # laptop, recorded Evo 2 scores
nextflow run main.nf -profile gpu                                             # GPU host with Docker
```

Needs Nextflow 24.10 or later (Java 17+), bcftools, tabix, bgzip and samtools. `resources/` is not committed:
CADD and the assay labels are licensed for non-commercial use. The replay profile has been run end to end; the gpu
profile uses the scoring code that produced the recorded H100 scores but has not been run end to end through
Nextflow on a GPU host.
