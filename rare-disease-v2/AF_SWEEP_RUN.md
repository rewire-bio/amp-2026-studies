# Running the frequency sweep

How to run the exploratory frequency sweep (PROTOCOL.md section 13) on the machine that holds the v6 working data
and the AVI SSD. The code, the protocol section, the tests and the chart code are on branch `af-sweep`; nothing in
this file has been run on real data yet.

The sweep reads a second build of the 200 SNV cases and the 500, with background variants below 5%, in `data/af5/`.
Every v6 file and result stays as it is: the scripts write to `data/af5/` and `results/af5/` through `RDV2_DATA` and
`RDV2_RESULTS`, and only the last steps write new files to `results/`.

Run every command from `rare-disease-v2/`, in order, in one shell (bash or zsh). Two functions stand for the
pinned environment, and for it pointed at the af5 copy:

```sh
run() { uv run -q --python 3.12 --with-requirements requirements.txt python "$@"; }
af5() { RDV2_DATA=data/af5 RDV2_RESULTS=results/af5 run "$@"; }
```

## 1. Tools and inputs

Tools on `PATH`:

- `bcftools` (`backgrounds.py` reads the 1000 Genomes BCFs, `build_cases.py` reads ClinVar, `gnomad_af.py` splits
  and queries the gnomAD rows) and `tabix` (CADD, AVI, GPN-MSA). On a Mac: `brew install bcftools samtools htslib`.
- `uv`, with the versions pinned in `requirements.txt` and Python 3.12 (v6 ran on 3.12.14). The 0.1% and 1% rows can
  only reproduce v6 exactly in that environment.
- `samtools` and SpliceAI are not needed for this sweep. In the SNV studies SpliceAI is CADD's precomputed
  `spliceai_max` (PROTOCOL.md section 6); `run_spliceai.sh` scores indels only, and the sweep ranks no indels.

Inputs, with paths relative to `rare-disease-v2/`. `..` is the repository root, and `../..` is the folder that
holds the repository and its sibling repositories (`common.BLOG`).

| Input | Path | Notes |
|---|---|---|
| 1000 Genomes genotypes, main set | `../rare-disease/data/kgp/chr*.bcf` | committed |
| 1000 Genomes genotypes, replication set | `data/rep/kgp/chr*.bcf` | committed |
| ClinVar VCF, fileDate 2026-10-11 | `../rare-disease/data/clinvar.vcf.gz` | not committed; background P/LP removal |
| gnomAD cache | `data/gnomad_subset/` (sample metadata, pedigree, `dense/`, `counts_chr*.tsv`, `dense/chr*.positions.txt`) | not committed; v6's cache, so only new positions are fetched |
| gnomAD full counts, people, modes, pools | `data/gnomad_full.tsv`, `data/people.tsv`, `data/gene_moi.tsv`, `data/pool/` | committed |
| MANE, PanelApp, regions | `../rare-disease/data/MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz`, `panelapp_285.json`, `region_genes.tsv` | committed |
| AVI table | `../rare-disease/data/avi/alphagenome_variant_impact_score_snvs.tsv.gz` and its `.tbi` | on the SSD; the `avi` link must resolve |
| AlphaMissense | `../rare-disease/data/am/AlphaMissense_hg38.tsv.gz` | not committed |
| REVEL v1.3 | `../../evo2-cloud-operations/data/raw/current-predictors/revel-v1.3_all_chromosomes.zip` | not committed |
| v6 CADD and REVEL scores | `data/snv/scores_cadd.tsv`, `scores_revel.tsv`; `data/replication/scores_cadd.tsv`, `scores_revel.tsv` | not committed; needed for the v6 results and reused by the lookups |
| Evo 2 scores | `evo2/run*/scores*.tsv` | committed; cover the v6 lists only |
| Network | CADD (kircherlab.bihealth.org), GPN-MSA and GPN-Star (Hugging Face), MyVariant.info, gnomAD on Google Cloud Storage | remote lookups |

Check them before starting:

```sh
git switch af-sweep
for t in bcftools tabix uv; do command -v $t >/dev/null || echo "missing tool: $t"; done
for f in ../rare-disease/data/clinvar.vcf.gz \
         ../rare-disease/data/avi/alphagenome_variant_impact_score_snvs.tsv.gz \
         ../rare-disease/data/am/AlphaMissense_hg38.tsv.gz \
         ../../evo2-cloud-operations/data/raw/current-predictors/revel-v1.3_all_chromosomes.zip \
         data/gnomad_subset/gnomad.genomes.v3.1.2.hgdp_1kg_subset_sample_meta.tsv.bgz \
         data/gnomad_subset/1kGP.3202_samples.pedigree_info.txt \
         data/snv/scores_cadd.tsv data/snv/scores_revel.tsv \
         data/replication/scores_cadd.tsv data/replication/scores_revel.tsv; do
  [ -s "$f" ] || echo "missing input: $f"
done
ls ../rare-disease/data/kgp/chr*.bcf data/rep/kgp/chr*.bcf | wc -l   # 46
ls data/gnomad_subset/counts_chr*.tsv | wc -l                        # 23
```

## 2. Check that v6 reproduces here

Run the four v6 analyses the sweep is compared with, into a scratch folder, and compare them with the committed
results. If they differ, stop: the environment is not v6's, and the sweep's reproduction check would fail for the
same reason.

```sh
mkdir -p results/af5/v6check logs
for s in snv replication; do for a in main bg1pct; do
  RDV2_RESULTS=results/af5/v6check run case_benchmark.py --study $s --analysis $a > logs/af5_v6check_${s}_${a}.log
  cmp results/af5/v6check/${s}_${a}.json results/${s}_${a}.json && echo "$s $a: identical"
done; done
```

## 3. Set up `data/af5/`

Copies of the inputs the build reads, so the build changes nothing under `data/` outside `data/af5/`. The gnomAD
cache is cloned (`cp -c` makes APFS clones, which take no space until changed), because the fetch adds part files
and redoes the counts of the chromosomes it touches.

```sh
mkdir -p data/af5/main data/af5/rep results/af5
cp data/people.tsv data/gene_moi.tsv data/gnomad_full.tsv data/af5/
cp -R data/pool data/af5/pool
cp data/rep/samples50.txt data/af5/rep/
ln -s ../../rep/kgp data/af5/rep/kgp
cp -Rc data/gnomad_subset data/af5/gnomad_subset
echo "$(date '+%F %T') af5 set up" >> logs/af_sweep_progress.log
```

## 4. Rebuild at the 5% cap

Each person's variants with 1000 Genomes AF below 5%, their family-left-out frequencies, then the cases: the same
people, diagnostics and seeds, with background grpmax below 5% and a `kgp_af` column. The fetch reads only the
positions v6 did not fetch (about an hour for all of v6's; less here). For a long step, add a line to
`logs/af_sweep_progress.log` about every 10 minutes, for example from a second shell with
`while sleep 600; do echo "$(date '+%F %T') $(tail -1 logs/af5_gnomad.log)" >> logs/af_sweep_progress.log; done`.

```sh
af5 backgrounds.py main rep --cap 0.05 > logs/af5_backgrounds.log
for step in full fetch counts af; do
  af5 gnomad_af.py $step >> logs/af5_gnomad.log 2>&1
  echo "$(date '+%F %T') gnomad_af.py $step done" >> logs/af_sweep_progress.log
done
af5 build_cases.py --study snv --cap 0.05 > logs/af5_build_snv.log
af5 build_cases.py --study replication --cap 0.05 > logs/af5_build_replication.log
echo "$(date '+%F %T') af5 cases built" >> logs/af_sweep_progress.log
```

`build_cases.py` writes `data/af5/<study>/cases.tsv` and `case_variants.tsv`, and its report to
`results/af5/build_<study>.md`. Its checks must all say yes. The snv study must be built first: the replication
study leaves out the snv study's diagnostic positions.

## 5. Score the case variants

Every lookup score for every key of the af5 lists. These runs leave `RDV2_DATA` unset, so the v6 score files in
`data/snv/` and `data/indel/` are reused and only new keys are looked up (the v6 replication keys are looked up
again from the same sources; step 6 checks that the values match). EVEE is not needed. SpliceAI comes with CADD.
Evo 2 is not run (it needs a paid GPU), so the sweep reports it at 1% and below only.

```sh
for s in snv replication; do
  run lookup_scores.py --variants data/af5/$s/case_variants.tsv --out-dir data/af5/$s \
      cadd revel am avi gpnmsa gpnstar > logs/af5_lookup_$s.log 2>&1
  echo "$(date '+%F %T') $s scored" >> logs/af_sweep_progress.log
done
```

Each source writes `data/af5/<study>/scores_<source>.tsv` (GPN-Star writes `scores_gpnstar_v.tsv` and
`scores_gpnstar_m.tsv`).

## 6. Check the af5 build against v6

Below the 1% caps (grpmax and 1000 Genomes AF), every af5 case must hold exactly the v6 case's rows, and every v6
case variant must have the same scores in both builds.

```sh
run check_af5.py > logs/af5_identity.log; echo "exit $?"
```

It writes `results/af5_identity.md`. Exit 0 means identical. Otherwise the file lists each differing case, the rows
found in one build only, and the other build's frequencies for the same person and key, which show the cause (a
changed frequency, a 1000 Genomes AF on the other side of 1%, or a variant missing from one build). Stop and report
any difference before going on: the sweep's 0.1% and 1% rows can only reproduce v6 if this check passes.

## 7. The sweep

```sh
for s in snv replication; do
  run case_benchmark.py --study $s --analysis sweep > logs/bench_${s}_sweep.log 2>&1
  echo "$(date '+%F %T') $s sweep done" >> logs/af_sweep_progress.log
done
grep "reproduces" logs/bench_*_sweep.log
```

It writes `results/snv_sweep.json`, `results/replication_sweep.json` and `results/analysis_sweep.md` (rewritten
after each study, so it holds both after the second).

## 8. The 0.1% and 1% reproduction check

The sweep compares its 0.1% and 1% rows with `results/<study>_main.json` and `results/<study>_bg1pct.json` on
every value both compute: cases, list length, coverage, ties, every score and contrast in all cases and in the
gnomAD strata, and the gnomAD split. The result is under `reproduction` in each sweep JSON and at the end of each
study in `results/analysis_sweep.md`. All four lines must say identical with 0 differ and 0 missing:

```sh
grep -A1 "against \`" results/analysis_sweep.md
```

A difference here, after step 2 and step 6 passed, is a bug in the sweep code: report it with the examples the
JSON lists.

## 9. Figures

```sh
uv run -q --with matplotlib python figures/slide_charts.py --brand gel --only af_sweep af_sweep_gain
uv run -q --with matplotlib python figures/slide_charts.py --only af_sweep af_sweep_gain
```

They are written to `figures/out/`, where the slide charts live: `af_sweep_gel.png` and `.svg`,
`af_sweep_gain_gel.png` and `.svg` (and the same without `_gel` in the rewire.bio palette).

## 10. Tests

```sh
run tests/test_sweep.py
run tests/test_rule_ties.py
run tests/test_invariance.py replication
```

`test_invariance.py` reruns `replication_main.json` in place; it must leave it unchanged (`git status results/`).

## 11. What to commit

Commit `data/af5/` except the copies that `.gitignore` leaves out (`gnomad_subset/`, `pool/`, `people.tsv`,
`gene_moi.tsv`, `rep/kgp`, `rep/samples50.txt`) and the CADD and REVEL score files (not redistributable). Commit
`results/af5/build_*.md`, `results/af5_identity.md`, `results/*_sweep.json`, `results/analysis_sweep.md` and the
logs. `git status` must show no change to any v6 file under `data/` or `results/`. Then add the paragraphs to the
paper and the README, and tag `amp-2026-v7`.
