# gnomAD overlap of the first run's 200 SNV cases, with the new frequencies

Made by `gnomad_af.py overlap`. Issue #1's table, redone. Bins use the overall gnomAD v3.1.2 genome AF; "not in gnomAD" means AF 0 (not found, or no allele left after leaving people out). The spiked-in (diagnostic) variants are not carried by the person, so they keep the full frequency. Full-frequency rows count distinct variants, as in the issue; leave-out rows count person variants (the same variant in two people can have different frequencies).

| Variants | Frequency | n | Not in gnomAD | Below 0.01% | 0.01 to 0.1% | 0.1% or more |
|---|---|---|---|---|---|---|
| 200 cases, spiked-in | full | 200 | 76% | 22% | 2% | 0% |
| 200 cases, background | full | 13,145 | 0% | 16% | 26% | 57% |
| 200 cases, background, per person variant | full | 15,216 | 0% | 14% | 24% | 61% |
| 200 cases, background | family left out | 15,216 | 2% | 12% | 24% | 61% |
| 200 cases, background | all 1000 Genomes left out | 15,216 | 4% | 12% | 24% | 60% |

## Background variants per case left after the grpmax filter

Median (and range) over the 200 cases of the background SNVs kept, by the frequency used for grpmax.

| Frequency | Before the filter | grpmax below 0.1% | grpmax below 1% |
|---|---|---|---|
| full | 233.5 (172 to 604) | 49 (13 to 96) | 139 (82 to 287) |
| family left out | 233.5 (172 to 604) | 50.5 (15 to 97) | 139.5 (82 to 287) |
| all 1000 Genomes left out | 233.5 (172 to 604) | 50 (16 to 97) | 140 (81 to 291) |

All 15,216 background person variants of the 200 cases join to `data/main/person_variant_af.tsv`.
