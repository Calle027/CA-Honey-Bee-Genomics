# California Honey Bee Genomics

Analysis scripts for population genomic analyses of California honey bee samples.

This repository is organized so that shared metadata, sample filtering, LD-pruned SNP sets, and KING relatedness objects are created once, then reused by focused figure/analysis scripts.

## Repository Structure

```text
honeybee_scripts/
  00_config.R                     # Shared paths, metadata, colors, helpers, sample map
  00_create_essential_objects.R   # LD pruning, pruned SNP table, KING/IBS0 objects
  01_inbreeding.R                 # Inbreeding statistics and figure
  02_pcair_mclust.R               # PC-AiR, related/unrelated status, Mclust figure
  03_relatedness_ibs0.R           # IBS0 heatmap and relatedness panel figures
  04_tree.R                       # Rooted ML tree figure
  05_nucleotide_diversity.R       # Windowed nucleotide diversity summaries and figure
  06_phylogeny_admixture_mclust.R # Tree with K=7 ADMIXTURE bars and Mclust tiles
  07_supplemental_mclust_figures.R # Supplemental Mclust-colored PC-AiR and support trees

analysis_outputs/                 # Generated tables and cached R objects
figures/                          # Generated figures
```

The original exploratory script, `genomics.R`, is retained as a reference while the modular scripts become the reproducible workflow.

## Required Inputs

The scripts currently expect these files in the project root:

- `final_cleaned.gds`
- `HB_fam.xlsx`, with a `raw_data` sheet containing at least `ID` and `pop`
- `gtr_ascR4_nonsyn_constrained_varsites.contree`, for the tree figure
- `admixture_input.7.Q`, for the ADMIXTURE-annotated tree
- `*_maf05.windowed.pi` files, for nucleotide diversity
- `genome_coverage.tsv`, for the manuscript sample metadata table

Large genomic inputs and generated outputs should generally not be committed to GitHub unless there is a specific reason to publish them.

## Running The Analysis In RStudio

Open this folder in RStudio, then open and click **Source** on the scripts below in order:

1. `honeybee_scripts/00_create_essential_objects.R`
2. `honeybee_scripts/01_inbreeding.R`
3. `honeybee_scripts/02_pcair_mclust.R`
4. `honeybee_scripts/03_relatedness_ibs0.R`
5. `honeybee_scripts/04_tree.R`
6. `honeybee_scripts/05_nucleotide_diversity.R`
7. `honeybee_scripts/06_phylogeny_admixture_mclust.R`
8. `honeybee_scripts/07_supplemental_mclust_figures.R`
9. `honeybee_scripts/08_manuscript_metadata.R`

Run `00_create_essential_objects.R` first. It creates the LD-pruned SNP set and KING objects used by the downstream scripts.
Run `02_pcair_mclust.R` before `06_phylogeny_admixture_mclust.R`; the annotated tree needs `analysis_outputs/pcair_scores_with_mclust.csv`.
Run `08_manuscript_metadata.R` after `01_inbreeding.R` and `02_pcair_mclust.R`; it creates `manuscript_sample_metadata.csv`.

The scripts detect the project folder from the active RStudio document, so they do not require terminal arguments. If needed, you can still override the project path inside RStudio before sourcing a script:

```r
Sys.setenv(CA_HONEYBEE_ROOT = "/path/to/CA Honey Bee Genomics")
```

## Outputs

Generated tables and cached objects are written to `analysis_outputs/`. Generated figures are written to `figures/`.

Key cached files include:

- `analysis_outputs/snpset_ld_pruned_ids.rds`
- `analysis_outputs/king_ibd_mom.rds`
- `analysis_outputs/king_matrix.rds`
- `analysis_outputs/sample_map.csv`
- `manuscript_sample_metadata.csv`

## R Packages

The workflow uses packages including:

`gdsfmt`, `SeqArray`, `SNPRelate`, `GENESIS`, `GWASTools`, `readxl`, `readr`, `dplyr`, `tibble`, `tidyr`, `ggplot2`, `ragg`, `ComplexHeatmap`, `circlize`, `RColorBrewer`, `patchwork`, `ggstar`, `ggnewscale`, `colorspace`, `mclust`, `ape`, `treeio`, `ggtree`, `FSA`, and `rcompanion`.

## Notes

This is an initial cleanup pass. The next useful step is to add a dependency file, such as `renv.lock`, once the package versions used for the manuscript figures are finalized.
