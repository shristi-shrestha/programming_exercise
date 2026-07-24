# Programming Exercise

A 9-question analysis exercise using gene-expression and metadata tables for
Black and White cohorts of a high-grade serous ovarian cancer (HGSC) study.
The exercise covers gene-list filtering, PCA, subtype-proportion comparisons
by race, and supervised classification of a k-means-derived molecular subtype
(`ClusterK4_kmeans`), including how well a model trained on one cohort
transfers to the other.

Full instructions: [docs/Instructions.docx](docs/Instructions.docx).
Write-up of results and answers: [docs/Report.Rmd](docs/Report.Rmd) (rendered:
[docs/Report.html](docs/Report.html)).

## Project structure

```
data/     Input data files (see Data below)
src/      Analysis scripts, numbered in run order (00-setup.R is shared/sourced by the rest)
results/  Figures and tables written by the scripts in src/
docs/     Exercise instructions and the results report
```

## Running the analysis

From the project root:

```r
Rscript src/run_all.R
```

This sources each numbered script in `src/` in order; each writes its
figures/tables to `results/`. Scripts can also be run individually (each
sources `src/00-setup.R` itself for shared paths, helpers, and constants).

## Data

The `data/` directory holds the raw input files used by the scripts (gene
lists, per-cohort metadata, and per-cohort expression matrices). Source
files: https://drive.google.com/drive/u/1/folders/11EEWzD87G6qnyJHXqRaSKFv0Jgla8Ra2
