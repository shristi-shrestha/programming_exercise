# 00-setup.R
# ---------------------------------------------------------------------------
# Shared libraries, paths, and helper functions for the programming exercise.
# Every numbered script sources this file first, so all the data-loading and
# harmonization logic lives here in one place (single source of truth).
#
# Assumed working directory: the `programming_exercise/` folder.
# Run scripts as, e.g.:  Rscript src/01-explore.R   (from that folder)
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})

set.seed(42)  # reproducibility for anything stochastic downstream

# --- Paths -----------------------------------------------------------------
DATA_DIR    <- "data"
RESULTS_DIR <- "results"
if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)

# The four data files we use repeatedly.
PATHS <- list(
  genelist_globalmad = file.path(DATA_DIR, "supp_table_1_GlobalMAD_genelist.csv"),
  genelist_common    = file.path(DATA_DIR, "supp_table_2_CommonGenes_genelist.csv"),
  meta_black         = file.path(DATA_DIR, "supp_table_3_main_black_metadata_table.tsv"),
  meta_white         = file.path(DATA_DIR, "supp_table_4_main_white_metadata_table.tsv"),
  expr_black         = file.path(DATA_DIR, "supp_table_6_black_expr.tsv"),
  expr_white         = file.path(DATA_DIR, "supp_table_7_white_expr.tsv")
)

# --- Constants -------------------------------------------------------------
# ClusterK4_kmeans is the prediction target. It is stored as integers 1-4;
# we relabel to syntactically-valid factor levels ("C1".."C4") because some
# model engines (randomForest, yardstick class metrics) prefer valid names.
CLUSTER_COL    <- "ClusterK4_kmeans"
CLUSTER_LEVELS <- c("C1", "C2", "C3", "C4")

# A consistent color scheme reused across every plot so a subtype/race always
# maps to the same color (house convention).
CLUSTER_COLORS <- setNames(RColorBrewer::brewer.pal(4, "Set1"), CLUSTER_LEVELS)
RACE_COLORS    <- setNames(RColorBrewer::brewer.pal(3, "Dark2")[1:2], c("Black", "White"))

# external_HGSCsubtype_estimate: a second, TCGA-signature-based subtype call
# (independent of this exercise's own ClusterK4_kmeans), only defined for a
# subset of samples. Own palette -- its 4 levels are not guaranteed to line up
# 1:1 with the ClusterK4_kmeans integer IDs, so we don't reuse CLUSTER_COLORS.
TCGA_SUBTYPE_LEVELS <- c("C1.MES", "C2.IMM", "C4.DIF", "C5.PRO")
TCGA_SUBTYPE_COLORS <- setNames(RColorBrewer::brewer.pal(4, "Set1"), TCGA_SUBTYPE_LEVELS)

# ---------------------------------------------------------------------------
# read_expr(): read one expression table.
#
# The files were written by R with row names but no label for the first
# (gene) column, so the header row has one fewer field than the data rows.
# read.table(header = TRUE, row.names = 1) resolves this off-by-one for us.
# Expression column names carry a "Sample_" prefix that the metadata `ID`
# column lacks, so we strip it here to make the two joinable.
# Returns a numeric matrix: genes (rows) x samples (columns).
# ---------------------------------------------------------------------------
read_expr <- function(path) {
  m <- read.table(path, header = TRUE, row.names = 1, sep = "\t",
                  check.names = FALSE, quote = "\"")
  colnames(m) <- sub("^Sample_", "", colnames(m))
  as.matrix(m)
}

# ---------------------------------------------------------------------------
# read_genelist(): read a one-column gene list CSV (column header is "x").
# Returns a character vector of gene symbols.
# ---------------------------------------------------------------------------
read_genelist <- function(path) {
  read.csv(path, stringsAsFactors = FALSE)[[1]]
}

# ---------------------------------------------------------------------------
# load_metadata(): read a metadata table, tag it with `race`, and tidy the
# target column. `race` is NOT a column in the data -- it is defined by which
# table a sample came from -- so we attach it explicitly here.
# ---------------------------------------------------------------------------
load_metadata <- function(path, race) {
  read.table(path, header = TRUE, sep = "\t", check.names = FALSE,
             quote = "\"", stringsAsFactors = FALSE) %>%
    mutate(
      race = race,
      # Relabel integer clusters 1-4 -> factor C1-C4 (NA stays NA).
      cluster = factor(paste0("C", .data[[CLUSTER_COL]]),
                       levels = CLUSTER_LEVELS)
    )
}

# ---------------------------------------------------------------------------
# eligible_samples(): the samples we are told to use downstream are those with
# ran_in_way_pipeline == TRUE. For modeling we additionally need a non-missing
# cluster label. `require_label = FALSE` keeps unlabeled samples (used by PCA,
# which does not need the label to compute components).
# Returns the metadata rows whose IDs are also present in the expression matrix.
# ---------------------------------------------------------------------------
eligible_samples <- function(meta, expr, require_label = TRUE) {
  keep <- meta$ran_in_way_pipeline %in% TRUE & meta$ID %in% colnames(expr)
  if (require_label) keep <- keep & !is.na(meta$cluster)
  meta[keep, , drop = FALSE]
}

# ---------------------------------------------------------------------------
# feature_genes(): the feature set shared across cohorts. Because the Black and
# White expression tables were quantified against slightly different gene sets,
# and we want Black/White samples to live in the SAME feature space (needed to
# put them in one PCA and to apply a White-trained model to Black samples), we
# intersect the requested gene list with the genes present in BOTH tables.
# ---------------------------------------------------------------------------
feature_genes <- function(genelist, expr_list) {
  common <- Reduce(intersect, lapply(expr_list, rownames))
  intersect(genelist, common)
}

# ---------------------------------------------------------------------------
# make_model_df(): assemble a tidy modeling data frame (one row per sample):
#   columns = feature genes (numeric) + `cluster` (factor outcome) + `race`.
# `expr` is genes x samples; we subset to `genes`, subset/align to the given
# metadata rows, transpose to samples x genes, and bind the labels.
# ---------------------------------------------------------------------------
make_model_df <- function(expr, meta, genes) {
  genes <- intersect(genes, rownames(expr))
  x <- t(expr[genes, meta$ID, drop = FALSE])          # samples x genes
  df <- as.data.frame(x, check.names = FALSE)
  df$cluster <- meta$cluster #outcome
  df$race    <- meta$race #for book keeping , not as predictor
  rownames(df) <- meta$ID
  df
}

# ---------------------------------------------------------------------------
# save_fig(): thin ggsave wrapper enforcing the house defaults
# (inches, dpi 300, scale multiplier) and echoing the path.
# ---------------------------------------------------------------------------
save_fig <- function(plot, filename, width, height, scale = 1.5, dpi = 300) {
  path <- file.path(RESULTS_DIR, filename)
  ggsave(path, plot = plot, width = width, height = height,
         units = "in", scale = scale, dpi = dpi)
  message("  wrote ", path)
  invisible(path)
}
