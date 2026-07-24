# 01-explore.R
# ---------------------------------------------------------------------------
# Exercise Q1: Filter the expression tables to the GlobalMAD genes and report
#   - how many genes / samples before and after filtering, and
#   - how many genes are shared between the GlobalMAD and CommonGenes lists.
#
# Outputs: results/01-explore_counts.csv  (a small tidy summary table)
# ---------------------------------------------------------------------------

source("src/00-setup.R")

# --- Load raw inputs -------------------------------------------------------
expr_black <- read_expr(PATHS$expr_black)
expr_white <- read_expr(PATHS$expr_white)
dim(expr_black)#18641   272
dim(expr_white)# 18851   316

globalmad <- read_genelist(PATHS$genelist_globalmad)   # supp_table_1 (features)#4355
common    <- read_genelist(PATHS$genelist_common)      # supp_table_2 #8360

# --- Filter each expression table to the GlobalMAD gene list ---------------
# Not every listed gene is quantified in a given table, so we keep the genes
# that are both in the list AND present in that table (an intersection filter).
filter_to_genes <- function(expr, genes) {
  expr[intersect(rownames(expr), genes), , drop = FALSE]
}
expr_black_filt <- filter_to_genes(expr_black, globalmad)
expr_white_filt <- filter_to_genes(expr_white, globalmad)
dim(expr_black_filt) #4355  272
dim(expr_white_filt) # 4355  316

# --- Q1a: genes / samples before vs after ----------------------------------
counts <- tibble::tribble(
  ~cohort,  ~stage,             ~n_genes,               ~n_samples,
  "Black",  "before_filter",    nrow(expr_black),       ncol(expr_black),
  "Black",  "after_filter",     nrow(expr_black_filt),  ncol(expr_black_filt),
  "White",  "before_filter",    nrow(expr_white),       ncol(expr_white),
  "White",  "after_filter",     nrow(expr_white_filt),  ncol(expr_white_filt)
)

# --- Q1b: gene-list overlap ------------------------------------------------
# "Shared genes" = the set intersection of the two feature lists.
n_globalmad <- length(unique(globalmad))
n_common    <- length(unique(common))
n_shared    <- length(intersect(globalmad, common))

message("== Q1: filtering counts ==")
print(counts)
message(sprintf(
  "\n== Q1: gene-list overlap ==\nGlobalMAD genes: %d | CommonGenes: %d | shared: %d (%.1f%% of GlobalMAD)",
  n_globalmad, n_common, n_shared, 100 * n_shared / n_globalmad))
# == Q1: gene-list overlap ==
# GlobalMAD genes: 4355 | CommonGenes: 8360 | shared: 4355 (100.0% of GlobalMAD)


# Note: filtering changes the number of GENES (rows) but not SAMPLES (columns);
# the sample counts are printed unchanged before/after to make that explicit.

# --- Persist a tidy summary for the writeup --------------------------------
summary_tbl <- counts %>%
  mutate(
    n_globalmad_genes = n_globalmad,
    n_common_genes    = n_common,
    n_shared_genes    = n_shared
  )
write.csv(summary_tbl, file.path(RESULTS_DIR, "01-explore_counts.csv"),
          row.names = FALSE)
message("\n  wrote ", file.path(RESULTS_DIR, "01-explore_counts.csv"))
