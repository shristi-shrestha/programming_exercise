# 02c-pca-metadata-corr.R
# ---------------------------------------------------------------------------
# Diagnostic (not a graded exercise question): for the combined GlobalMAD PCA
# from 02-pca.R, test PC1-PC5 against every metadata column available for
# Black and White samples -- not just race/ClusterK4_kmeans -- to see what
# else associates with the PCA structure. This is the kind of check flagged
# under "other data cleaning/validation steps" in the instructions: technical
# covariates (sequencing `version`, `resequenced`, `low_qual` flags) leaking
# into early PCs would point to a batch effect rather than biology.
#
# Association metric: for each (PC, covariate) pair, fit
# lm(PC ~ factor(covariate)) and report R^2 (proportion of that PC's
# variance explained by the covariate) and the model F-test p-value. This is
# the same logic as PVCA -- categorical covariates, so ANOVA R^2 rather than
# a numeric correlation coefficient.
#
# `REMOVE_*` / `resequenced` / `low_qual` / `failed_seq` are flag columns
# where NA means "not flagged", not "unknown" -- recoded to an explicit
# FALSE/"not_flagged" level before testing. All other categorical columns
# (cluster labels, version, race) keep NA as genuine missingness and those
# rows are simply dropped for that test (lm's default).
#
# Outputs:
#   results/02c-pca_metadata_corr.png   (covariate x PC heatmap of R^2)
#   results/02c-pca_metadata_corr.csv   (long-format R^2 / p-value / n table)
# ---------------------------------------------------------------------------

source("src/00-setup.R")
suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)   # colorRamp2
  library(ggplotify)  # as.ggplot
})

N_PC <- 5

# --- Rebuild the combined GlobalMAD PCA exactly as in 02-pca.R -------------
expr_black <- read_expr(PATHS$expr_black)
expr_white <- read_expr(PATHS$expr_white)
globalmad  <- read_genelist(PATHS$genelist_globalmad)
genes      <- feature_genes(globalmad, list(expr_black, expr_white))

meta <- bind_rows(
  eligible_samples(load_metadata(PATHS$meta_black, "Black"), expr_black, require_label = FALSE),
  eligible_samples(load_metadata(PATHS$meta_white, "White"), expr_white, require_label = FALSE)
)

x <- rbind(
  t(expr_black[genes, meta$ID[meta$race == "Black"], drop = FALSE]),
  t(expr_white[genes, meta$ID[meta$race == "White"], drop = FALSE])
)
meta <- meta[match(rownames(x), meta$ID), ]

pca    <- prcomp(x, center = TRUE, scale. = TRUE)
scores <- as.data.frame(pca$x[, 1:N_PC])
colnames(scores) <- paste0("PC", 1:N_PC)

# --- Recode flag columns: NA -> explicit "not flagged" level ---------------
flag_cols <- grep("^REMOVE_", names(meta), value = TRUE)
flag_cols <- union(flag_cols, intersect(c("resequenced", "low_qual", "failed_seq"), names(meta)))

recode_flag <- function(v) {
  if (is.logical(v)) return(factor(ifelse(is.na(v), "FALSE", "TRUE")))
  lvl <- unique(v[!is.na(v)])
  if (length(lvl) == 0) return(factor(rep(NA, length(v))))  # no information at all
  factor(ifelse(is.na(v), paste0("not_", lvl[1]), v))
}
for (col in flag_cols) meta[[col]] <- recode_flag(meta[[col]])

# --- Covariates to test: every metadata column except sample identifiers ---
# and `cluster` (a re-labeled duplicate of ClusterK4_kmeans already tested).
exclude_cols <- c("ID", "Sample_ID", "cluster")
covariate_cols <- setdiff(names(meta), exclude_cols)

# --- lm(PC ~ factor(covariate)) -> R^2 and p-value for one (PC, covariate) -
pc_covariate_r2 <- function(pc_values, covariate) {
  d <- data.frame(pc = pc_values, cov = covariate)
  d <- d[!is.na(d$cov), , drop = FALSE]
  if (nlevels(droplevels(factor(d$cov))) < 2 || nrow(d) < 3) {
    return(c(r2 = NA_real_, p = NA_real_, n = nrow(d)))
  }
  fit <- lm(pc ~ factor(cov), data = d)
  c(r2 = summary(fit)$r.squared,
    p  = anova(fit)[1, "Pr(>F)"],
    n  = nrow(d))
}

results <- expand.grid(covariate = covariate_cols, PC = colnames(scores),
                       stringsAsFactors = FALSE)
stats <- t(mapply(function(cov, pc) pc_covariate_r2(scores[[pc]], meta[[cov]]),
                  results$covariate, results$PC))
results <- cbind(results, as.data.frame(stats))
write.csv(results, file.path(RESULTS_DIR, "02c-pca_metadata_corr.csv"), row.names = FALSE)
message("  wrote ", file.path(RESULTS_DIR, "02c-pca_metadata_corr.csv"))

# --- Heatmap: covariate (rows) x PC (cols), fill = R^2 ----------------------
mat_df <- tidyr::pivot_wider(results[, c("covariate", "PC", "r2")],
                             names_from = PC, values_from = r2)
mat <- as.matrix(mat_df[, colnames(scores), drop = FALSE])
rownames(mat) <- mat_df$covariate

pval_lookup <- results$p
names(pval_lookup) <- paste(results$covariate, results$PC)
stars <- function(p) if (is.na(p)) "" else if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else ""

col_fun <- colorRamp2(seq(0, max(mat, na.rm = TRUE), length.out = 9),
                      rev(RColorBrewer::brewer.pal(9, "RdBu")))

h <- Heatmap(mat, name = "R2",
            col = col_fun, na_col = "grey85",
            cluster_columns = FALSE, cluster_rows = FALSE,
            row_names_side = "left", column_names_rot = 0,
            column_title = sprintf("PC1-%d vs. all metadata (Black + White combined)", N_PC),
            cell_fun = function(j, i, x, y, width, height, fill) {
              v <- mat[i, j]
              if (is.na(v)) return(invisible())
              lab <- sprintf("%.2f%s", v, stars(pval_lookup[paste(rownames(mat)[i], colnames(mat)[j])]))
              grid::grid.text(lab, x, y, gp = grid::gpar(fontsize = 9,
                                                         col = if (v > 0.5) "white" else "black"))
            })
p <- as.ggplot(grid::grid.grabExpr(draw(h, padding = unit(c(2, 20, 2, 2), "mm")))) +
  labs(caption = "02c-pca-metadata-corr.R")

ggsave(file.path(RESULTS_DIR, "02c-pca_metadata_corr.png"), plot = p,
      width = 7, height = 8, units = "in", scale = 1.5, dpi = 300)
message("  wrote ", file.path(RESULTS_DIR, "02c-pca_metadata_corr.png"))
