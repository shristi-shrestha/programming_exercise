# 02b-pca-hvg.R
# ---------------------------------------------------------------------------
# Secondary/robustness PCA (not one of the graded exercise questions): same
# combined Black + White samples as 02-pca.R, but features are the top 500
# most-variable genes computed on this cohort instead of the curated
# GlobalMAD list. Purpose: check whether the GlobalMAD PCA structure (esp.
# the by-cluster grouping) holds up under a purely data-driven feature
# choice, and whether an HVG-driven PCA instead surfaces cohort/batch
# structure by race (see 02-pca.R header for why GlobalMAD is primary).
#
# Outputs:
#   results/02b-pca_hvg_by_cluster.png
#   results/02b-pca_hvg_by_race.png
#   results/02b-pca_hvg_scree.png
#   results/02b-pca_hvg_variance.csv
#   results/02b-pca_hvg_genelist.csv   (the 500 genes used, for reference)
# ---------------------------------------------------------------------------

source("src/00-setup.R")
suppressPackageStartupMessages(library(ggalt))  # geom_encircle for group hulls

N_HVG <- 500

# --- Assemble the combined matrix over ALL genes shared by both cohorts ----
expr_black <- read_expr(PATHS$expr_black)
expr_white <- read_expr(PATHS$expr_white)

common_genes <- Reduce(intersect, list(rownames(expr_black), rownames(expr_white)))

meta <- bind_rows(
  eligible_samples(load_metadata(PATHS$meta_black, "Black"), expr_black, require_label = FALSE),
  eligible_samples(load_metadata(PATHS$meta_white, "White"), expr_white, require_label = FALSE)
)

x_common <- rbind(
  t(expr_black[common_genes, meta$ID[meta$race == "Black"], drop = FALSE]),
  t(expr_white[common_genes, meta$ID[meta$race == "White"], drop = FALSE])
)
meta <- meta[match(rownames(x_common), meta$ID), ]   # keep row order aligned to x_common

# --- Select the top N_HVG most-variable genes across all combined samples --
# Variance is computed on the same pooled Black+White samples PCA will run
# on (not per-race), so the ranking reflects overall spread rather than
# either cohort alone -- this is also exactly the setting where high-variance
# genes are most likely to be cohort/batch effects rather than biology.
gene_var  <- apply(x_common, 2, var)
hvg_genes <- names(sort(gene_var, decreasing = TRUE))[1:N_HVG]
write.csv(tibble(gene = hvg_genes),
          file.path(RESULTS_DIR, "02b-pca_hvg_genelist.csv"), row.names = FALSE)

x <- x_common[, hvg_genes, drop = FALSE]

# --- Run PCA (same centering/scaling convention as 02-pca.R) ---------------
pca <- prcomp(x, center = TRUE, scale. = TRUE)

var_explained <- pca$sdev^2 / sum(pca$sdev^2)
scores <- as.data.frame(pca$x[, 1:3])
scores$cluster <- meta$cluster
scores$race    <- meta$race

pc_lab <- function(i) sprintf("PC%d (%.1f%%)", i, 100 * var_explained[i])

message("== Q2 (secondary check): variance explained, top ", N_HVG, " HVG ==")
for (i in 1:3) message(sprintf("  PC%d: %.1f%%", i, 100 * var_explained[i]))

# --- Plot helper: identical convention to 02-pca.R's pca_panel -------------
pca_panel <- function(df, xpc, ypc, color_var, palette) {
  ggplot(df, aes(x = .data[[xpc]], y = .data[[ypc]], color = .data[[color_var]])) +
    geom_point(size = 1.8, alpha = 0.8) +
    geom_encircle(aes(group = .data[[color_var]], fill = .data[[color_var]]),
                  alpha = 0.12, s_shape = 0.5, expand = 0.02, na.rm = TRUE) +
    scale_color_manual(values = palette, na.value = "grey70",
                       name = color_var, drop = FALSE) +
    scale_fill_manual(values = palette, na.value = "grey70", guide = "none",
                      drop = FALSE) +
    labs(x = pc_lab(as.integer(sub("PC", "", xpc))),
         y = pc_lab(as.integer(sub("PC", "", ypc)))) +
    theme_bw(base_size = 12) +
    coord_fixed() +
    theme(aspect.ratio = 1,
          axis.line = element_line(color = "black", linewidth = 0.8))
}

pc_pairs <- list(c("PC1", "PC2"), c("PC2", "PC3"), c("PC1", "PC3"))

build_panels <- function(color_var, palette, title) {
  panels <- lapply(pc_pairs, function(p) pca_panel(scores, p[1], p[2], color_var, palette))
  wrap_plots(panels, nrow = 1, guides = "collect") +
    plot_annotation(title = title, caption = "02b-pca-hvg.R") &
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5, face = "bold"))
}

p_cluster <- build_panels("cluster", CLUSTER_COLORS,
                          sprintf("Combined PCA (top %d HVG features), colored by ClusterK4_kmeans", N_HVG))
p_race    <- build_panels("race", RACE_COLORS,
                          sprintf("Combined PCA (top %d HVG features), colored by race", N_HVG))

save_fig(p_cluster, "02b-pca_hvg_by_cluster.png", width = 11, height = 4.5)
save_fig(p_race,    "02b-pca_hvg_by_race.png",    width = 11, height = 4.5)

# --- Scree plot (variance explained, first 10 PCs) -------------------------
scree_df <- tibble(PC = factor(paste0("PC", 1:10), levels = paste0("PC", 1:10)),
                   variance = 100 * var_explained[1:10])
p_scree <- ggplot(scree_df, aes(PC, variance)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = sprintf("%.1f", variance)), vjust = -0.4, size = 3) +
  labs(x = NULL, y = "% variance explained",
       title = sprintf("Scree plot (top %d HVG features, top 10 PCs)", N_HVG),
       caption = "02b-pca-hvg.R") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
save_fig(p_scree, "02b-pca_hvg_scree.png", width = 6, height = 4)

# --- Persist variance table -------------------------------------------------
write.csv(
  tibble(PC = paste0("PC", 1:10),
         prop_variance = var_explained[1:10],
         pct_variance  = 100 * var_explained[1:10]),
  file.path(RESULTS_DIR, "02b-pca_hvg_variance.csv"), row.names = FALSE)
message("  wrote ", file.path(RESULTS_DIR, "02b-pca_hvg_variance.csv"))
