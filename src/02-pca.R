# 02-pca.R
# ---------------------------------------------------------------------------
# Exercise Q2: PCA of Black + White expression together on the GlobalMAD
# feature set. Plot PC1v2, PC2v3, PC1v3, colored by (a) ClusterK4_kmeans and
# (b) race; report variance explained by PC1-3. Also colored by (c)
# external_HGSCsubtype_estimate as a secondary check (not requested by the
# instructions, but a second, independently-derived subtype call worth
# comparing against ClusterK4_kmeans).
#
# Outputs:
#   results/02-pca_by_cluster.png      (3 PC pairs, colored by subtype)
#   results/02-pca_by_race.png         (3 PC pairs, colored by race)
#   results/02-pca_by_tcga_subtype.png (3 PC pairs, colored by external_HGSCsubtype_estimate)
#   results/02-pca_scree.png           (variance explained)
#   results/02-pca_variance.csv        (proportion of variance per PC)
# ---------------------------------------------------------------------------

source("src/00-setup.R")
suppressPackageStartupMessages(library(ggalt))  # geom_encircle for group hulls

# --- Assemble the combined, harmonized matrix ------------------------------
expr_black <- read_expr(PATHS$expr_black)
expr_white <- read_expr(PATHS$expr_white)
globalmad  <- read_genelist(PATHS$genelist_globalmad)

# Feature space shared by both cohorts (see feature_genes() rationale).
genes <- feature_genes(globalmad, list(expr_black, expr_white))
length(genes) #4355

# Use the downstream-eligible samples (ran_in_way_pipeline == TRUE). We do NOT
# require a cluster label here -- PCA is unsupervised -- but unlabeled points
# simply show up grey in the by-cluster plot.
meta <- bind_rows(
  eligible_samples(load_metadata(PATHS$meta_black, "Black"), expr_black, require_label = FALSE),
  eligible_samples(load_metadata(PATHS$meta_white, "White"), expr_white, require_label = FALSE)
)

# samples x genes matrix, columns ordered to match `meta`.
x <- rbind(
  t(expr_black[genes, meta$ID[meta$race == "Black"], drop = FALSE]),
  t(expr_white[genes, meta$ID[meta$race == "White"], drop = FALSE])
)
meta <- meta[match(rownames(x), meta$ID), ]   # keep row order aligned to x

# --- Run PCA ---------------------------------------------------------------
# center = TRUE, scale. = TRUE puts every gene on equal footing (unit variance)
# so that a handful of very highly-expressed genes don't dominate the axes --
# the standard choice for PCA on expression features.
pca <- prcomp(x, center = TRUE, scale. = TRUE)

var_explained <- pca$sdev^2 / sum(pca$sdev^2)
scores <- as.data.frame(pca$x[, 1:3])
scores$cluster      <- meta$cluster
scores$race         <- meta$race
scores$tcga_subtype <- factor(meta$external_HGSCsubtype_estimate, levels = TCGA_SUBTYPE_LEVELS)

pc_lab <- function(i) sprintf("PC%d (%.1f%%)", i, 100 * var_explained[i])

message("== Q2: variance explained ==")
for (i in 1:3) message(sprintf("  PC%d: %.1f%%", i, 100 * var_explained[i]))

# --- Plot helper: one PC-pair scatter, colored by a chosen grouping --------
# No per-sample text labels: with hundreds of samples they would be unreadable,
# so we lean on color + translucent convex hulls (geom_encircle) instead.
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

build_panels <- function(color_var, palette, title, df = scores) {
  panels <- lapply(pc_pairs, function(p) pca_panel(df, p[1], p[2], color_var, palette))
  wrap_plots(panels, nrow = 1, guides = "collect") +
    plot_annotation(title = title, caption = "02-pca.R") &
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5, face = "bold"))
}

p_cluster <- build_panels("cluster", CLUSTER_COLORS,
                          "Combined PCA (GlobalMAD features), colored by ClusterK4_kmeans")
p_race    <- build_panels("race", RACE_COLORS,
                          "Combined PCA (GlobalMAD features), colored by race")

# external_HGSCsubtype_estimate is only defined for a subset of samples (a
# separate TCGA-signature-based subtype call); the rest are dropped here
# rather than plotted grey, so the panel only shows the labeled subset.
# p_tcga    <- build_panels("tcga_subtype", TCGA_SUBTYPE_COLORS,
#                           "Combined PCA (GlobalMAD features), colored by external_HGSCsubtype_estimate",
#                           df = scores[!is.na(scores$tcga_subtype), ])

save_fig(p_cluster, "02-pca_by_cluster.png", width = 11, height = 4.5)
save_fig(p_race,    "02-pca_by_race.png",    width = 11, height = 4.5)
# save_fig(p_tcga,    "02-pca_by_tcga_subtype.png", width = 11, height = 4.5)

# --- Scree plot (variance explained, first 10 PCs) -------------------------
scree_df <- tibble(PC = factor(paste0("PC", 1:10), levels = paste0("PC", 1:10)),
                   variance = 100 * var_explained[1:10])
p_scree <- ggplot(scree_df, aes(PC, variance)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = sprintf("%.1f", variance)), vjust = -0.4, size = 3) +
  labs(x = NULL, y = "% variance explained",
       title = "Scree plot (top 10 PCs)", caption = "02-pca.R") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
save_fig(p_scree, "02-pca_scree.png", width = 6, height = 4)

# --- Persist variance table ------------------------------------------------
write.csv(
  tibble(PC = paste0("PC", 1:10),
         prop_variance = var_explained[1:10],
         pct_variance  = 100 * var_explained[1:10]),
  file.path(RESULTS_DIR, "02-pca_variance.csv"), row.names = FALSE)
message("  wrote ", file.path(RESULTS_DIR, "02-pca_variance.csv"))
