# 01b-expr-qc.R
# ---------------------------------------------------------------------------
# QC diagnostic (not a graded exercise question): are per-sample expression
# distributions on comparable scales between the Black and White cohorts
# before they're combined for 02-pca.R? This is the same concern the source
# paper's methods address with per-sample scaling and a log10(x+1) transform
# before pooling cohorts/platforms into one PCA -- here we only have two
# cohorts off (presumably) the same RNA-seq pipeline, but it's worth checking
# rather than assuming.
#
# Uses the same eligible-sample, GlobalMAD-filtered matrix that feeds
# 02-pca.R. Raw values are heavily right-skewed (max ~1e6, median ~1e3), so
# we log10(x+1) transform before plotting, same as the source paper.
#
# Outputs:
#   results/01b-expr_qc_density_box.png
# ---------------------------------------------------------------------------

source("src/00-setup.R")

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

x_log <- log10(x + 1)

long_df <- as.data.frame(x_log, check.names = FALSE) %>%
  mutate(ID = rownames(x_log), race = meta$race) %>%
  pivot_longer(cols = -c(ID, race), names_to = "gene", values_to = "log10_expr")

# --- Density: one line per sample, colored by cohort ------------------------
p_density <- ggplot(long_df, aes(log10_expr, group = ID, color = race)) +
  geom_density(alpha = 0.05, linewidth = 0.2) +
  scale_color_manual(values = RACE_COLORS, name = "race") +
  labs(x = "log10(expression + 1)", y = "Density",
      title = "Per-sample expression density (GlobalMAD genes)") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# --- Boxplot: one box per sample, ordered by cohort then median ------------
# No x-axis tick labels: with 500+ samples individual sample IDs would be
# unreadable regardless of angle, so we drop them and keep the fill legend
# as the only cohort indicator.
sample_order <- long_df %>%
  group_by(ID, race) %>%
  summarize(med = median(log10_expr), .groups = "drop") %>%
  arrange(race, med) %>%
  pull(ID)
long_df$ID <- factor(long_df$ID, levels = sample_order)

p_box <- ggplot(long_df, aes(ID, log10_expr, fill = race)) +
  geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.3, linewidth = 0.15) +
  scale_fill_manual(values = RACE_COLORS, name = "race") +
  labs(x = sprintf("Sample (n = %d, ordered by cohort then median)", length(sample_order)),
      y = "log10(expression + 1)",
      title = "Per-sample expression boxplots (GlobalMAD genes)") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
       plot.title = element_text(hjust = 0.5, face = "bold"))

p <- (p_density / p_box) + plot_annotation(caption = "01b-expr-qc.R")

save_fig(p, "01b-expr_qc_density_box.png", width = 9, height = 9)
