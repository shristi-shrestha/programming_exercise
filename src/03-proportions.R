# 03-proportions.R
# ---------------------------------------------------------------------------
# Exercise Q3: How do the proportions of ClusterK4_kmeans differ between Black
# and White individuals?
#
# Outputs:
#   results/03-cluster_proportions.csv
#   results/03-cluster_proportions.png  (stacked bar of within-race proportions)
# ---------------------------------------------------------------------------

source("src/00-setup.R")

expr_black <- read_expr(PATHS$expr_black)
expr_white <- read_expr(PATHS$expr_white)

# Restrict to the analysis population: ran_in_way_pipeline == TRUE with a
# non-missing cluster label (the same population used for modeling).
meta <- bind_rows(
  eligible_samples(load_metadata(PATHS$meta_black, "Black"), expr_black),
  eligible_samples(load_metadata(PATHS$meta_white, "White"), expr_white)
)

# --- Counts and within-race proportions ------------------------------------
prop_tbl <- meta %>%
  count(race, cluster, name = "n") %>%
  group_by(race) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()

# Wide view for the writeup (one row per cluster, a column per race).
wide <- prop_tbl %>%
  select(race, cluster, proportion) %>%
  pivot_wider(names_from = race, values_from = proportion, values_fill = 0)

message("== Q3: ClusterK4_kmeans proportions by race ==")
print(wide)

# --- Test for association between race and subtype -------------------------
# A chi-square test on the contingency table asks whether subtype membership
# is independent of race; a significant result means the mixes differ.
tab <- table(meta$race, meta$cluster)
chisq <- suppressWarnings(chisq.test(tab))
message(sprintf("\nChi-square test of race x cluster: X2 = %.2f, df = %d, p = %.3g",
                chisq$statistic, chisq$parameter, chisq$p.value))

# --- Plot ------------------------------------------------------------------
p <- ggplot(prop_tbl, aes(x = race, y = proportion, fill = cluster)) +
  geom_col(position = "fill", width = 0.7) +
  geom_text(aes(label = scales::percent(proportion, accuracy = 1)),
            position = position_fill(vjust = 0.5), size = 3, color = "white") +
  scale_fill_manual(values = CLUSTER_COLORS, name = "ClusterK4_kmeans") +
  scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = "Proportion of individuals",
       title = "ClusterK4_kmeans composition by race",
       caption = "03-proportions.R") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.line = element_line(color = "black", linewidth = 0.8))
save_fig(p, "03-cluster_proportions.png", width = 5, height = 5)

# --- Persist ---------------------------------------------------------------
write.csv(prop_tbl, file.path(RESULTS_DIR, "03-cluster_proportions.csv"),
          row.names = FALSE)
message("  wrote ", file.path(RESULTS_DIR, "03-cluster_proportions.csv"))
