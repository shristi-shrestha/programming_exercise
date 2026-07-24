# 05-model-combined.R
# ---------------------------------------------------------------------------
# Exercise Q6: Train on ALL individuals (Black + White) but hold out 20% for a
# final evaluation. Use a stratified split for both cross-validation and the
# final held-out set. Evaluate on the held-out set and explain why this model
# does better/worse than the White-only model from Q4-5.
#
# Key choice: we stratify by the COMBINATION of race and subtype. Q3 shows the
# subtype mix differs by race, and the PCA (Q2) shows race is a major axis of
# variation -- so a split that is stratified on subtype alone could still end
# up race-imbalanced. Stratifying on race x cluster keeps both balanced in the
# train / test / CV folds.
#
# Outputs:
#   results/05-cv_metrics.csv
#   results/05-heldout_metrics.csv
#   results/05-confusion_heldout.png
# ---------------------------------------------------------------------------

source("src/00-setup.R")
source("src/model-utils.R")

# Most ofthe mechanics (recipe, candidate models, tune_grid, CV-ROC-AUC selection)
# are identical to 04-model-white.R

# --- Combined harmonized modeling frame ------------------------------------
expr_black <- read_expr(PATHS$expr_black)
expr_white <- read_expr(PATHS$expr_white)
globalmad  <- read_genelist(PATHS$genelist_globalmad)
genes      <- feature_genes(globalmad, list(expr_black, expr_white))

meta_white <- eligible_samples(load_metadata(PATHS$meta_white, "White"), expr_white)
meta_black <- eligible_samples(load_metadata(PATHS$meta_black, "Black"), expr_black)

# 309 White rows + 262 Black rows = 571 total, all in the same 4,355-gene feature space so they can be modeled together.
all_df <- bind_rows(
  make_model_df(expr_white, meta_white, genes),
  make_model_df(expr_black, meta_black, genes)
)
# Combined stratification key: subtype within race (e.g. "White_C1").
# stratify on race and subtype jointly, subtype alone could still leave a split race-imbalanced
all_df$strata <- factor(paste(all_df$race, all_df$cluster, sep = "_"))

message(sprintf("Combined modeling frame: %d samples (%d White, %d Black) x %d genes",
                nrow(all_df), sum(all_df$race == "White"),
                sum(all_df$race == "Black"), length(genes)))

# --- Stratified 80/20 split + stratified CV on the training 80% ------------
set.seed(42)
split    <- initial_split(all_df, prop = 0.8, strata = strata)
train_df <- training(split)
test_df  <- testing(split)
folds    <- vfold_cv(train_df, v = 5, strata = strata)

# `strata` was only for splitting; drop it before modeling so it is not a
# predictor (race is already given an ID role inside the recipe).
train_df$strata <- NULL
test_df$strata  <- NULL

# --- Tune both candidates, choose by CV ROC AUC ----------------------------
models <- candidate_models(make_recipe(train_df), n_predictors = length(genes))

results <- lapply(names(models), function(nm) {
  message("\n--- Tuning ", nm, " ---")
  out <- tune_and_finalize(models[[nm]], folds, train_df, select_metric = "roc_auc")
  out$cv_metrics$model <- nm
  out
})
names(results) <- names(models)

cv_compare <- bind_rows(lapply(results, `[[`, "cv_metrics"))
message("\n== Q6: cross-validated performance (combined training set) ==")
print(cv_compare %>% filter(.metric %in% c("roc_auc", "bal_accuracy", "kap")))
write.csv(cv_compare, file.path(RESULTS_DIR, "05-cv_metrics.csv"), row.names = FALSE)

best_auc <- sapply(results, function(r)
  r$cv_metrics$mean[r$cv_metrics$.metric == "roc_auc"])
chosen <- names(which.max(best_auc))
message(sprintf("\nChosen model: %s (CV ROC AUC = %.3f)", chosen, max(best_auc)))
final_fit <- results[[chosen]]$fit

# --- Final evaluation on the 20% held-out set ------------------------------
heldout_eval <- evaluate(final_fit, test_df, label = "Combined 20% held-out")

# Also break the held-out performance down by race, to see whether pooling the
# cohorts actually fixed the cross-cohort gap seen in Q5.
by_race <- lapply(split(seq_len(nrow(test_df)), test_df$race), function(idx) {
  ev <- evaluate(final_fit, test_df[idx, , drop = FALSE],
                 label = paste0("Held-out: ", test_df$race[idx][1]))
  ev$metrics
})
heldout_metrics <- bind_rows(heldout_eval$metrics, bind_rows(by_race)) %>%
  mutate(model = chosen)
write.csv(heldout_metrics, file.path(RESULTS_DIR, "05-heldout_metrics.csv"),
          row.names = FALSE)

p_cm <- autoplot(heldout_eval$conf_mat, type = "heatmap") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(title = sprintf("%s on combined 20%% held-out set", chosen),
       caption = "05-model-combined.R") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
save_fig(p_cm, "05-confusion_heldout.png", width = 5, height = 4.5)
