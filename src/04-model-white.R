# 04-model-white.R
# ---------------------------------------------------------------------------
# Exercise Q4-Q5:
#   Q4. Train a supervised model of ClusterK4_kmeans from the WHITE cohort's
#       gene expression, using stratified cross-validation, and evaluate it.
#   Q5. Apply the trained model (no retraining) to the BLACK cohort and
#       evaluate.
#
# We fit two candidates (regularized multinomial + random forest), pick the
# better by cross-validated ROC AUC, and report both CV performance and the
# out-of-cohort (Black) performance.
#
# Outputs:
#   results/04-cv_metrics.csv       (CV metrics for both candidate models)
#   results/04-eval_metrics.csv     (white-CV vs. black-apply, chosen model)
#   results/04-confusion_black.png  (confusion matrix on the Black cohort)
#   results/04-model.rds            (fitted workflow, reused by nothing but kept)
# ---------------------------------------------------------------------------

source("src/00-setup.R")
source("src/model-utils.R")

# --- Build the harmonized feature set and modeling frames ------------------
expr_black <- read_expr(PATHS$expr_black)
expr_white <- read_expr(PATHS$expr_white)
globalmad  <- read_genelist(PATHS$genelist_globalmad)
genes      <- feature_genes(globalmad, list(expr_black, expr_white))

meta_white <- eligible_samples(load_metadata(PATHS$meta_white, "White"), expr_white)
meta_black <- eligible_samples(load_metadata(PATHS$meta_black, "Black"), expr_black)

white_df <- make_model_df(expr_white, meta_white, genes)
black_df <- make_model_df(expr_black, meta_black, genes)   # external test set

message(sprintf("Training on %d White samples x %d genes; external test = %d Black samples",
                nrow(white_df), length(genes), nrow(black_df)))
# Training on 309 White samples x 4355 genes; external test = 262 Black samples

# --- Stratified cross-validation -------------------------------------------
# Stratify folds by `cluster`: the subtypes are imbalanced (Q3), so stratifying
# keeps each fold's class mix representative and avoids folds missing a class.
# (Within the White-only model, race is constant, so cluster is the thing to
# stratify on here.)
set.seed(42)
folds <- vfold_cv(white_df, v = 5, strata = cluster) #defaults to v = 10

# --- Tune + fit both candidates, compare by CV ROC AUC ---------------------
# candidate_models() wraps that recipe into two full workflows each with its own hyperparameter grid
models <- candidate_models(make_recipe(white_df), n_predictors = length(genes))

# runs tune_grid() across all 5 folds × every grid combination 
# try every hyperparameter combination, on every CV fold
# picks the best hyperparameters by CV ROC AUC, finalizes that configuration on the entire White training set.
results <- lapply(names(models), function(nm) {
  message("\n--- Tuning ", nm, " ---")
  out <- tune_and_finalize(models[[nm]], folds, white_df, select_metric = "roc_auc")
  out$cv_metrics$model <- nm
  out
})
names(results) <- names(models)

# Both candidates' CV metrics are written out for comparison
cv_compare <- bind_rows(lapply(results, `[[`, "cv_metrics"))
message("\n== Q4: cross-validated performance (White cohort) ==")
print(cv_compare %>% filter(.metric %in% c("roc_auc", "bal_accuracy", "kap")))
write.csv(cv_compare, file.path(RESULTS_DIR, "04-cv_metrics.csv"), row.names = FALSE)

# Choose the model with the best CV ROC AUC as the final classifier.
best_auc <- sapply(results, function(r)
  r$cv_metrics$mean[r$cv_metrics$.metric == "roc_auc"])
chosen <- names(which.max(best_auc))
message(sprintf("\nChosen model: %s (CV ROC AUC = %.3f)", chosen, max(best_auc)))
# Chosen model: random_forest (CV ROC AUC = 0.978)

final_fit <- results[[chosen]]$fit #winning workflow refit on all white samples, result is one fitted model + its CV metrics per candidate.
saveRDS(final_fit, file.path(RESULTS_DIR, "04-model.rds")) #saved to as to not retrain later


# --- Q4: report the chosen model's CV metrics as its White performance -----
white_cv <- results[[chosen]]$cv_metrics %>%
  transmute(eval_set = "White (5-fold CV)", .metric, .estimate = mean)

# --- Q5: apply to Black cohort (no retraining) -----------------------------
# generates class predictions
# and class probabilities for all 262 Black samples, computes the class
# metrics (accuracy, balanced accuracy, kappa, macro F1) and multiclass ROC
# AUC against the true cluster labels, and builds a confusion matrix.
# every metric drops from White CV to Black apply
# the confusion matrix shows the drop is concentrated in C4 (rare in Black, common in White training) — 
# consistent with a distribution-shift / class-imbalance-transfer problem rather than a generically weaker model.
black_eval <- evaluate(final_fit, black_df, label = "Black (external apply)")

# --- Persist combined evaluation + confusion matrix ------------------------
eval_metrics <- bind_rows(
  white_cv,
  black_eval$metrics %>% select(eval_set, .metric, .estimate)
) %>% mutate(model = chosen)
write.csv(eval_metrics, file.path(RESULTS_DIR, "04-eval_metrics.csv"), row.names = FALSE)

p_cm <- autoplot(black_eval$conf_mat, type = "heatmap") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(title = sprintf("%s applied to Black cohort", chosen),
       caption = "04-model-white.R") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
save_fig(p_cm, "04-confusion_black.png", width = 5, height = 4.5)
