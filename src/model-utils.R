# model-utils.R
# ---------------------------------------------------------------------------
# Shared tidymodels helpers for the two supervised tasks (Q4-5 and Q6):
#   * two candidate classifiers (regularized multinomial vs. random forest),
#   * a tune-and-select routine using stratified cross-validation, and
#   * a single evaluation routine (metrics + confusion matrix).
# Sourced by 04-model-white.R and 05-model-combined.R.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidymodels)   # parsnip, recipes, rsample, tune, yardstick, workflows
  library(glmnet)
  library(randomForest)
})
tidymodels_prefer()

# Metrics reported everywhere. accuracy is intuitive but misleading under class
# imbalance (Q3 shows the subtypes are imbalanced), so we lead with metrics that
# are imbalance-aware: balanced accuracy, Cohen's kappa, macro-averaged F1, and
# multiclass ROC AUC (Hand-Till). mn_log_loss rewards calibrated probabilities.
CLASS_METRICS <- metric_set(accuracy, bal_accuracy, kap, f_meas)

# ---------------------------------------------------------------------------
# make_recipe(): outcome ~ all gene features. `race` is carried in the frame
# for stratification/inspection but must not be a predictor, so we tag it with
# an "ID" role. step_zv drops zero-variance genes; step_normalize centers and
# scales -- essential for glmnet's penalty to treat genes comparably (and
# harmless for random forest).
# ---------------------------------------------------------------------------
make_recipe <- function(train_df) {
  recipe(cluster ~ ., data = train_df) %>%
    update_role(race, new_role = "ID") %>% #race given an ID role so it’s never a predictor.
    step_zv(all_predictors()) %>% # drop zero-variance genes
    step_normalize(all_predictors())
}

# ---------------------------------------------------------------------------
# candidate_models(): the two workflows to compare, each paired with a tuning
# grid. Returned as a named list so the caller can loop over them.
#
#   elastic_net : multinomial logistic regression with an L1/L2 penalty
#                 (glmnet). Well suited to p >> n expression data -- built-in
#                 regularization + feature selection, natively multiclass,
#                 and the coefficients stay interpretable.
#   random_forest : non-linear tree ensemble; captures interactions, needs no
#                 scaling, and gives a strong non-linear baseline to compare
#                 the linear model against.
# ---------------------------------------------------------------------------
candidate_models <- function(rec, n_predictors) {
  en_spec <- multinom_reg(penalty = tune(), mixture = tune()) %>%
    set_engine("glmnet")
  en_grid <- grid_regular(
    penalty(range = c(-4, 0)),          # 1e-4 .. 1 on the log10 scale
    mixture(range = c(0, 1)),           # 0 = ridge .. 1 = lasso
    levels = c(penalty = 10, mixture = 5))

  rf_spec <- rand_forest(mtry = tune(), min_n = tune(), trees = 500) %>%
    set_engine("randomForest") %>%
    set_mode("classification")
  rf_grid <- grid_regular(
    mtry(range = c(20, min(200, n_predictors))),
    min_n(range = c(2, 10)),
    levels = c(mtry = 5, min_n = 3))

  list(
    elastic_net   = list(wf = workflow() %>% add_recipe(rec) %>% add_model(en_spec),
                         grid = en_grid),
    random_forest = list(wf = workflow() %>% add_recipe(rec) %>% add_model(rf_spec),
                         grid = rf_grid)
  )
}

# ---------------------------------------------------------------------------
# tune_and_finalize(): cross-validate one workflow over its grid, pick the best
# hyperparameters by `select_metric`, and refit the finalized workflow on the
# full training frame. Returns the fitted workflow plus its best CV metrics.
# ---------------------------------------------------------------------------
tune_and_finalize <- function(model, resamples, train_df, select_metric = "roc_auc") {
  # roc_auc needs class probabilities; include it alongside the class metrics.
  mset <- metric_set(accuracy, bal_accuracy, kap, f_meas, roc_auc, mn_log_loss)
  tuned <- tune_grid(
    model$wf, resamples = resamples, grid = model$grid,
    metrics = mset,
    control = control_grid(save_pred = FALSE))

  best <- select_best(tuned, metric = select_metric)   # carries a .config id
  final_wf <- finalize_workflow(model$wf, best)
  fit <- fit(final_wf, data = train_df)

  # Keep only the winning configuration's cross-validated metrics.
  cv_metrics <- collect_metrics(tuned) %>%
    filter(.config == best$.config) %>% 
    select(.metric, mean, std_err,
           dplyr::any_of(c("penalty", "mixture", "mtry", "min_n")))

  list(fit = fit, best = best, cv_metrics = cv_metrics, tuned = tuned)
}

# ---------------------------------------------------------------------------
# evaluate(): apply a fitted workflow to `new_df` and compute the class metrics,
# multiclass ROC AUC, and a confusion matrix against the true `cluster`.
# `label` tags the output so results from different eval sets stay identifiable.
# ---------------------------------------------------------------------------
evaluate <- function(fit, new_df, label) {
  pred_class <- predict(fit, new_df, type = "class")
  pred_prob  <- predict(fit, new_df, type = "prob")
  res <- bind_cols(select(new_df, cluster), pred_class, pred_prob)

  prob_cols <- grep("^\\.pred_C", names(res), value = TRUE)
  metrics <- bind_rows(
    CLASS_METRICS(res, truth = cluster, estimate = .pred_class,
                  estimator = "macro"),
    roc_auc(res, truth = cluster, all_of(prob_cols), estimator = "hand_till")
  ) %>%
    mutate(eval_set = label, .before = 1)

  cm <- conf_mat(res, truth = cluster, estimate = .pred_class)

  message("\n== Evaluation: ", label, " ==")
  print(metrics)
  print(cm)
  list(metrics = metrics, conf_mat = cm, predictions = res)
}
