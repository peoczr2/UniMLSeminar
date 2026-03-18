fit_btgq <- function(
    sim_data,
    num_trees,
    seed,
    budget = 0.5,
    target_quota = 0.5,
    validation_data = NULL,
    use_validation = TRUE,
    method_name = "BTGQ") {
  fit_data <- sim_data$train
  final_fit_data <- fit_data
  validation_X <- NULL
  validation_target_group <- NULL

  if (use_validation && !is.null(sim_data$btgq_valid_train)) {
    fit_data <- sim_data$btgq_valid_train
    final_fit_data <- sim_data$train
  }

  if (is.null(validation_data) && use_validation && !is.null(sim_data$valid)) {
    validation_data <- sim_data$valid
  }

  if (!is.null(validation_data)) {
    validation_X <- validation_data$X
    validation_target_group <- validation_data$Z[, "Z1"]
  }

  tuned <- beat::tune_btgq_causal_forest(
    fit_data$X,
    fit_data$Y,
    fit_data$W,
    target.group = fit_data$Z[, "Z1"],
    budget = budget,
    target.quota = target_quota,
    validation.X = validation_X,
    validation.target.group = validation_target_group,
    num.trees = num_trees,
    seed = seed
  )

  final_forest <- tuned$forest
  if (use_validation) {
    final_forest <- beat::btgq_causal_forest(
      final_fit_data$X,
      final_fit_data$Y,
      final_fit_data$W,
      target.group = final_fit_data$Z[, "Z1"],
      btgq.lambda = tuned$selected.lambda,
      num.trees = num_trees,
      seed = seed
    )
  }

  score <- as.numeric(predict(final_forest, sim_data$test$X)$predictions)

  list(
    method = method_name,
    score = score,
    twin_score = score,
    selected_lambda = tuned$selected.lambda,
    lambda_trace = transform(
      tuned$search.trace,
      targeted_group_demo = achieved.quota
    )
  )
}
