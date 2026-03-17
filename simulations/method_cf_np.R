fit_cf_np <- function(sim_data, num_trees, seed) {
  fit <- beat:::causal_forest(
    sim_data$train$X,
    sim_data$train$Y,
    sim_data$train$W,
    num.trees = num_trees,
    seed = seed
  )

  score <- as.numeric(beat:::predict.causal_forest(fit, sim_data$test$X)$predictions)

  list(
    method = "CF-NP",
    score = score,
    twin_score = score
  )
}
