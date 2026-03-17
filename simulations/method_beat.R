fit_beat <- function(sim_data, num_trees, beat_penalty, seed) {
  fit <- beat::balanced_causal_forest(
    sim_data$train$X,
    sim_data$train$Y,
    sim_data$train$W,
    target.weights = as.matrix(sim_data$train$Z),
    target.weight.penalty = beat_penalty,
    num.trees = num_trees,
    seed = seed
  )

  score <- as.numeric(predict(fit, sim_data$test$X)$predictions)

  list(
    method = "BEAT",
    score = score,
    twin_score = score
  )
}
