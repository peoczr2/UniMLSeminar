fit_cf_fd <- function(sim_data, num_trees, seed) {
  fit <- beat::causal_forest(
    sim_data$train$X_full,
    sim_data$train$Y,
    sim_data$train$W,
    num.trees = num_trees,
    seed = seed
  )

  twin_full <- cbind(sim_data$test$X, flip_protected(sim_data$test$Z))

  list(
    method = "CF-FD",
    score = as.numeric(predict(fit, sim_data$test$X_full)$predictions),
    twin_score = as.numeric(predict(fit, twin_full)$predictions)
  )
}
