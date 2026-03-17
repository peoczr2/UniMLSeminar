fit_debiased <- function(sim_data, num_trees, seed) {
  debias_forest <- beat::multi_regression_forest(
    sim_data$train$Z,
    sim_data$train$X,
    num.trees = num_trees,
    seed = seed
  )

  x_hat_train <- predict(debias_forest)$predictions
  x_hat_test <- predict(debias_forest, sim_data$test$Z)$predictions
  x_hat_twin <- predict(debias_forest, flip_protected(sim_data$test$Z))$predictions

  x_debiased_train <- sim_data$train$X - x_hat_train
  x_debiased_test <- sim_data$test$X - x_hat_test
  x_debiased_twin <- sim_data$test$X - x_hat_twin

  fit <- beat::causal_forest(
    x_debiased_train,
    sim_data$train$Y,
    sim_data$train$W,
    num.trees = num_trees,
    seed = seed + 1
  )

  list(
    method = "Debiased",
    score = as.numeric(predict(fit, x_debiased_test)$predictions),
    twin_score = as.numeric(predict(fit, x_debiased_twin)$predictions)
  )
}
