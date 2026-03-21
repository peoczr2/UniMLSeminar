args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args[grep(file_arg, args)])

if (length(script_path) != 1 || !nzchar(script_path)) {
  stop("Run this script with Rscript scripts/run_btgq_max_smoke_test.R", call. = FALSE)
}

script_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = TRUE))
repo_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
lib_dir <- file.path(repo_root, ".r-library")

if (!dir.exists(lib_dir)) {
  stop("Local R library not found. Run Rscript scripts/setup_beat.R first.", call. = FALSE)
}

.libPaths(c(lib_dir, .libPaths()))

suppressPackageStartupMessages(library(beat))

set.seed(7)

n_train <- 180
n_valid <- 60
n_test <- 120
n <- n_train + n_valid + n_test

X_cont <- matrix(rnorm(n * 4), nrow = n, ncol = 4)
X_disc <- matrix(rbinom(n * 3, 1, 0.35), nrow = n, ncol = 3)
X <- cbind(X_cont, X_disc)
target_group <- rbinom(n, 1, plogis(0.7 * X_cont[, 2] - 0.4 * X_cont[, 3]))
tau <- 0.3 + 1.0 * pmax(X[, 1], 0) + 0.6 * X[, 5] + 0.8 * target_group
W <- rbinom(n, 1, 0.5)
Y0 <- 0.5 * X[, 1] - X[, 2] + 0.8 * X[, 4] + 0.4 * target_group + rnorm(n)
Y <- Y0 + tau * W

train_idx <- seq_len(n_train)
valid_idx <- seq.int(n_train + 1, n_train + n_valid)
test_idx <- seq.int(n_train + n_valid + 1, n)

baseline_fit <- btgq_max_causal_forest(
  X[train_idx, , drop = FALSE],
  Y[train_idx],
  W[train_idx],
  target.group = target_group[train_idx],
  btgq.lambda = 0,
  num.trees = 200
)

tuned <- tune_btgq_max_causal_forest(
  X[train_idx, , drop = FALSE],
  Y[train_idx],
  W[train_idx],
  target.group = target_group[train_idx],
  budget.max = 0.30,
  target.quota = 0.65,
  validation.X = X[valid_idx, , drop = FALSE],
  validation.Y = Y[valid_idx],
  validation.W = W[valid_idx],
  validation.target.group = target_group[valid_idx],
  num.trees = 200,
  max.search.iter = 6
)

validation_stats <- {
  pred <- predict(baseline_fit, X[valid_idx, , drop = FALSE])$predictions
  stats <- beat:::btgq_max_policy_stats(pred, target_group[valid_idx], budget.max = 0.30, viability.epsilon = 1e-08)
  c(quota = stats$quota, budget_actual = stats$budget.actual, viable_share = stats$viable.share)
}

tuned_stats <- {
  pred <- predict(tuned$forest, X[valid_idx, , drop = FALSE])$predictions
  stats <- beat:::btgq_max_policy_stats(pred, target_group[valid_idx], budget.max = 0.30, viability.epsilon = 1e-08)
  c(quota = stats$quota, budget_actual = stats$budget.actual, viable_share = stats$viable.share)
}

test_pred <- predict(tuned$forest, X[test_idx, , drop = FALSE])$predictions

cat("BTGQ-Max smoke test succeeded.\n")
cat("Selected lambda:", tuned$selected.lambda, "\n")
cat("Quota attained:", tuned$quota.attained, "\n")
cat("Baseline quota / budget / viable share:\n")
print(validation_stats)
cat("Tuned quota / budget / viable share:\n")
print(tuned_stats)
cat("Held-out predictions summary:\n")
print(summary(test_pred))
cat("\nSearch trace:\n")
print(tuned$search.trace)
