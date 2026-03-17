args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args[grep(file_arg, args)])

if (length(script_path) != 1 || !nzchar(script_path)) {
  stop("Run this script with Rscript scripts/run_btgq_smoke_test.R", call. = FALSE)
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

n_train <- 240
n_test <- 120
n <- n_train + n_test

X_cont <- matrix(rnorm(n * 4), nrow = n, ncol = 4)
X_disc <- matrix(rbinom(n * 3, 1, 0.35), nrow = n, ncol = 3)
X <- cbind(X_cont, X_disc)
target_group <- rbinom(n, 1, plogis(0.7 * X_cont[, 2] - 0.4 * X_cont[, 3]))
tau <- 0.3 + 1.0 * pmax(X[, 1], 0) + 0.6 * X[, 5] + 0.8 * target_group
W <- rbinom(n, 1, 0.5)
Y0 <- 0.5 * X[, 1] - X[, 2] + 0.8 * X[, 4] + 0.4 * target_group + rnorm(n)
Y <- Y0 + tau * W

train_idx <- seq_len(n_train)
test_idx <- seq.int(n_train + 1, n)

baseline_fit <- btgq_causal_forest(
  X[train_idx, , drop = FALSE],
  Y[train_idx],
  W[train_idx],
  target.group = target_group[train_idx],
  btgq.lambda = 0,
  num.trees = 200
)

tuned <- tune_btgq_causal_forest(
  X[train_idx, , drop = FALSE],
  Y[train_idx],
  W[train_idx],
  target.group = target_group[train_idx],
  budget = 0.30,
  target.quota = 0.65,
  num.trees = 200,
  max.search.iter = 6
)

baseline_stats <- {
  pred <- predict(baseline_fit)$predictions
  top_n <- max(1, ceiling(0.30 * length(pred)))
  idx <- order(-pred, seq_along(pred))[seq_len(top_n)]
  mean(target_group[train_idx][idx])
}

tuned_stats <- {
  pred <- predict(tuned$forest)$predictions
  top_n <- max(1, ceiling(0.30 * length(pred)))
  idx <- order(-pred, seq_along(pred))[seq_len(top_n)]
  mean(target_group[train_idx][idx])
}

test_pred <- predict(tuned$forest, X[test_idx, , drop = FALSE])$predictions

cat("BTGQ smoke test succeeded.\n")
cat("Selected lambda:", tuned$selected.lambda, "\n")
cat("Quota attained:", tuned$quota.attained, "\n")
cat("Baseline top-budget quota:", baseline_stats, "\n")
cat("Tuned top-budget quota:", tuned_stats, "\n")
cat("Held-out predictions summary:\n")
print(summary(test_pred))
cat("\nSearch trace:\n")
print(tuned$search.trace)
