args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args[grep(file_arg, args)])

if (length(script_path) != 1 || !nzchar(script_path)) {
  stop("Run this script with Rscript scripts/run_beat_smoke_test.R", call. = FALSE)
}

script_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = TRUE))
repo_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
lib_dir <- file.path(repo_root, ".r-library")

if (!dir.exists(lib_dir)) {
  stop("Local R library not found. Run Rscript scripts/setup_beat.R first.", call. = FALSE)
}

.libPaths(c(lib_dir, .libPaths()))

suppressPackageStartupMessages(library(beat))

set.seed(42)

n_train <- 200
n_test <- 100
n <- n_train + n_test

X_cont <- matrix(rnorm(n * 4), nrow = n, ncol = 4)
X_disc <- matrix(rbinom(n * 3, 1, 0.3), nrow = n, ncol = 3)
X <- cbind(X_cont, X_disc)
Z <- rbinom(n, 1, 1 / (1 + exp(-X_cont[, 2])))
tau <- -1 + pmax(X[, 1], 0) + X[, 2] + abs(X[, 3]) + X[, 5]
W <- rbinom(n, 1, 0.5)
Y_r <- X[, 1] - 2 * X[, 2] + X[, 4] + 3 * Z + runif(n)
Y <- Y_r + tau * W

train_idx <- seq_len(n_train)
test_idx <- seq.int(n_train + 1, n)

fit <- balanced_causal_forest(
  X[train_idx, , drop = FALSE],
  Y[train_idx],
  W[train_idx],
  target.weights = as.matrix(Z[train_idx]),
  target.weight.penalty = 10,
  num.trees = 200
)

pred <- predict(fit, X[test_idx, , drop = FALSE])$predictions

cat("Smoke test succeeded.\n")
cat("Predictions summary:\n")
print(summary(pred))
