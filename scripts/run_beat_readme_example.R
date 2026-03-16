args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args[grep(file_arg, args)])

if (length(script_path) != 1 || !nzchar(script_path)) {
  stop("Run this script with Rscript scripts/run_beat_readme_example.R", call. = FALSE)
}

script_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = TRUE))
repo_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
lib_dir <- file.path(repo_root, ".r-library")
output_dir <- file.path(repo_root, "outputs")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib_dir, .libPaths()))

suppressPackageStartupMessages(library(beat))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))

set.seed(1)

n1 <- 1000
n2 <- 1000
p_continuous <- 4
p_discrete <- 3
n <- n1 + n2

X_cont <- matrix(rnorm(n * p_continuous), n, p_continuous)
X_disc <- matrix(rbinom(n * p_discrete, 1, 0.3), n, p_discrete)
X <- cbind(X_cont, X_disc)
Z <- rbinom(n, 1, 1 / (1 + exp(-X_cont[, 2])))
tau <- -1 + pmax(X[, 1], 0) + X[, 2] + abs(X[, 3]) + X[, 5]
W <- rbinom(n, 1, 0.5)
Y_r <- X[, 1] - 2 * X[, 2] + X[, 4] + 3 * Z + runif(n)
Y <- Y_r + tau * W

train_data <- data.frame(
  Y = Y[1:n1],
  Z = Z[1:n1],
  W = W[1:n1],
  X = X[1:n1, ],
  tau = tau[1:n1],
  Y_r = Y_r[1:n1]
)

test_data <- data.frame(
  Y = Y[(n1 + 1):(n1 + n2)],
  Z = Z[(n1 + 1):(n1 + n2)],
  W = W[(n1 + 1):(n1 + n2)],
  X = X[(n1 + 1):(n1 + n2), ],
  tau = tau[(n1 + 1):(n1 + n2)],
  Y_r = Y_r[(n1 + 1):(n1 + n2)]
)

X_train <- train_data[, 4:10]
W_train <- train_data$W
Z_train <- train_data[, 2]
Y_train <- train_data$Y
Y_r_train <- train_data$Y_r
X_test <- test_data[, 4:10]
Z_test <- test_data$Z

num_trees <- 2000
my_penalty <- 10

fit_causal_beat <- balanced_causal_forest(
  X_train,
  Y_train,
  W_train,
  target.weights = as.matrix(Z_train),
  target.weight.penalty = my_penalty,
  num.trees = num_trees
)

cbt_causal_train <- predict(fit_causal_beat)$predictions
cbt_causal_test <- predict(fit_causal_beat, X_test)$predictions

fit_regression_beat <- balanced_regression_forest(
  X_train,
  Y_r_train,
  target.weights = as.matrix(Z_train),
  target.weight.penalty = my_penalty,
  num.trees = num_trees
)

cbt_regression_train <- predict(fit_regression_beat)$predictions
cbt_regression_test <- predict(fit_regression_beat, X_test)$predictions

dat_plot <- data.table(
  cbt_causal = cbt_causal_test,
  cbt_regr = cbt_regression_test,
  true_causal = test_data$tau,
  true_reg = test_data$Y_r,
  Z = as.factor(Z_test)
)

make_density_plot <- function(data, xvar, title_text) {
  ggplot(data, aes_string(x = xvar, color = "Z", fill = "Z")) +
    geom_density(alpha = 0.2) +
    geom_vline(
      data = data[, .(mean_value = mean(get(xvar))), by = Z],
      aes(xintercept = mean_value, color = Z),
      linetype = "dashed",
      linewidth = 0.7,
      show.legend = FALSE
    ) +
    labs(title = title_text, x = xvar, y = "Density") +
    theme_minimal()
}

p1 <- make_density_plot(dat_plot, "true_causal", "true causal")
p2 <- make_density_plot(dat_plot, "cbt_causal", "cbt causal")
p3 <- make_density_plot(dat_plot, "true_reg", "true regression")
p4 <- make_density_plot(dat_plot, "cbt_regr", "cbt regression")

png(file.path(output_dir, "beat_readme_example.png"), width = 1600, height = 1200, res = 150)
grid.arrange(p1, p2, p3, p4, ncol = 2)
dev.off()

summary_table <- dat_plot[, .(
  mean_true_causal = mean(true_causal),
  mean_cbt_causal = mean(cbt_causal),
  mean_true_reg = mean(true_reg),
  mean_cbt_reg = mean(cbt_regr)
), by = Z]

cat("README example completed.\n")
cat("Causal prediction summary:\n")
print(summary(cbt_causal_test))
cat("\nRegression prediction summary:\n")
print(summary(cbt_regression_test))
cat("\nGroup means by protected attribute Z:\n")
print(summary_table)
cat("\nPlot saved to:", file.path(output_dir, "beat_readme_example.png"), "\n")
