run_scenario_2_low_corr_tau_z <- function(repo_root, output_dir, n_train, n_test, num_trees, beat_penalty, target_share, seed) {
  set.seed(seed)

  n_total <- n_train + n_test

  x1 <- rnorm(n_total, mean = 0, sd = 1)
  x2 <- rnorm(n_total, mean = 0, sd = 1)
  x3 <- rnorm(n_total, mean = 0, sd = 1)
  x4 <- rnorm(n_total, mean = 0, sd = 1)
  x5 <- rbinom(n_total, size = 1, prob = 0.3)
  x6 <- rbinom(n_total, size = 1, prob = 0.3)
  x7 <- rbinom(n_total, size = 1, prob = 0.3)
  x8 <- rbinom(n_total, size = 1, prob = 0.3)
  x9 <- rbinom(n_total, size = 1, prob = 0.3)
  x10 <- rbinom(n_total, size = 1, prob = 0.3)
  X <- cbind(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10)
  colnames(X) <- paste0("X", seq_len(ncol(X)))

  z1_prob <- plogis(1 * x2)
  z1 <- rbinom(n_total, size = 1, prob = z1_prob)
  z2 <- rnorm(n_total, mean = 0, sd = 1)
  z3 <- rnorm(n_total, mean = 0, sd = 1)
  z4 <- rbinom(n_total, size = 1, prob = 0.3)
  Z <- cbind(Z1 = z1, Z2 = z2, Z3 = z3, Z4 = z4)

  tau_raw <- pmax(x1, 0) - 0.5 * x3 + z1
  tau <- as.numeric(scale(tau_raw))

  W <- rbinom(n_total, size = 1, prob = 0.5)
  Y <- tau * W + rnorm(n_total, mean = 0, sd = 1)

  train_idx <- seq_len(n_train)
  test_idx <- seq.int(n_train + 1, n_total)
  train_perm <- sample(train_idx)
  valid_size <- max(1L, floor(0.2 * n_train))
  valid_idx <- train_perm[seq_len(valid_size)]
  fit_idx <- train_perm[-seq_len(valid_size)]

  sim_data <- list(
    train = list(
      X = X[train_idx, , drop = FALSE],
      Z = Z[train_idx, , drop = FALSE],
      X_full = cbind(X[train_idx, , drop = FALSE], Z[train_idx, , drop = FALSE]),
      W = W[train_idx],
      Y = Y[train_idx],
      tau = tau[train_idx]
    ),
    valid = list(
      X = X[valid_idx, , drop = FALSE],
      Z = Z[valid_idx, , drop = FALSE],
      X_full = cbind(X[valid_idx, , drop = FALSE], Z[valid_idx, , drop = FALSE]),
      W = W[valid_idx],
      Y = Y[valid_idx],
      tau = tau[valid_idx]
    ),
    btgq_valid_train = list(
      X = X[fit_idx, , drop = FALSE],
      Z = Z[fit_idx, , drop = FALSE],
      X_full = cbind(X[fit_idx, , drop = FALSE], Z[fit_idx, , drop = FALSE]),
      W = W[fit_idx],
      Y = Y[fit_idx],
      tau = tau[fit_idx]
    ),
    test = list(
      X = X[test_idx, , drop = FALSE],
      Z = Z[test_idx, , drop = FALSE],
      X_full = cbind(X[test_idx, , drop = FALSE], Z[test_idx, , drop = FALSE]),
      W = W[test_idx],
      Y = Y[test_idx],
      tau = tau[test_idx]
    )
  )

  method_results <- list(
    fit_cf_fd(sim_data, num_trees, seed),
    fit_cf_np(sim_data, num_trees, seed + 10),
    fit_debiased(sim_data, num_trees, seed + 20),
    fit_beat(sim_data, num_trees, beat_penalty, seed + 30),
    fit_btgq(sim_data, num_trees, seed + 40, budget = 0.5, target_quota = 0.5, use_validation = FALSE, method_name = "BTGQ"),
    fit_btgq(sim_data, num_trees, seed + 50, budget = 0.5, target_quota = 0.5, use_validation = TRUE, method_name = "BTGQ_VALID"),
    fit_btgq_dumb(sim_data, num_trees, seed + 60, budget = 0.5, target_quota = 0.5, use_validation = FALSE, method_name = "BTGQ_DUMB")
  )

  btgq_plot_path <- file.path(output_dir, sprintf("scenario_2_low_corr_tau_z_btgq_lambda_trace_seed_%s.png", seed))
  plot_btgq_lambda_trace(
    method_results[[5]]$lambda_trace,
    btgq_plot_path,
    "Scenario 2 BTGQ: Lambda vs Targeted Group Demo",
    target_quota = 0.5
  )

  btgq_valid_plot_path <- file.path(output_dir, sprintf("scenario_2_low_corr_tau_z_btgq_valid_lambda_trace_seed_%s.png", seed))
  plot_btgq_lambda_trace(
    method_results[[6]]$lambda_trace,
    btgq_valid_plot_path,
    "Scenario 2 BTGQ_VALID: Lambda vs Targeted Group Demo",
    target_quota = 0.5
  )

  btgq_dumb_plot_path <- file.path(output_dir, sprintf("scenario_2_low_corr_tau_z_btgq_dumb_lambda_trace_seed_%s.png", seed))
  plot_btgq_lambda_trace(
    method_results[[7]]$lambda_trace,
    btgq_dumb_plot_path,
    "Scenario 2 BTGQ_DUMB: Lambda vs Targeted Group Demo",
    target_quota = 0.5
  )

  cf_fd_targeted <- target_top_share(method_results[[1]]$score, target_share)
  cf_fd_imbalance <- imbalance_metric(cf_fd_targeted, scenario_imbalance_input(sim_data$test))

  metrics_df <- do.call(
    rbind,
    lapply(method_results, function(result) {
      evaluate_method(sim_data$test, result, target_share, cf_fd_imbalance)
    })
  )
  metrics_df$scenario <- "scenario_2_low_corr_tau_z"

  scores_df <- build_scores_frame(sim_data$test, method_results)
  plot_path <- file.path(output_dir, "scenario_2_low_corr_tau_z.png")
  plot_method_scores(scores_df, plot_path, "Scenario 2: Low corr, tau = f(Z)")

  metrics_df
}

if (sys.nframe() == 0) {
  run_scenario_script(run_scenario_2_low_corr_tau_z, "scenario_2_low_corr_tau_z")
}
