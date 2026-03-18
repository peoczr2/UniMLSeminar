run_scenario_4_low_corr_tau_z_x2_opposite_sign <- function(repo_root, output_dir, n_train, n_test, num_trees, beat_penalty, target_share, seed) {
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

  tau_raw <- pmax(x1, 0) - 0.5 * x3 + z1 - x2
  tau <- as.numeric(scale(tau_raw))

  W <- rbinom(n_total, size = 1, prob = 0.5)
  Y <- tau * W + rnorm(n_total, mean = 0, sd = 1)

  train_idx <- seq_len(n_train)
  test_idx <- seq.int(n_train + 1, n_total)

  sim_data <- list(
    train = list(
      X = X[train_idx, , drop = FALSE],
      Z = Z[train_idx, , drop = FALSE],
      X_full = cbind(X[train_idx, , drop = FALSE], Z[train_idx, , drop = FALSE]),
      W = W[train_idx],
      Y = Y[train_idx],
      tau = tau[train_idx]
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
    fit_beat(sim_data, num_trees, beat_penalty, seed + 30)
  )

  cf_fd_targeted <- target_top_share(method_results[[1]]$score, target_share)
  cf_fd_imbalance <- imbalance_metric(cf_fd_targeted, scenario_imbalance_input(sim_data$test))

  metrics_df <- do.call(
    rbind,
    lapply(method_results, function(result) {
      evaluate_method(sim_data$test, result, target_share, cf_fd_imbalance)
    })
  )
  metrics_df$scenario <- "scenario_4_low_corr_tau_z_x2_opposite_sign"

  scores_df <- build_scores_frame(sim_data$test, method_results)
  plot_path <- file.path(output_dir, "scenario_4_low_corr_tau_z_x2_opposite_sign.png")
  plot_method_scores(scores_df, plot_path, "Scenario 4: Low corr, tau = f(Z, -X2)")

  metrics_df
}

if (sys.nframe() == 0) {
  run_scenario_script(
    run_scenario_4_low_corr_tau_z_x2_opposite_sign,
    "scenario_4_low_corr_tau_z_x2_opposite_sign"
  )
}
