run_scenario_4_low_corr_tau_z_x2_opposite_sign <- function(repo_root, output_dir, n_train, n_test, num_trees, beat_penalty, target_share, seed, methods_to_run = NULL) {
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

  if (is.null(methods_to_run)) {
    methods_to_run <- c("CF-FD", "CF-NP", "DEBIASED", "BEAT", "BTGQ", "BTGQ_VALID", "BTGQ_MAX")
  }
  methods_to_run <- toupper(methods_to_run)
  include_method <- function(name) name %in% methods_to_run

  cf_fd_result <- fit_cf_fd(sim_data, num_trees, seed)
  method_results <- list()
  if (include_method("CF-FD")) method_results[[length(method_results) + 1L]] <- cf_fd_result
  if (include_method("CF-NP")) method_results[[length(method_results) + 1L]] <- fit_cf_np(sim_data, num_trees, seed + 10)
  if (include_method("DEBIASED")) method_results[[length(method_results) + 1L]] <- fit_debiased(sim_data, num_trees, seed + 20)
  if (include_method("BEAT")) method_results[[length(method_results) + 1L]] <- fit_beat(sim_data, num_trees, beat_penalty, seed + 30)
  if (include_method("BTGQ")) method_results[[length(method_results) + 1L]] <- fit_btgq(sim_data, num_trees, seed + 40, budget = 0.5, target_quota = 0.5, use_validation = FALSE, method_name = "BTGQ")
  if (include_method("BTGQ_VALID")) method_results[[length(method_results) + 1L]] <- fit_btgq(sim_data, num_trees, seed + 50, budget = 0.5, target_quota = 0.5, use_validation = TRUE, method_name = "BTGQ_VALID")
  if (include_method("BTGQ_MAX")) method_results[[length(method_results) + 1L]] <- fit_btgq_max(sim_data, num_trees, seed + 60, budget_max = 0.5, target_quota = 0.5, use_validation = FALSE, method_name = "BTGQ_MAX")

  result_names <- vapply(method_results, `[[`, character(1), "method")
  if ("BTGQ" %in% result_names) {
    btgq_plot_path <- file.path(output_dir, sprintf("scenario_4_low_corr_tau_z_x2_opposite_sign_btgq_lambda_trace_seed_%s.png", seed))
    plot_btgq_lambda_trace(
      method_results[[which(result_names == "BTGQ")[1]]]$lambda_trace,
      btgq_plot_path,
      "Scenario 4 BTGQ: Lambda vs Targeted Group Demo",
      target_quota = 0.5
    )
  }

  if ("BTGQ_VALID" %in% result_names) {
    btgq_valid_plot_path <- file.path(output_dir, sprintf("scenario_4_low_corr_tau_z_x2_opposite_sign_btgq_valid_lambda_trace_seed_%s.png", seed))
    plot_btgq_lambda_trace(
      method_results[[which(result_names == "BTGQ_VALID")[1]]]$lambda_trace,
      btgq_valid_plot_path,
      "Scenario 4 BTGQ_VALID: Lambda vs Targeted Group Demo",
      target_quota = 0.5
    )
  }

  if ("BTGQ_MAX" %in% result_names) {
    btgq_max_plot_path <- file.path(output_dir, sprintf("scenario_4_low_corr_tau_z_x2_opposite_sign_btgq_max_lambda_trace_seed_%s.png", seed))
    plot_btgq_lambda_trace(
      method_results[[which(result_names == "BTGQ_MAX")[1]]]$lambda_trace,
      btgq_max_plot_path,
      "Scenario 4 BTGQ_MAX: Lambda vs Targeted Group Demo",
      target_quota = 0.5,
      color_by = "budget.actual",
      color_label = "Budget size"
    )
  }

  cf_fd_targeted <- target_top_share(cf_fd_result$score, target_share)
  cf_fd_imbalance <- imbalance_metric(cf_fd_targeted, scenario_imbalance_input(sim_data$test))

  metrics_df <- do.call(
    rbind,
    lapply(method_results, function(result) {
      evaluate_method(sim_data$test, result, target_share, cf_fd_imbalance)
    })
  )
  metrics_df$scenario <- "scenario_4_low_corr_tau_z_x2_opposite_sign"

  if (length(method_results) > 0) {
    scores_df <- build_scores_frame(sim_data$test, method_results)
    plot_path <- file.path(output_dir, "scenario_4_low_corr_tau_z_x2_opposite_sign.png")
    plot_method_scores(scores_df, plot_path, "Scenario 4: Low corr, tau = f(Z, -X2)")
  }

  metrics_df
}

if (sys.nframe() == 0) {
  run_scenario_script(
    run_scenario_4_low_corr_tau_z_x2_opposite_sign,
    "scenario_4_low_corr_tau_z_x2_opposite_sign"
  )
}
