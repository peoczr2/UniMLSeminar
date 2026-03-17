run_scenario_4_low_corr_tau_z_x2_opposite_sign <- function(repo_root, output_dir, n_train, n_test, num_trees, beat_penalty, target_share, seed) {
  set.seed(seed)

  features <- make_base_features(n_train + n_test, corr_strength = 1)
  tau <- as.numeric(scale(pmax(features$X[, 1], 0) - 0.5 * features$X[, 3] + features$Z[, "Z1"] - features$X[, 2]))
  sim_data <- make_standard_train_test(features$X, features$Z, tau, n_train, n_test)

  method_results <- list(
    fit_cf_fd(sim_data, num_trees, seed),
    fit_cf_np(sim_data, num_trees, seed + 10),
    fit_debiased(sim_data, num_trees, seed + 20),
    fit_beat(sim_data, num_trees, beat_penalty, seed + 30)
  )

  cf_fd_targeted <- target_top_share(method_results[[1]]$score, target_share)
  cf_fd_imbalance <- imbalance_metric(cf_fd_targeted, sim_data$test$Z[, "Z1"])

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
