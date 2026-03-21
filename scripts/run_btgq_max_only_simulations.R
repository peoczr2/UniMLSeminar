args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args[grep(file_arg, args)])

if (length(script_path) != 1 || !nzchar(script_path)) {
  stop("Run this script with Rscript scripts/run_btgq_max_only_simulations.R", call. = FALSE)
}

repo_root <- normalizePath(
  file.path(dirname(normalizePath(script_path, winslash = "/", mustWork = TRUE)), ".."),
  winslash = "/",
  mustWork = TRUE
)

source(file.path(repo_root, "simulations", "helpers.R"))
load_simulation_runtime(repo_root)
source(file.path(repo_root, "simulations", "scenario_1_high_corr_tau_z.R"))
source(file.path(repo_root, "simulations", "scenario_2_low_corr_tau_z.R"))
source(file.path(repo_root, "simulations", "scenario_3_low_corr_tau_z_x2_same_sign.R"))
source(file.path(repo_root, "simulations", "scenario_4_low_corr_tau_z_x2_opposite_sign.R"))
ensure_local_library(repo_root)
suppressPackageStartupMessages(library(beat))

config <- parse_cli_args(list(
  reps = 3,
  n_train = 2000,
  n_test = 1000,
  num_trees = 250,
  target_share = 0.5,
  beat_penalty = 10,
  output_dir = file.path(repo_root, "outputs", "simulations")
))

dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)

scenario_runners <- list(
  run_scenario_1_high_corr_tau_z,
  run_scenario_2_low_corr_tau_z,
  run_scenario_3_low_corr_tau_z_x2_same_sign,
  run_scenario_4_low_corr_tau_z_x2_opposite_sign
)
scenario_names <- c(
  "scenario_1_high_corr_tau_z",
  "scenario_2_low_corr_tau_z",
  "scenario_3_low_corr_tau_z_x2_same_sign",
  "scenario_4_low_corr_tau_z_x2_opposite_sign"
)

updated_rows <- list()
for (i in seq_along(scenario_runners)) {
  for (rep_id in seq_len(config$reps)) {
    seed <- 1000 + rep_id + i * 100
    result <- scenario_runners[[i]](
      repo_root = repo_root,
      output_dir = config$output_dir,
      n_train = as.integer(config$n_train),
      n_test = as.integer(config$n_test),
      num_trees = as.integer(config$num_trees),
      beat_penalty = config$beat_penalty,
      target_share = config$target_share,
      seed = seed,
      methods_to_run = c("BTGQ_MAX")
    )
    result$replication <- rep_id
    result$n_train <- as.integer(config$n_train)
    result$n_test <- as.integer(config$n_test)
    result$num_trees <- as.integer(config$num_trees)
    result$beat_penalty <- config$beat_penalty
    result$target_share <- config$target_share
    result$scenario <- scenario_names[[i]]
    updated_rows[[length(updated_rows) + 1L]] <- result
  }
}

new_rows <- do.call(rbind, updated_rows)
existing_results <- read.csv(file.path(config$output_dir, "simulation_results.csv"))
remaining_results <- existing_results[existing_results$method != "BTGQ_MAX", , drop = FALSE]
results_df <- rbind(remaining_results, new_rows)

mean_or_na <- function(x) {
  if (all(is.na(x))) {
    return(NA_real_)
  }
  mean(x, na.rm = TRUE)
}

summary_df <- aggregate(
  results_df[, c("efficiency", "raw_imbalance", "imbalance", "targeted_group_demo", "budget_actual", "maximum_viable_ceiling_z0", "maximum_viable_ceiling_z1", "delta_policy")],
  by = list(scenario = results_df$scenario, method = results_df$method),
  FUN = mean_or_na
)

write.csv(results_df, file.path(config$output_dir, "simulation_results.csv"), row.names = FALSE)
write.csv(summary_df, file.path(config$output_dir, "simulation_summary.csv"), row.names = FALSE)

cat("BTGQ-Max-only simulation refresh completed.\n")
print(subset(summary_df, method == "BTGQ_MAX"))
