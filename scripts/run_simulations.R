args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args[grep(file_arg, args)])

if (length(script_path) != 1 || !nzchar(script_path)) {
  stop("Run this script with Rscript scripts/run_simulations.R", call. = FALSE)
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

config <- parse_cli_args(list(
  reps = 3,
  n_train = 2000,
  n_test = 1000,
  num_trees = 500,
  target_share = 0.5,
  beat_penalty = 10,
  methods = "all",
  use_cache = TRUE,
  output_dir = file.path(repo_root, "outputs", "simulations")
))

dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
ensure_local_library(repo_root)
selected_methods <- parse_method_selection(config$methods)

suppressPackageStartupMessages(library(beat))

all_results <- list()
scenario_runners <- list(
  scenario_1_high_corr_tau_z = run_scenario_1_high_corr_tau_z,
  scenario_2_low_corr_tau_z = run_scenario_2_low_corr_tau_z,
  scenario_3_low_corr_tau_z_x2_same_sign = run_scenario_3_low_corr_tau_z_x2_same_sign,
  scenario_4_low_corr_tau_z_x2_opposite_sign = run_scenario_4_low_corr_tau_z_x2_opposite_sign
)

scenario_counter <- 1L
for (scenario_name in names(scenario_runners)) {
  for (rep_id in seq_len(config$reps)) {
    seed <- 1000 + rep_id + scenario_counter * 100
    result <- scenario_runners[[scenario_name]](
      repo_root = repo_root,
      output_dir = config$output_dir,
      n_train = as.integer(config$n_train),
      n_test = as.integer(config$n_test),
      num_trees = as.integer(config$num_trees),
      beat_penalty = config$beat_penalty,
      target_share = config$target_share,
      seed = seed,
      selected_methods = selected_methods,
      use_cache = isTRUE(config$use_cache)
    )
    metrics_df <- result
    metrics_df$replication <- rep_id
    metrics_df$n_train <- as.integer(config$n_train)
    metrics_df$n_test <- as.integer(config$n_test)
    metrics_df$num_trees <- as.integer(config$num_trees)
    metrics_df$beat_penalty <- config$beat_penalty
    metrics_df$target_share <- config$target_share
    all_results[[length(all_results) + 1L]] <- metrics_df
  }
  scenario_counter <- scenario_counter + 1L
}

results_df <- do.call(rbind, all_results)
summary_df <- aggregate(
  cbind(efficiency, raw_imbalance, imbalance, targeted_group_demo, maximum_viable_ceiling_z0, maximum_viable_ceiling_z1, delta_policy) ~ scenario + method,
  data = results_df,
  FUN = function(x) mean(x, na.rm = TRUE)
)

raw_path <- file.path(config$output_dir, "simulation_results.csv")
summary_path <- file.path(config$output_dir, "simulation_summary.csv")

write.csv(results_df, raw_path, row.names = FALSE)
write.csv(summary_df, summary_path, row.names = FALSE)

cat("Simulation run completed.\n")
cat("Raw results:", raw_path, "\n")
cat("Summary:", summary_path, "\n\n")
print(summary_df)
