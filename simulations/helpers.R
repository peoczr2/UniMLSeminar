get_repo_root_from_script <- function(script_path) {
  normalizePath(
    file.path(dirname(normalizePath(script_path, winslash = "/", mustWork = TRUE)), ".."),
    winslash = "/",
    mustWork = TRUE
  )
}

ensure_local_library <- function(repo_root) {
  lib_dir <- file.path(repo_root, ".r-library")
  if (!dir.exists(lib_dir)) {
    stop("Local R library not found. Run Rscript scripts/setup_beat.R first.", call. = FALSE)
  }

  .libPaths(c(lib_dir, .libPaths()))
  invisible(lib_dir)
}

parse_cli_args <- function(defaults) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0) {
    return(defaults)
  }

  parsed <- defaults
  for (arg in args) {
    if (!startsWith(arg, "--") || !grepl("=", arg, fixed = TRUE, useBytes = TRUE)) {
      next
    }

    key_value <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    if (length(key_value) != 2) {
      next
    }

    key <- gsub("-", "_", key_value[1], fixed = TRUE)
    value <- key_value[2]
    if (!key %in% names(parsed)) {
      next
    }

    template <- parsed[[key]]
    if (is.numeric(template)) {
      parsed[[key]] <- as.numeric(value)
    } else if (is.logical(template)) {
      parsed[[key]] <- tolower(value) %in% c("1", "true", "t", "yes", "y")
    } else {
      parsed[[key]] <- value
    }
  }

  parsed
}

load_simulation_runtime <- function(repo_root) {
  source(file.path(repo_root, "simulations", "helpers.R"))
  source(file.path(repo_root, "simulations", "metrics.R"))
  source(file.path(repo_root, "simulations", "method_cf_fd.R"))
  source(file.path(repo_root, "simulations", "method_cf_np.R"))
  source(file.path(repo_root, "simulations", "method_debiased.R"))
  source(file.path(repo_root, "simulations", "method_beat.R"))
}

default_scenario_config <- function(repo_root, scenario_name) {
  parse_cli_args(list(
    n_train = 2000,
    n_test = 1000,
    num_trees = 500,
    target_share = 0.5,
    beat_penalty = 10,
    seed = 123,
    output_dir = file.path(repo_root, "outputs", "simulations", scenario_name)
  ))
}

run_scenario_script <- function(run_fun, scenario_name) {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script_path <- sub(file_arg, "", args[grep(file_arg, args)])

  if (length(script_path) != 1 || !nzchar(script_path)) {
    stop(paste("Run this script with Rscript simulations", paste0("/", scenario_name, ".R")), call. = FALSE)
  }

  repo_root <- get_repo_root_from_script(script_path)
  load_simulation_runtime(repo_root)
  ensure_local_library(repo_root)
  suppressPackageStartupMessages(library(beat))

  config <- default_scenario_config(repo_root, scenario_name)
  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)

  metrics_df <- run_fun(
    repo_root = repo_root,
    output_dir = config$output_dir,
    n_train = as.integer(config$n_train),
    n_test = as.integer(config$n_test),
    num_trees = as.integer(config$num_trees),
    beat_penalty = config$beat_penalty,
    target_share = config$target_share,
    seed = as.integer(config$seed)
  )

  output_path <- file.path(config$output_dir, paste0(scenario_name, "_metrics.csv"))
  write.csv(metrics_df, output_path, row.names = FALSE)

  cat("Scenario run completed.\n")
  cat("Metrics:", output_path, "\n\n")
  print(metrics_df)
}

make_base_features <- function(n, corr_strength) {
  x_cont <- matrix(rnorm(n * 4), nrow = n, ncol = 4)
  x_disc <- matrix(rbinom(n * 6, 1, 0.3), nrow = n, ncol = 6)
  x <- cbind(x_cont, x_disc)
  colnames(x) <- paste0("X", seq_len(ncol(x)))

  z1_prob <- plogis(corr_strength * x_cont[, 2])
  z1 <- rbinom(n, 1, z1_prob)
  z2 <- rnorm(n)
  z3 <- rnorm(n)
  z4 <- rbinom(n, 1, 0.3)
  z <- cbind(Z1 = z1, Z2 = z2, Z3 = z3, Z4 = z4)

  list(X = x, Z = z)
}

make_standard_train_test <- function(x, z, tau, n_train, n_test) {
  n <- n_train + n_test
  w <- rbinom(n, 1, 0.5)
  y <- tau * w + rnorm(n)

  train_idx <- seq_len(n_train)
  test_idx <- seq.int(n_train + 1, n)

  list(
    train = list(
      X = x[train_idx, , drop = FALSE],
      Z = z[train_idx, , drop = FALSE],
      X_full = cbind(x[train_idx, , drop = FALSE], z[train_idx, , drop = FALSE]),
      W = w[train_idx],
      Y = y[train_idx],
      tau = tau[train_idx]
    ),
    test = list(
      X = x[test_idx, , drop = FALSE],
      Z = z[test_idx, , drop = FALSE],
      X_full = cbind(x[test_idx, , drop = FALSE], z[test_idx, , drop = FALSE]),
      W = w[test_idx],
      Y = y[test_idx],
      tau = tau[test_idx]
    )
  )
}

flip_protected <- function(z) {
  z_flipped <- z
  for (j in seq_len(ncol(z_flipped))) {
    unique_values <- sort(unique(z_flipped[, j]))
    if (length(unique_values) == 2 && identical(unique_values, c(0, 1))) {
      z_flipped[, j] <- 1 - z_flipped[, j]
    } else {
      z_mean <- mean(z_flipped[, j], na.rm = TRUE)
      z_sd <- stats::sd(z_flipped[, j], na.rm = TRUE)

      if (isTRUE(all.equal(z_sd, 0))) {
        next
      }

      z_std <- (z_flipped[, j] - z_mean) / z_sd
      z_std_shifted <- ifelse(z_std >= 0, z_std - 1, z_std + 1)
      z_flipped[, j] <- z_std_shifted * z_sd + z_mean
    }
  }
  z_flipped
}

plot_method_scores <- function(scores_df, output_path, title_text) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 is not installed; skipping plot generation.", call. = FALSE)
    return(invisible(NULL))
  }

  plot_obj <- ggplot2::ggplot(
    scores_df,
    ggplot2::aes(x = score, color = factor(Z1), fill = factor(Z1))
  ) +
    ggplot2::geom_density(alpha = 0.2) +
    ggplot2::facet_wrap(~method, scales = "free") +
    ggplot2::labs(
      title = title_text,
      x = "Predicted score",
      y = "Density",
      color = "Z1",
      fill = "Z1"
    ) +
    ggplot2::theme_minimal()

  ggplot2::ggsave(output_path, plot_obj, width = 12, height = 8, dpi = 150)
  invisible(plot_obj)
}
