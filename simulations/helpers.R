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

parse_method_selection <- function(methods) {
  if (is.null(methods) || length(methods) == 0L) {
    return(NULL)
  }

  methods <- trimws(as.character(methods))
  methods <- methods[nzchar(methods)]
  if (length(methods) == 0L || any(tolower(methods) == "all")) {
    return(NULL)
  }

  unique(trimws(unlist(strsplit(paste(methods, collapse = ","), ",", fixed = TRUE), use.names = FALSE)))
}

sanitize_cache_component <- function(value) {
  value <- as.character(value)
  value <- gsub("[^A-Za-z0-9._-]+", "_", value)
  gsub("^_|_$", "", value)
}

simulation_method_specs <- function() {
  specs <- list(
    "CF-FD" = function(sim_data, num_trees, beat_penalty, seed) {
      fit_cf_fd(sim_data, num_trees, seed)
    },
    "CF-NP" = function(sim_data, num_trees, beat_penalty, seed) {
      fit_cf_np(sim_data, num_trees, seed)
    },
    "Debiased" = function(sim_data, num_trees, beat_penalty, seed) {
      fit_debiased(sim_data, num_trees, seed)
    },
    "BEAT" = function(sim_data, num_trees, beat_penalty, seed) {
      fit_beat(sim_data, num_trees, beat_penalty, seed)
    },
    "BTGQ" = function(sim_data, num_trees, beat_penalty, seed) {
      fit_btgq(sim_data, num_trees, seed, budget = 0.5, target_quota = 0.5, use_validation = FALSE, method_name = "BTGQ")
    },
    "BTGQ_VALID" = function(sim_data, num_trees, beat_penalty, seed) {
      fit_btgq(sim_data, num_trees, seed + 50, budget = 0.5, target_quota = 0.5, use_validation = TRUE, method_name = "BTGQ_VALID")
    },
    "BTGQ_DUMB" = function(sim_data, num_trees, beat_penalty, seed) {
      fit_btgq_dumb(sim_data, num_trees, seed + 60, budget = 0.5, target_quota = 0.5, use_validation = FALSE, method_name = "BTGQ_DUMB")
    }
  )

  specs
}

method_cache_path <- function(cache_dir, cache_key, method_name) {
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  file.path(
    cache_dir,
    paste0(
      sanitize_cache_component(cache_key),
      "__",
      sanitize_cache_component(method_name),
      ".rds"
    )
  )
}

collect_method_results <- function(
    sim_data,
    num_trees,
    beat_penalty,
    seed,
    cache_dir,
    cache_key,
    selected_methods = NULL,
    use_cache = TRUE) {
  specs <- simulation_method_specs()
  method_names <- names(specs)

  if (!is.null(selected_methods)) {
    unknown_methods <- setdiff(selected_methods, method_names)
    if (length(unknown_methods) > 0L) {
      stop(
        paste("Unknown method(s):", paste(unknown_methods, collapse = ", ")),
        call. = FALSE
      )
    }
    method_names <- unique(c("CF-FD", selected_methods))
  }

  method_results <- vector("list", length(method_names))

  for (i in seq_along(method_names)) {
    method_name <- method_names[[i]]
    cache_path <- method_cache_path(cache_dir, cache_key, method_name)

    if (use_cache && file.exists(cache_path)) {
      method_results[[i]] <- readRDS(cache_path)
      next
    }

    method_result <- specs[[method_name]](sim_data, num_trees, beat_penalty, seed)
    method_results[[i]] <- method_result
    if (use_cache) {
      saveRDS(method_result, cache_path)
    }
  }

  method_results
}

get_method_result <- function(method_results, method_name) {
  if (is.null(method_results) || length(method_results) == 0L) {
    return(NULL)
  }

  method_names <- vapply(method_results, function(x) x$method, character(1L))
  idx <- match(method_name, method_names)
  if (is.na(idx)) {
    return(NULL)
  }
  method_results[[idx]]
}

scenario_plot_dir <- function(output_dir, scenario_name) {
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  if (basename(output_dir) == scenario_name) {
    return(output_dir)
  }
  file.path(output_dir, scenario_name)
}

load_simulation_runtime <- function(repo_root) {
  source(file.path(repo_root, "simulations", "helpers.R"))
  source(file.path(repo_root, "simulations", "metrics.R"))
  source(file.path(repo_root, "simulations", "method_cf_fd.R"))
  source(file.path(repo_root, "simulations", "method_cf_np.R"))
  source(file.path(repo_root, "simulations", "method_debiased.R"))
  source(file.path(repo_root, "simulations", "method_beat.R"))
  source(file.path(repo_root, "simulations", "method_btgq.R"))
  source(file.path(repo_root, "simulations", "method_btgq_dumb.R"))
}

default_scenario_config <- function(repo_root, scenario_name) {
  parse_cli_args(list(
    n_train = 2000,
    n_test = 1000,
    num_trees = 500,
    target_share = 0.5,
    beat_penalty = 10,
    seed = 123,
    methods = "all",
    use_cache = TRUE,
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
  selected_methods <- parse_method_selection(config$methods)

  metrics_df <- run_fun(
    repo_root = repo_root,
    output_dir = config$output_dir,
    n_train = as.integer(config$n_train),
    n_test = as.integer(config$n_test),
    num_trees = as.integer(config$num_trees),
    beat_penalty = config$beat_penalty,
    target_share = config$target_share,
    seed = as.integer(config$seed),
    selected_methods = selected_methods,
    use_cache = isTRUE(config$use_cache)
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

plot_btgq_lambda_trace <- function(lambda_trace, output_path, title_text, target_quota = 0.5) {
  if (is.null(lambda_trace) || nrow(lambda_trace) == 0) {
    return(invisible(NULL))
  }

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 is not installed; skipping plot generation.", call. = FALSE)
    return(invisible(NULL))
  }

  plot_data <- transform(lambda_trace, lambda_plot = sign(lambda) * log10(abs(lambda) + 1))

  plot_obj <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = lambda_plot, y = targeted_group_demo, group = 1)
  ) +
    ggplot2::geom_point(color = "#1b6ca8", size = 2) +
    ggplot2::geom_hline(
      yintercept = target_quota,
      linetype = "dashed",
      color = "#c1121f"
    ) +
    ggplot2::labs(
      title = title_text,
      x = "signed log10(|lambda| + 1)",
      y = "Targeted Group Demo"
    ) +
    ggplot2::theme_minimal()

  if (nrow(lambda_trace) >= 2) {
    plot_obj <- plot_obj + ggplot2::geom_line(color = "#1b6ca8")
  }

  ggplot2::ggsave(output_path, plot_obj, width = 8, height = 5, dpi = 150)
  invisible(plot_obj)
}
