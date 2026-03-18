args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args[grep(file_arg, args)])

if (length(script_path) != 1 || !nzchar(script_path)) {
  stop("Run this script with Rscript scripts/reinstall_beat.R", call. = FALSE)
}

script_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = TRUE))
repo_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
lib_dir <- file.path(repo_root, ".r-library")
beat_dir <- file.path(repo_root, "vendor", "beat")

if (!dir.exists(lib_dir)) {
  stop("Local R library not found. Run Rscript scripts/setup_beat.R first.", call. = FALSE)
}

if (!dir.exists(beat_dir)) {
  stop("Expected vendored BEAT source at vendor/beat", call. = FALSE)
}

.libPaths(c(lib_dir, .libPaths()))

if (!requireNamespace("remotes", quietly = TRUE)) {
  stop("Package 'remotes' is required in .r-library. Run Rscript scripts/setup_beat.R first.", call. = FALSE)
}

remotes::install_local(
  beat_dir,
  upgrade = "never",
  dependencies = FALSE,
  force = TRUE,
  lib = lib_dir
)

message("BEAT reinstalled into: ", lib_dir)
message("Next step: Rscript scripts/run_beat_smoke_test.R")
