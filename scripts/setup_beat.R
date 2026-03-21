args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args[grep(file_arg, args)])

if (length(script_path) != 1 || !nzchar(script_path)) {
  stop("Run this script with Rscript scripts/setup_beat.R", call. = FALSE)
}

script_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = TRUE))
repo_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
lib_dir <- file.path(repo_root, ".r-library")
beat_dir <- file.path(repo_root, "vendor", "beat")

if (!dir.exists(beat_dir)) {
  stop("Expected vendored BEAT source at vendor/beat", call. = FALSE)
}

if (getRversion() > "4.2.3") {
  stop(
    paste(
      "BEAT upstream documents support only for R 4.2.3 or earlier.",
      "Install R 4.2.3 before running this setup."
    ),
    call. = FALSE
  )
}

dir.create(lib_dir, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib_dir, .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"))

install_if_missing <- function(packages) {
  installed <- rownames(installed.packages(lib.loc = .libPaths()))
  missing <- setdiff(packages, installed)
  if (length(missing) > 0) {
    install.packages(missing, lib = lib_dir)
  }
}

package_version_in_paths <- function(package) {
  for (lib_path in .libPaths()) {
    if (requireNamespace(package, quietly = TRUE, lib.loc = lib_path)) {
      return(as.character(utils::packageVersion(package, lib.loc = lib_path)))
    }
  }

  NULL
}

install_pinned_if_needed <- function(package, version) {
  installed_version <- package_version_in_paths(package)

  if (!is.null(installed_version) && identical(installed_version, version)) {
    message(
      "Using existing ", package, " ", version,
      " from local library path."
    )
    return(invisible(FALSE))
  }

  if (!is.null(installed_version)) {
    message(
      package, " ", installed_version,
      " is installed, but BEAT requires ", version,
      ". Reinstalling pinned version."
    )
  }

  remotes::install_version(
    package,
    version = version,
    upgrade = "never",
    dependencies = FALSE,
    lib = lib_dir
  )

  invisible(TRUE)
}

install_if_missing(c("remotes"))

install_pinned_if_needed("RcppEigen", "0.3.3.7.0")
install_pinned_if_needed("RcppArmadillo", "0.11.4.0.1")

install_if_missing(c(
  "Rcpp",
  "data.table",
  "DiceKriging",
  "lmtest",
  "Matrix",
  "Rcpp",
  "sandwich",
  "doParallel",
  "foreach",
  "iterators",
  "testthat",
  "arrangements",
  "ROI",
  "stringi"
))

remotes::install_local(
  beat_dir,
  upgrade = "never",
  dependencies = FALSE,
  force = TRUE,
  lib = lib_dir
)

message("BEAT installed into: ", lib_dir)
message("Next step: Rscript scripts/run_beat_smoke_test.R")
