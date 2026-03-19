#' BTGQ Dumb Causal Forest
#'
#' Trains a causal forest using a simplified BTGQ split criterion that rewards
#' target-group concentration without multiplying by the estimated treatment effect.
#'
#' @param X The covariates used in the causal regression.
#' @param Y The outcome.
#' @param W The treatment assignment.
#' @param target.group Binary indicator for the favored group. Must be a numeric or
#'   logical vector of length `nrow(X)` with values in `{0, 1}`.
#' @param btgq.lambda Reward weight used by the BTGQ split rule.
#' @param Y.hat Optional outcome nuisance estimates. If NULL they are estimated.
#' @param W.hat Optional treatment nuisance estimates. If NULL they are estimated.
#' @inheritParams balanced_causal_forest
#'
#' @return A trained BTGQ dumb causal forest object.
#'
#' @import data.table
#' @import Rcpp
#' @export
btgq_dumb_causal_forest <- function(X,
                                    Y,
                                    W,
                                    target.group,
                                    btgq.lambda = 0,
                                    Y.hat = NULL,
                                    W.hat = NULL,
                                    num.trees = 2000,
                                    sample.weights = NULL,
                                    clusters = NULL,
                                    equalize.cluster.weights = FALSE,
                                    sample.fraction = 0.5,
                                    mtry = min(ceiling(sqrt(ncol(X)) + 20), ncol(X)),
                                    min.node.size = 5,
                                    honesty = TRUE,
                                    honesty.fraction = 0.5,
                                    honesty.prune.leaves = TRUE,
                                    alpha = 0.05,
                                    imbalance.penalty = 0,
                                    stabilize.splits = TRUE,
                                    ci.group.size = 2,
                                    compute.oob.predictions = TRUE,
                                    orthog.boosting = FALSE,
                                    num.threads = NULL,
                                    seed = runif(1, 0, .Machine$integer.max)) {
  has.missing.values <- validate_X(X, allow.na = TRUE)
  validate_sample_weights(sample.weights, X)
  Y <- validate_observations(Y, X)
  W <- validate_observations(W, X)
  target.group <- validate_btgq_target_group(target.group, X)
  stopifnot(is.numeric(btgq.lambda), length(btgq.lambda) == 1)

  clusters <- validate_clusters(clusters, X)
  samples.per.cluster <- validate_equalize_cluster_weights(equalize.cluster.weights, clusters, sample.weights)
  num.threads <- validate_num_threads(num.threads)

  args.orthog <- list(
    X = X,
    num.trees = max(50, num.trees / 4),
    sample.weights = sample.weights,
    clusters = clusters,
    equalize.cluster.weights = equalize.cluster.weights,
    sample.fraction = sample.fraction,
    mtry = mtry,
    min.node.size = 5,
    honesty = TRUE,
    honesty.fraction = 0.5,
    honesty.prune.leaves = honesty.prune.leaves,
    alpha = alpha,
    imbalance.penalty = imbalance.penalty,
    ci.group.size = 1,
    seed = seed
  )

  if (is.null(Y.hat) && !orthog.boosting) {
    forest.Y <- do.call(regression_forest, c(Y = list(Y), args.orthog))
    Y.hat <- predict(forest.Y)$predictions
  } else if (is.null(Y.hat) && orthog.boosting) {
    forest.Y <- do.call(boosted_regression_forest, c(Y = list(Y), args.orthog))
    Y.hat <- predict(forest.Y)$predictions
  } else if (length(Y.hat) == 1) {
    Y.hat <- rep(Y.hat, nrow(X))
  } else if (length(Y.hat) != nrow(X)) {
    stop("Y.hat has incorrect length.")
  }

  if (is.null(W.hat) && !orthog.boosting) {
    forest.W <- do.call(regression_forest, c(Y = list(W), args.orthog))
    W.hat <- predict(forest.W)$predictions
  } else if (is.null(W.hat) && orthog.boosting) {
    forest.W <- do.call(boosted_regression_forest, c(Y = list(W), args.orthog))
    W.hat <- predict(forest.W)$predictions
  } else if (length(W.hat) == 1) {
    W.hat <- rep(W.hat, nrow(X))
  } else if (length(W.hat) != nrow(X)) {
    stop("W.hat has incorrect length.")
  }

  Y.centered <- Y - Y.hat
  W.centered <- W - W.hat
  data <- create_train_matrices(X, outcome = Y.centered, treatment = W.centered, sample.weights = sample.weights)

  args <- list(
    num.trees = num.trees,
    btgq.target.group = as.numeric(target.group),
    btgq.lambda = btgq.lambda,
    clusters = clusters,
    samples.per.cluster = samples.per.cluster,
    sample.fraction = sample.fraction,
    mtry = mtry,
    min.node.size = min.node.size,
    honesty = honesty,
    honesty.fraction = honesty.fraction,
    honesty.prune.leaves = honesty.prune.leaves,
    alpha = alpha,
    imbalance.penalty = imbalance.penalty,
    stabilize.splits = stabilize.splits,
    ci.group.size = ci.group.size,
    compute.oob.predictions = compute.oob.predictions,
    num.threads = num.threads,
    seed = seed,
    reduced.form.weight = 0
  )

  forest <- do.call.rcpp(btgq_dumb_causal_train, c(data, args))

  class(forest) <- c("btgq_dumb_causal_forest", "btgq_causal_forest", "causal_forest", "grf")
  forest[["ci.group.size"]] <- ci.group.size
  forest[["X.orig"]] <- X
  forest[["Y.orig"]] <- Y
  forest[["W.orig"]] <- W
  forest[["Y.hat"]] <- Y.hat
  forest[["W.hat"]] <- W.hat
  forest[["target.group"]] <- as.numeric(target.group)
  forest[["btgq.lambda"]] <- btgq.lambda
  forest[["clusters"]] <- clusters
  forest[["equalize.cluster.weights"]] <- equalize.cluster.weights
  forest[["sample.weights"]] <- sample.weights
  forest[["has.missing.values"]] <- has.missing.values

  forest
}

#' Tune BTGQ dumb lambda using a budget and quota target
#'
#' @param budget Budget fraction used for the top-k policy evaluation.
#' @param target.quota Desired share of the target group within the budgeted slice.
#' @param validation.X Optional data set used to evaluate the quota during the
#'   search. If NULL, the training set is used as before.
#' @param validation.target.group Optional target-group indicator for
#'   `validation.X`. Required when `validation.X` is supplied.
#' @param lambda.lower Lower bound for the BTGQ lambda search.
#' @param lambda.upper Optional upper bound for the BTGQ lambda search. If NULL, it
#'   is expanded automatically by doubling.
#' @param quota.tol Acceptable absolute difference between achieved and target quota.
#' @param max.search.iter Maximum number of doubling/binary-search iterations.
#' @inheritParams btgq_dumb_causal_forest
#'
#' @return A list containing the selected forest, search trace, and quota metrics.
#' @export
tune_btgq_dumb_causal_forest <- function(X,
                                         Y,
                                         W,
                                         target.group,
                                         budget,
                                         target.quota,
                                         validation.X = NULL,
                                         validation.target.group = NULL,
                                         lambda.lower = 0,
                                         lambda.upper = NULL,
                                         quota.tol = 0.01,
                                         max.search.iter = 12,
                                         Y.hat = NULL,
                                         W.hat = NULL,
                                         num.trees = 2000,
                                         sample.weights = NULL,
                                         clusters = NULL,
                                         equalize.cluster.weights = FALSE,
                                         sample.fraction = 0.5,
                                         mtry = min(ceiling(sqrt(ncol(X)) + 20), ncol(X)),
                                         min.node.size = 5,
                                         honesty = TRUE,
                                         honesty.fraction = 0.5,
                                         honesty.prune.leaves = TRUE,
                                         alpha = 0.05,
                                         imbalance.penalty = 0,
                                         stabilize.splits = TRUE,
                                         ci.group.size = 2,
                                         orthog.boosting = FALSE,
                                         num.threads = NULL,
                                         seed = runif(1, 0, .Machine$integer.max)) {
  validate_btgq_budget_quota(budget, target.quota, quota.tol, max.search.iter)
  target.group <- validate_btgq_target_group(target.group, X)

  if (is.null(validation.X)) {
    validation.X <- X
    validation.target.group <- target.group
  } else {
    validate_X(validation.X, allow.na = TRUE)
    validation.target.group <- validate_btgq_target_group(validation.target.group, validation.X)
  }

  fit_once <- function(lambda, Y.hat.arg, W.hat.arg) {
    btgq_dumb_causal_forest(
      X = X,
      Y = Y,
      W = W,
      target.group = target.group,
      btgq.lambda = lambda,
      Y.hat = Y.hat.arg,
      W.hat = W.hat.arg,
      num.trees = num.trees,
      sample.weights = sample.weights,
      clusters = clusters,
      equalize.cluster.weights = equalize.cluster.weights,
      sample.fraction = sample.fraction,
      mtry = mtry,
      min.node.size = min.node.size,
      honesty = honesty,
      honesty.fraction = honesty.fraction,
      honesty.prune.leaves = honesty.prune.leaves,
      alpha = alpha,
      imbalance.penalty = imbalance.penalty,
      stabilize.splits = stabilize.splits,
      ci.group.size = ci.group.size,
      compute.oob.predictions = TRUE,
      orthog.boosting = orthog.boosting,
      num.threads = num.threads,
      seed = seed
    )
  }

  search.trace <- data.frame(
    lambda = numeric(0),
    achieved.quota = numeric(0),
    quota.gap = numeric(0),
    phase = character(0),
    stringsAsFactors = FALSE
  )
  fitted.models <- list()

  evaluate_fit <- function(forest) {
    validation.predictions <- predict(forest, validation.X)$predictions
    stats <- btgq_policy_stats(validation.predictions, validation.target.group, budget)
    list(quota = stats$quota, stats = stats)
  }

  record_fit <- function(forest, phase) {
    eval <- evaluate_fit(forest)
    key <- formatC(forest$btgq.lambda, digits = 17, format = "fg", flag = "#")
    fitted.models[[key]] <<- list(forest = forest, quota = eval$quota, stats = eval$stats)
    search.trace <<- rbind(
      search.trace,
      data.frame(
        lambda = forest$btgq.lambda,
        achieved.quota = eval$quota,
        quota.gap = eval$quota - target.quota,
        phase = phase,
        stringsAsFactors = FALSE
      )
    )
    eval
  }

  baseline <- fit_once(lambda.lower, Y.hat, W.hat)
  baseline.eval <- record_fit(baseline, "baseline")
  nuisance.Y.hat <- baseline$Y.hat
  nuisance.W.hat <- baseline$W.hat

  lambda_key <- function(lambda) {
    formatC(lambda, digits = 17, format = "fg", flag = "#")
  }

  seen_lambda <- function(lambda) {
    lambda_key(lambda) %in% names(fitted.models)
  }

  evaluate_lambda <- function(lambda, phase) {
    if (seen_lambda(lambda)) {
      return(invisible(NULL))
    }
    fit <- fit_once(lambda, nuisance.Y.hat, nuisance.W.hat)
    record_fit(fit, phase)
  }

  coarse_magnitudes <- as.vector(outer(c(0.1, 0.3, 1, 3), 10^(0:max.search.iter)))
  coarse_lambdas <- unique(c(0, sort(c(-coarse_magnitudes, coarse_magnitudes))))
  if (!is.null(lambda.upper)) {
    coarse_lambdas <- coarse_lambdas[coarse_lambdas <= lambda.upper]
  }
  if (lambda.lower != 0) {
    coarse_lambdas <- unique(c(lambda.lower, coarse_lambdas))
  }

  for (lambda in coarse_lambdas) {
    evaluate_lambda(lambda, "coarse")
  }

  coarse.trace <- search.trace[search.trace$phase == "coarse", , drop = FALSE]
  if (nrow(coarse.trace) == 0) {
    evaluate_lambda(lambda.lower, "coarse")
    coarse.trace <- search.trace[search.trace$phase == "coarse", , drop = FALSE]
  }

  refine.centers <- coarse.trace$lambda[order(abs(coarse.trace$quota.gap), coarse.trace$lambda)]
  refine.centers <- unique(refine.centers[seq_len(min(3L, length(refine.centers)))])

  make_refine_grid <- function(center) {
    if (center == 0) {
      grid <- c(0, -0.01, -0.03, -0.1, -0.3, -1, 0.01, 0.03, 0.1, 0.3, 1)
    } else {
      grid <- center * c(0.25, 0.5, 0.75, 0.85, 1, 1.15, 1.3, 1.5, 2)
    }
    grid <- unique(grid)
    if (!is.null(lambda.upper)) {
      grid <- grid[grid <= lambda.upper]
    }
    grid
  }

  for (center in refine.centers) {
    refine.grid <- make_refine_grid(center)
    for (lambda in refine.grid) {
      evaluate_lambda(lambda, "refine")
    }
  }

  best <- btgq_best_search_result(fitted.models, target.quota)
  list(
    forest = best$forest,
    selected.lambda = best$forest$btgq.lambda,
    achieved.quota = best$quota,
    target.quota = target.quota,
    budget = budget,
    quota.attained = abs(best$quota - target.quota) <= quota.tol,
    search.trace = search.trace
  )
}
