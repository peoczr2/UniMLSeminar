#' BTGQ Max Causal Forest
#'
#' BTGQ-Max keeps the BTGQ local split rule and adds a global selection loop
#' that shrinks the effective budget to the viable positive-prediction pool.
#'
#' @param budget.max Maximum administrative budget ceiling used in policy
#'   evaluation.
#' @param target.quota Minimum target-group share required inside the selected
#'   budget slice.
#' @param validation.X Optional validation set used for lambda search.
#' @param validation.Y Optional validation outcomes used for the efficiency
#'   objective.
#' @param validation.W Optional validation treatments used for the efficiency
#'   objective.
#' @param validation.target.group Optional validation target-group indicator.
#' @param lambda.lower Lower bound for the lambda search.
#' @param lambda.upper Optional upper bound for the lambda search.
#' @param viability.epsilon Minimum positive prediction treated as viable.
#' @param max.search.iter Maximum number of coarse/refine search iterations.
#' @inheritParams btgq_causal_forest
#'
#' @return A trained BTGQ-Max causal forest object.
#' @export
btgq_max_causal_forest <- function(X,
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

  forest <- do.call.rcpp(btgq_max_causal_train, c(data, args))

  class(forest) <- c("btgq_max_causal_forest", "btgq_causal_forest", "causal_forest", "grf")
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

#' Tune BTGQ-Max lambda using a budget ceiling and quota target
#'
#' @inheritParams btgq_causal_forest
#' @param budget.max Maximum administrative budget ceiling used for the policy
#'   slice. The actual budget may shrink below this ceiling if the viable pool is
#'   smaller.
#' @param target.quota Minimum target-group share required inside the selected
#'   slice.
#' @param validation.X Optional data set used to evaluate the quota during the
#'   search. If NULL, the training set is used as before.
#' @param validation.target.group Optional target-group indicator for
#'   `validation.X`. Required when `validation.X` is supplied.
#' @param lambda.lower Lower bound for the BTGQ-Max lambda search.
#' @param lambda.upper Optional upper bound for the BTGQ-Max lambda search. If
#'   NULL, it is expanded automatically by doubling.
#' @param viability.epsilon Small positive threshold used to define the viable
#'   pool.
#' @param max.search.iter Maximum number of doubling/binary-search iterations.
#'
#' @return A list containing the selected forest, search trace, and quota metrics.
#' @export
tune_btgq_max_causal_forest <- function(X,
                                        Y,
                                        W,
                                        target.group,
                                        budget.max,
                                        target.quota,
                                        validation.X = NULL,
                                        validation.Y = NULL,
                                        validation.W = NULL,
                                        validation.target.group = NULL,
                                        lambda.lower = -1e12,
                                        lambda.upper = NULL,
                                        viability.epsilon = 1e-08,
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
  validate_btgq_max_budget_quota(budget.max, target.quota, viability.epsilon, max.search.iter)
  target.group <- validate_btgq_target_group(target.group, X)

  if (is.null(validation.X)) {
    validation.X <- X
    validation.target.group <- target.group
    validation.Y <- Y
    validation.W <- W
  } else {
    validate_X(validation.X, allow.na = TRUE)
    validation.target.group <- validate_btgq_target_group(validation.target.group, validation.X)
    if (is.null(validation.Y) || length(validation.Y) != nrow(validation.X)) {
      stop("validation.Y must be provided and match validation.X when validation.X is supplied.")
    }
    if (is.null(validation.W) || length(validation.W) != nrow(validation.X)) {
      stop("validation.W must be provided and match validation.X when validation.X is supplied.")
    }
  }

  fit_once <- function(lambda, Y.hat.arg, W.hat.arg) {
    btgq_max_causal_forest(
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
    policy.value = numeric(0),
    efficiency = numeric(0),
    viable.share = numeric(0),
    budget.actual = numeric(0),
    quota.gap = numeric(0),
    phase = character(0),
    stringsAsFactors = FALSE
  )
  fitted.models <- list()

  evaluate_fit <- function(forest) {
    validation.predictions <- predict(forest, validation.X)$predictions
    stats <- btgq_max_policy_stats(validation.predictions, validation.target.group, budget.max, viability.epsilon)
    targeted <- integer(length(validation.predictions))
    if (length(stats$selected) > 0) {
      targeted[stats$selected] <- 1L
    }
    list(
      quota = stats$quota,
      policy.value = stats$policy.value,
      efficiency = btgq_max_ipsefficiency(targeted, validation.Y, validation.W),
      viable.share = stats$viable.share,
      budget.actual = stats$budget.actual,
      stats = stats
    )
  }

  record_fit <- function(forest, phase) {
    eval <- evaluate_fit(forest)
    key <- formatC(forest$btgq.lambda, digits = 17, format = "fg", flag = "#")
    fitted.models[[key]] <<- list(
      forest = forest,
      quota = eval$quota,
      policy.value = eval$policy.value,
      efficiency = eval$efficiency,
      viable.share = eval$viable.share,
      budget.actual = eval$budget.actual,
      stats = eval$stats
    )
    search.trace <<- rbind(
      search.trace,
      data.frame(
        lambda = forest$btgq.lambda,
        achieved.quota = eval$quota,
        policy.value = eval$policy.value,
        efficiency = eval$efficiency,
        viable.share = eval$viable.share,
        budget.actual = eval$budget.actual,
        quota.gap = eval$quota - target.quota,
        phase = phase,
        stringsAsFactors = FALSE
      )
    )
    eval
  }

  baseline <- fit_once(0, Y.hat, W.hat)
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
  coarse_lambdas <- unique(c(0, sort(coarse_magnitudes), -sort(coarse_magnitudes)))
  lower_bound <- if (is.null(lambda.lower)) -Inf else lambda.lower
  upper_bound <- if (is.null(lambda.upper)) Inf else lambda.upper
  coarse_lambdas <- coarse_lambdas[coarse_lambdas >= lower_bound & coarse_lambdas <= upper_bound]

  for (lambda in coarse_lambdas) {
    evaluate_lambda(lambda, "coarse")
  }

  coarse.trace <- search.trace[search.trace$phase == "coarse", , drop = FALSE]
  if (nrow(coarse.trace) == 0) {
    evaluate_lambda(lambda.lower, "coarse")
    coarse.trace <- search.trace[search.trace$phase == "coarse", , drop = FALSE]
  }

  refine.centers <- coarse.trace$lambda[order(-coarse.trace$achieved.quota, -coarse.trace$policy.value, coarse.trace$lambda)]
  refine.centers <- unique(refine.centers[seq_len(min(3L, length(refine.centers)))])

  make_refine_grid <- function(center) {
    if (center == 0) {
      grid <- c(0, -0.01, -0.03, -0.1, -0.3, -1, 0.01, 0.03, 0.1, 0.3, 1)
    } else {
      grid <- center * c(0.25, 0.5, 0.75, 0.85, 1, 1.15, 1.3, 1.5, 2)
    }
    grid <- unique(grid)
    grid <- grid[grid >= lower_bound & grid <= upper_bound]
    grid
  }

  for (center in refine.centers) {
    refine.grid <- make_refine_grid(center)
    for (lambda in refine.grid) {
      evaluate_lambda(lambda, "refine")
    }
  }

  best <- btgq_max_best_search_result(fitted.models, target.quota)
  list(
    forest = best$forest,
    selected.lambda = best$forest$btgq.lambda,
    achieved.quota = best$quota,
    policy.value = best$policy.value,
    viable.share = best$viable.share,
    budget.actual = best$budget.actual,
    target.quota = target.quota,
    budget.max = budget.max,
    quota.attained = best$quota >= target.quota,
    search.trace = search.trace
  )
}

validate_btgq_max_budget_quota <- function(budget.max, target.quota, viability.epsilon, max.search.iter) {
  if (!is.numeric(budget.max) || length(budget.max) != 1 || budget.max <= 0 || budget.max > 1) {
    stop("budget.max must be a scalar in (0, 1].")
  }
  if (!is.numeric(target.quota) || length(target.quota) != 1 || target.quota < 0 || target.quota > 1) {
    stop("target.quota must be a scalar in [0, 1].")
  }
  if (!is.numeric(viability.epsilon) || length(viability.epsilon) != 1 || viability.epsilon < 0) {
    stop("viability.epsilon must be a non-negative scalar.")
  }
  if (!is.numeric(max.search.iter) || length(max.search.iter) != 1 || max.search.iter < 1) {
    stop("max.search.iter must be a positive scalar.")
  }
}

btgq_max_policy_stats <- function(predictions, target.group, budget.max, viability.epsilon) {
  n <- length(predictions)
  if (n == 0) {
    return(list(
      selected = integer(0),
      quota = 0,
      policy.value = 0,
      viable.share = 0,
      budget.actual = 0,
      threshold = NA_real_
    ))
  }

  rank.order <- order(-predictions, seq_along(predictions))
  viable.n <- sum(predictions > viability.epsilon)
  ceiling.n <- max(0L, ceiling(budget.max * n))
  selected.n <- min(ceiling.n, viable.n)
  selected <- if (selected.n > 0) rank.order[seq_len(selected.n)] else integer(0)

  list(
    selected = selected,
    quota = if (selected.n > 0) mean(target.group[selected]) else 0,
    policy.value = if (selected.n > 0) sum(predictions[selected]) else 0,
    viable.share = viable.n / n,
    budget.actual = selected.n / n,
    threshold = if (selected.n > 0) predictions[selected[selected.n]] else NA_real_
  )
}

btgq_max_ipsefficiency <- function(targeted, y, w, treatment_probability = 0.5) {
  if (isTRUE(all.equal(mean(targeted), 0))) {
    return(NA_real_)
  }

  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv)
  } else {
    old_seed <- NULL
  }
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(1)
  random_targeted <- integer(length(targeted))
  random_targeted[sample.int(length(targeted), sum(targeted))] <- 1L
  target_share <- mean(targeted)
  if (isTRUE(all.equal(target_share, 0))) {
    return(NA_real_)
  }

  policy_value <- function(targeted_vec, y_vec, w_vec, treatment_probability = 0.5) {
    mean(ifelse(targeted_vec == 1L & w_vec == 1L, y_vec / treatment_probability, 0), na.rm = TRUE)
  }

  targeted_policy_value <- policy_value(targeted, y, w, treatment_probability) / target_share
  targeted_baseline <- policy_value(random_targeted, y, w, treatment_probability) / target_share
  100 * (targeted_policy_value - targeted_baseline)
}

btgq_max_best_search_result <- function(fitted.models, target.quota) {
  model.entries <- unname(fitted.models)
  quotas <- vapply(model.entries, function(entry) entry$quota, numeric(1))
  values <- vapply(model.entries, function(entry) entry$efficiency, numeric(1))
  budgets <- vapply(model.entries, function(entry) entry$budget.actual, numeric(1))
  lambdas <- vapply(model.entries, function(entry) entry$forest$btgq.lambda, numeric(1))
  values[is.na(values)] <- -Inf

  hit_idx <- which(quotas >= target.quota)
  if (length(hit_idx) > 0) {
    best.index <- hit_idx[order(-values[hit_idx], -quotas[hit_idx], budgets[hit_idx], lambdas[hit_idx])[1]]
    return(model.entries[[best.index]])
  }

  best.index <- order(-quotas, -values, budgets, lambdas)[1]
  model.entries[[best.index]]
}
