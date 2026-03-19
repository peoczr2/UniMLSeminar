## Rank individuals by score and convert the top `share` fraction into a 0/1 target indicator.
target_top_share <- function(scores, share) {
  n_target <- max(1L, floor(length(scores) * share))
  rank_idx <- order(scores, decreasing = TRUE)
  targeted <- integer(length(scores))
  targeted[rank_idx[seq_len(n_target)]] <- 1L
  targeted
}

## IPS-style policy value:
## estimate the mean outcome contributed by units that were both targeted and actually treated.
policy_value <- function(targeted, y, w, treatment_probability = 0.5) {
  mean(ifelse(targeted == 1L & w == 1L, y / treatment_probability, 0), na.rm = TRUE)
}

## Paper-style efficiency:
## compare the policy value of the learned ranking against a random ranking,
## then scale to the "among the targeted" interpretation used in the tables.
relative_efficiency <- function(targeted, y, w) {
  random_targeted <- integer(length(targeted))
  random_targeted[sample.int(length(targeted), sum(targeted))] <- 1L
  target_share <- mean(targeted)
  baseline <- policy_value(random_targeted, y, w)

  if (isTRUE(all.equal(target_share, 0))) {
    return(NA_real_)
  }

  targeted_policy_value <- policy_value(targeted, y, w) / target_share
  targeted_baseline <- baseline / target_share

  100 * (targeted_policy_value - targeted_baseline)
}

## Euclidean distance between the average protected attributes of targeted and non-targeted units.
## In this project we pass `Z1` only, so the metric reduces to the absolute difference in means.
imbalance_metric <- function(targeted, z_matrix) {
  targeted_idx <- targeted == 1L
  non_targeted_idx <- targeted == 0L

  if (!any(targeted_idx) || !any(non_targeted_idx)) {
    return(NA_real_)
  }

  z_targeted <- colMeans(z_matrix[targeted_idx, , drop = FALSE])
  z_non_targeted <- colMeans(z_matrix[non_targeted_idx, , drop = FALSE])
  sqrt(sum((z_targeted - z_non_targeted) ^ 2))
}

## Demographic composition of the targeted group only.
## With binary `Z1`, this is the share of targeted units that belong to the protected group.
targeted_group_demographic <- function(targeted, z_matrix) {
  targeted_idx <- targeted == 1L

  if (!any(targeted_idx)) {
    return(NA_real_)
  }

  mean(z_matrix[targeted_idx, 1L])
}

## Maximum viable ceiling:
## calculates the absolute maximum representation each group can achieve
## inside a fixed budget, without treating negative-CATE units.
maximum_viable_ceiling <- function(test_data, budget_fraction = 0.50) {
  viable_idx <- test_data$tau > 0

  if (!any(viable_idx)) {
    return(c(Z1_0 = NA_real_, Z1_1 = NA_real_))
  }

  total_population <- length(test_data$tau)
  budget_size <- floor(total_population * budget_fraction)

  if (budget_size <= 0) {
    return(c(Z1_0 = NA_real_, Z1_1 = NA_real_))
  }

  viable_group_0_count <- sum(test_data$Z[viable_idx, "Z1"] == 0, na.rm = TRUE)
  viable_group_1_count <- sum(test_data$Z[viable_idx, "Z1"] == 1, na.rm = TRUE)

  c(
    Z1_0 = min(1.0, viable_group_0_count / budget_size),
    Z1_1 = min(1.0, viable_group_1_count / budget_size)
  )
}

extract_ceiling_value <- function(ceiling_values, group_label) {
  value <- ceiling_values[[group_label]]
  if (is.null(value) || length(value) == 0) {
    return(NA_real_)
  }
  as.numeric(value)
}

## Scenario-specific protected-attribute input for imbalance:
## use `Z1` only to mirror the paper-style table comparisons.
scenario_imbalance_input <- function(test_data) {
  as.matrix(test_data$Z[, "Z1", drop = FALSE])
}

## Individual-fairness proxy:
## percentage of units whose allocation changes when protected attributes are perturbed.
delta_policy_metric <- function(targeted, targeted_twin) {
  100 * mean(targeted != targeted_twin)
}

## Bundle all metrics for one fitted method in one scenario.
## `CF-FD` imbalance is used as the normalization reference.
evaluate_method <- function(test_data, method_result, target_share, cf_fd_imbalance) {
  targeted <- target_top_share(method_result$score, target_share)
  targeted_twin <- target_top_share(method_result$twin_score, target_share)
  raw_imbalance <- imbalance_metric(targeted, scenario_imbalance_input(test_data))
  targeted_group_demo <- targeted_group_demographic(targeted, scenario_imbalance_input(test_data))
  ceiling_values <- maximum_viable_ceiling(test_data, budget_fraction = target_share)

  data.frame(
    method = method_result$method,
    efficiency = relative_efficiency(targeted, test_data$Y, test_data$W),
    raw_imbalance = raw_imbalance,
    imbalance = if (cf_fd_imbalance == 0) NA_real_ else 100 * raw_imbalance / cf_fd_imbalance,
    targeted_group_demo = targeted_group_demo,
    maximum_viable_ceiling_z0 = extract_ceiling_value(ceiling_values, "Z1_0"),
    maximum_viable_ceiling_z1 = extract_ceiling_value(ceiling_values, "Z1_1"),
    delta_policy = delta_policy_metric(targeted, targeted_twin),
    stringsAsFactors = FALSE
  )
}

## Keep the raw score table so the scenario scripts can plot score distributions by method.
build_scores_frame <- function(test_data, method_results) {
  do.call(
    rbind,
    lapply(method_results, function(result) {
      data.frame(
        method = result$method,
        score = result$score,
        Z1 = test_data$Z[, "Z1"],
        tau = test_data$tau,
        stringsAsFactors = FALSE
      )
    })
  )
}
