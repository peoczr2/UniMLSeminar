target_top_share <- function(scores, share) {
  n_target <- max(1L, floor(length(scores) * share))
  rank_idx <- order(scores, decreasing = TRUE)
  targeted <- integer(length(scores))
  targeted[rank_idx[seq_len(n_target)]] <- 1L
  targeted
}

policy_value <- function(targeted, y, w, treatment_probability = 0.5) {
  mean(ifelse(w == targeted, y / treatment_probability, 0), na.rm = TRUE)
}

relative_efficiency <- function(targeted, y, w) {
  random_targeted <- integer(length(targeted))
  random_targeted[sample.int(length(targeted), sum(targeted))] <- 1L
  baseline <- policy_value(random_targeted, y, w)

  if (isTRUE(all.equal(baseline, 0))) {
    return(NA_real_)
  }

  100 * (policy_value(targeted, y, w) / baseline - 1)
}

imbalance_metric <- function(targeted, z_binary) {
  abs(mean(targeted[z_binary == 1]) - mean(targeted[z_binary == 0]))
}

delta_policy_metric <- function(targeted, targeted_twin) {
  100 * mean(targeted != targeted_twin)
}

evaluate_method <- function(test_data, method_result, target_share, cf_fd_imbalance) {
  targeted <- target_top_share(method_result$score, target_share)
  targeted_twin <- target_top_share(method_result$twin_score, target_share)
  raw_imbalance <- imbalance_metric(targeted, test_data$Z[, "Z1"])

  data.frame(
    method = method_result$method,
    efficiency = relative_efficiency(targeted, test_data$Y, test_data$W),
    raw_imbalance = raw_imbalance,
    imbalance = if (cf_fd_imbalance == 0) NA_real_ else 100 * raw_imbalance / cf_fd_imbalance,
    delta_policy = delta_policy_metric(targeted, targeted_twin),
    stringsAsFactors = FALSE
  )
}

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
