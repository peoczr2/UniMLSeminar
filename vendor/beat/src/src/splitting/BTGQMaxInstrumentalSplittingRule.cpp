/*-------------------------------------------------------------------------------
  This file is part of generalized random forest (grf).

  grf is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  grf is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with grf. If not, see <http://www.gnu.org/licenses/>.
 #-------------------------------------------------------------------------------*/

#include <algorithm>
#include <cmath>

#include "BTGQMaxInstrumentalSplittingRule.h"

namespace grf {

namespace {

bool compute_child_tau(double weight_sum,
                       double outcome_sum,
                       double treatment_sum,
                       double instrument_sum,
                       double outcome_instrument_sum,
                       double treatment_instrument_sum,
                       double& tau) {
  if (std::abs(weight_sum) <= 1e-16) {
    return false;
  }

  double avg_outcome = outcome_sum / weight_sum;
  double avg_treatment = treatment_sum / weight_sum;
  double avg_instrument = instrument_sum / weight_sum;
  double avg_outcome_instrument = outcome_instrument_sum / weight_sum;
  double avg_treatment_instrument = treatment_instrument_sum / weight_sum;

  double numerator = avg_outcome_instrument - avg_outcome * avg_instrument;
  double denominator = avg_treatment_instrument - avg_treatment * avg_instrument;

  if (std::abs(denominator) <= 1e-10) {
    return false;
  }

  tau = numerator / denominator;
  return true;
}

double compute_btgq_uplift(double weight_sum_left,
                           double outcome_sum_left,
                           double treatment_sum_left,
                           double instrument_sum_left,
                           double outcome_instrument_sum_left,
                           double treatment_instrument_sum_left,
                           double target_group_sum_left,
                           double weight_sum_right,
                           double outcome_sum_right,
                           double treatment_sum_right,
                           double instrument_sum_right,
                           double outcome_instrument_sum_right,
                           double treatment_instrument_sum_right,
                           double target_group_sum_right) {
  double tau_left = 0.0;
  double tau_right = 0.0;
  bool left_valid = compute_child_tau(weight_sum_left,
                                      outcome_sum_left,
                                      treatment_sum_left,
                                      instrument_sum_left,
                                      outcome_instrument_sum_left,
                                      treatment_instrument_sum_left,
                                      tau_left);
  bool right_valid = compute_child_tau(weight_sum_right,
                                       outcome_sum_right,
                                       treatment_sum_right,
                                       instrument_sum_right,
                                       outcome_instrument_sum_right,
                                       treatment_instrument_sum_right,
                                       tau_right);

  if (!left_valid && !right_valid) {
    return 0.0;
  }

  bool choose_left = left_valid;
  if (left_valid && right_valid) {
    choose_left = tau_left >= tau_right;
  }

  double tau_target = choose_left ? tau_left : tau_right;
  if (tau_target <= 0.0) {
    return 0.0;
  }

  double weight_sum_target = choose_left ? weight_sum_left : weight_sum_right;
  double target_group_sum = choose_left ? target_group_sum_left : target_group_sum_right;
  double concentration = target_group_sum / weight_sum_target;

  return concentration * tau_target;
}

} // namespace

BTGQMaxInstrumentalSplittingRule::BTGQMaxInstrumentalSplittingRule(size_t max_num_unique_values,
                                                                   uint min_node_size,
                                                                   double alpha,
                                                                   double imbalance_penalty):
    min_node_size(min_node_size),
    alpha(alpha),
    imbalance_penalty(imbalance_penalty) {
  this->counter = new size_t[max_num_unique_values];
  this->weight_sums = new double[max_num_unique_values];
  this->sums = new double[max_num_unique_values];
  this->num_small_z = new size_t[max_num_unique_values];
  this->sums_z = new double[max_num_unique_values];
  this->sums_z_squared = new double[max_num_unique_values];
  this->outcome_sums = new double[max_num_unique_values];
  this->treatment_sums = new double[max_num_unique_values];
  this->outcome_instrument_sums = new double[max_num_unique_values];
  this->treatment_instrument_sums = new double[max_num_unique_values];
  this->target_group_sums = new double[max_num_unique_values];
}

BTGQMaxInstrumentalSplittingRule::~BTGQMaxInstrumentalSplittingRule() {
  if (counter != nullptr) {
    delete[] counter;
  }
  if (weight_sums != nullptr) {
    delete[] weight_sums;
  }
  if (sums != nullptr) {
    delete[] sums;
  }
  if (num_small_z != nullptr) {
    delete[] num_small_z;
  }
  if (sums_z != nullptr) {
    delete[] sums_z;
  }
  if (sums_z_squared != nullptr) {
    delete[] sums_z_squared;
  }
  if (outcome_sums != nullptr) {
    delete[] outcome_sums;
  }
  if (treatment_sums != nullptr) {
    delete[] treatment_sums;
  }
  if (outcome_instrument_sums != nullptr) {
    delete[] outcome_instrument_sums;
  }
  if (treatment_instrument_sums != nullptr) {
    delete[] treatment_instrument_sums;
  }
  if (target_group_sums != nullptr) {
    delete[] target_group_sums;
  }
}

bool BTGQMaxInstrumentalSplittingRule::find_best_split(const Data& data,
                                                       size_t node,
                                                       const std::vector<size_t>& possible_split_vars,
                                                       const Eigen::ArrayXXd& responses_by_sample,
                                                       const std::vector<std::vector<size_t>>& samples,
                                                       std::vector<size_t>& split_vars,
                                                       std::vector<double>& split_values,
                                                       std::vector<bool>& send_missing_left) {
  size_t num_samples = samples[node].size();

  double weight_sum_node = 0.0;
  double sum_node = 0.0;
  double sum_node_z = 0.0;
  double sum_node_z_squared = 0.0;
  double outcome_sum_node = 0.0;
  double treatment_sum_node = 0.0;
  double outcome_instrument_sum_node = 0.0;
  double treatment_instrument_sum_node = 0.0;
  double target_group_sum_node = 0.0;
  for (auto& sample : samples[node]) {
    double sample_weight = data.get_weight(sample);
    double outcome = data.get_outcome(sample);
    double treatment = data.get_treatment(sample);
    double z = data.get_instrument(sample);

    weight_sum_node += sample_weight;
    sum_node += sample_weight * responses_by_sample(sample);
    sum_node_z += sample_weight * z;
    sum_node_z_squared += sample_weight * z * z;
    outcome_sum_node += sample_weight * outcome;
    treatment_sum_node += sample_weight * treatment;
    outcome_instrument_sum_node += sample_weight * outcome * z;
    treatment_instrument_sum_node += sample_weight * treatment * z;
    target_group_sum_node += sample_weight * data.get_btgq_target_group(sample);
  }

  double size_node = sum_node_z_squared - sum_node_z * sum_node_z / weight_sum_node;
  double min_child_size = size_node * alpha;

  double mean_z_node = sum_node_z / weight_sum_node;
  size_t num_node_small_z = 0;
  for (auto& sample : samples[node]) {
    double z = data.get_instrument(sample);
    if (z < mean_z_node) {
      num_node_small_z++;
    }
  }

  size_t best_var = 0;
  double best_value = 0.0;
  double best_decrease = 0.0;
  bool best_send_missing_left = true;

  for (auto& var : possible_split_vars) {
    find_best_split_value(data,
                          node,
                          var,
                          num_samples,
                          weight_sum_node,
                          sum_node,
                          mean_z_node,
                          num_node_small_z,
                          sum_node_z,
                          sum_node_z_squared,
                          outcome_sum_node,
                          treatment_sum_node,
                          outcome_instrument_sum_node,
                          treatment_instrument_sum_node,
                          target_group_sum_node,
                          min_child_size,
                          best_value,
                          best_var,
                          best_decrease,
                          best_send_missing_left,
                          responses_by_sample,
                          samples);
  }

  if (best_decrease <= 0.0) {
    return true;
  }

  split_vars[node] = best_var;
  split_values[node] = best_value;
  send_missing_left[node] = best_send_missing_left;
  return false;
}

void BTGQMaxInstrumentalSplittingRule::find_best_split_value(const Data& data,
                                                             size_t node,
                                                             size_t var,
                                                             size_t num_samples,
                                                             double weight_sum_node,
                                                             double sum_node,
                                                             double mean_node_z,
                                                             size_t num_node_small_z,
                                                             double sum_node_z,
                                                             double sum_node_z_squared,
                                                             double outcome_sum_node,
                                                             double treatment_sum_node,
                                                             double outcome_instrument_sum_node,
                                                             double treatment_instrument_sum_node,
                                                             double target_group_sum_node,
                                                             double min_child_size,
                                                             double& best_value,
                                                             size_t& best_var,
                                                             double& best_decrease,
                                                             bool& best_send_missing_left,
                                                             const Eigen::ArrayXXd& responses_by_sample,
                                                             const std::vector<std::vector<size_t>>& samples) {
  std::vector<double> possible_split_values;
  std::vector<size_t> sorted_samples;
  data.get_all_values(possible_split_values, sorted_samples, samples[node], var);

  if (possible_split_values.size() < 2) {
    return;
  }

  size_t num_splits = possible_split_values.size() - 1;

  std::fill(counter, counter + num_splits, 0);
  std::fill(weight_sums, weight_sums + num_splits, 0);
  std::fill(sums, sums + num_splits, 0);
  std::fill(num_small_z, num_small_z + num_splits, 0);
  std::fill(sums_z, sums_z + num_splits, 0);
  std::fill(sums_z_squared, sums_z_squared + num_splits, 0);
  std::fill(outcome_sums, outcome_sums + num_splits, 0);
  std::fill(treatment_sums, treatment_sums + num_splits, 0);
  std::fill(outcome_instrument_sums, outcome_instrument_sums + num_splits, 0);
  std::fill(treatment_instrument_sums, treatment_instrument_sums + num_splits, 0);
  std::fill(target_group_sums, target_group_sums + num_splits, 0);

  size_t n_missing = 0;
  double weight_sum_missing = 0.0;
  double sum_missing = 0.0;
  double sum_z_missing = 0.0;
  double sum_z_squared_missing = 0.0;
  size_t num_small_z_missing = 0;
  double outcome_sum_missing = 0.0;
  double treatment_sum_missing = 0.0;
  double outcome_instrument_sum_missing = 0.0;
  double treatment_instrument_sum_missing = 0.0;
  double target_group_sum_missing = 0.0;

  size_t split_index = 0;
  for (size_t i = 0; i < num_samples - 1; i++) {
    size_t sample = sorted_samples[i];
    size_t next_sample = sorted_samples[i + 1];
    double sample_value = data.get(sample, var);
    double z = data.get_instrument(sample);
    double outcome = data.get_outcome(sample);
    double treatment = data.get_treatment(sample);
    double sample_weight = data.get_weight(sample);
    double target_group = data.get_btgq_target_group(sample);

    if (std::isnan(sample_value)) {
      weight_sum_missing += sample_weight;
      sum_missing += sample_weight * responses_by_sample(sample);
      ++n_missing;

      sum_z_missing += sample_weight * z;
      sum_z_squared_missing += sample_weight * z * z;
      outcome_sum_missing += sample_weight * outcome;
      treatment_sum_missing += sample_weight * treatment;
      outcome_instrument_sum_missing += sample_weight * outcome * z;
      treatment_instrument_sum_missing += sample_weight * treatment * z;
      target_group_sum_missing += sample_weight * target_group;
      if (z < mean_node_z) {
        ++num_small_z_missing;
      }
    } else {
      weight_sums[split_index] += sample_weight;
      sums[split_index] += sample_weight * responses_by_sample(sample);
      ++counter[split_index];

      sums_z[split_index] += sample_weight * z;
      sums_z_squared[split_index] += sample_weight * z * z;
      outcome_sums[split_index] += sample_weight * outcome;
      treatment_sums[split_index] += sample_weight * treatment;
      outcome_instrument_sums[split_index] += sample_weight * outcome * z;
      treatment_instrument_sums[split_index] += sample_weight * treatment * z;
      target_group_sums[split_index] += sample_weight * target_group;
      if (z < mean_node_z) {
        ++num_small_z[split_index];
      }
    }

    double next_sample_value = data.get(next_sample, var);
    if (sample_value != next_sample_value && !std::isnan(next_sample_value)) {
      ++split_index;
    }
  }

  size_t n_left = n_missing;
  double weight_sum_left = weight_sum_missing;
  double sum_left = sum_missing;
  double sum_left_z = sum_z_missing;
  double sum_left_z_squared = sum_z_squared_missing;
  size_t num_left_small_z = num_small_z_missing;
  double outcome_sum_left = outcome_sum_missing;
  double treatment_sum_left = treatment_sum_missing;
  double outcome_instrument_sum_left = outcome_instrument_sum_missing;
  double treatment_instrument_sum_left = treatment_instrument_sum_missing;
  double target_group_sum_left = target_group_sum_missing;

  for (bool send_left : {true, false}) {
    if (!send_left) {
      if (n_missing == 0) {
        break;
      }
      n_left = 0;
      weight_sum_left = 0.0;
      sum_left = 0.0;
      sum_left_z = 0.0;
      sum_left_z_squared = 0.0;
      num_left_small_z = 0;
      outcome_sum_left = 0.0;
      treatment_sum_left = 0.0;
      outcome_instrument_sum_left = 0.0;
      treatment_instrument_sum_left = 0.0;
      target_group_sum_left = 0.0;
    }

    for (size_t i = 0; i < num_splits; ++i) {
      if (i == 0 && !send_left) {
        continue;
      }

      n_left += counter[i];
      num_left_small_z += num_small_z[i];
      weight_sum_left += weight_sums[i];
      sum_left += sums[i];
      sum_left_z += sums_z[i];
      sum_left_z_squared += sums_z_squared[i];
      outcome_sum_left += outcome_sums[i];
      treatment_sum_left += treatment_sums[i];
      outcome_instrument_sum_left += outcome_instrument_sums[i];
      treatment_instrument_sum_left += treatment_instrument_sums[i];
      target_group_sum_left += target_group_sums[i];

      size_t num_left_large_z = n_left - num_left_small_z;
      if (num_left_small_z < min_node_size || num_left_large_z < min_node_size) {
        continue;
      }

      size_t n_right = num_samples - n_left;
      size_t num_right_small_z = num_node_small_z - num_left_small_z;
      size_t num_right_large_z = n_right - num_right_small_z;
      if (num_right_small_z < min_node_size || num_right_large_z < min_node_size) {
        break;
      }

      double size_left = sum_left_z_squared - sum_left_z * sum_left_z / weight_sum_left;
      if (size_left < min_child_size || (imbalance_penalty > 0.0 && size_left == 0)) {
        continue;
      }

      double weight_sum_right = weight_sum_node - weight_sum_left;
      double sum_right = sum_node - sum_left;
      double sum_right_z_squared = sum_node_z_squared - sum_left_z_squared;
      double sum_right_z = sum_node_z - sum_left_z;
      double size_right = sum_right_z_squared - sum_right_z * sum_right_z / weight_sum_right;

      if (size_right < min_child_size || (imbalance_penalty > 0.0 && size_right == 0)) {
        continue;
      }

      double decrease_left = sum_left * sum_left / weight_sum_left;
      double decrease_right = sum_right * sum_right / weight_sum_right;
      double decrease = decrease_left + decrease_right;
      decrease -= imbalance_penalty * (1.0 / size_left + 1.0 / size_right);

      double outcome_sum_right = outcome_sum_node - outcome_sum_left;
      double treatment_sum_right = treatment_sum_node - treatment_sum_left;
      double outcome_instrument_sum_right = outcome_instrument_sum_node - outcome_instrument_sum_left;
      double treatment_instrument_sum_right = treatment_instrument_sum_node - treatment_instrument_sum_left;
      double target_group_sum_right = target_group_sum_node - target_group_sum_left;

      double uplift = compute_btgq_uplift(weight_sum_left,
                                          outcome_sum_left,
                                          treatment_sum_left,
                                          sum_left_z,
                                          outcome_instrument_sum_left,
                                          treatment_instrument_sum_left,
                                          target_group_sum_left,
                                          weight_sum_right,
                                          outcome_sum_right,
                                          treatment_sum_right,
                                          sum_right_z,
                                          outcome_instrument_sum_right,
                                          treatment_instrument_sum_right,
                                          target_group_sum_right);

      decrease += data.get_btgq_lambda() * uplift;

      if (decrease > best_decrease) {
        best_value = possible_split_values[i];
        best_var = var;
        best_decrease = decrease;
        best_send_missing_left = send_left;
      }
    }
  }
}

} // namespace grf
