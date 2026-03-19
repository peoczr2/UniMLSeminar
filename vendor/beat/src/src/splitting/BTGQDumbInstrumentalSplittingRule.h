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

#ifndef GRF_BTGQDUMBINSTRUMENTALSPLITTINGRULE_H
#define GRF_BTGQDUMBINSTRUMENTALSPLITTINGRULE_H

#include "SplittingRule.h"

namespace grf {

class BTGQDumbInstrumentalSplittingRule final : public SplittingRule {
public:
  BTGQDumbInstrumentalSplittingRule(size_t max_num_unique_values,
                                    uint min_node_size,
                                    double alpha,
                                    double imbalance_penalty);
  ~BTGQDumbInstrumentalSplittingRule();

  bool find_best_split(const Data& data,
                       size_t node,
                       const std::vector<size_t>& possible_split_vars,
                       const Eigen::ArrayXXd& responses_by_sample,
                       const std::vector<std::vector<size_t>>& samples,
                       std::vector<size_t>& split_vars,
                       std::vector<double>& split_values,
                       std::vector<bool>& send_missing_left) override;

private:
  void find_best_split_value(const Data& data,
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
                             const std::vector<std::vector<size_t>>& samples);

  size_t min_node_size;
  double alpha;
  double imbalance_penalty;

  size_t* counter;
  double* weight_sums;
  double* sums;
  size_t* num_small_z;
  double* sums_z;
  double* sums_z_squared;
  double* outcome_sums;
  double* treatment_sums;
  double* outcome_instrument_sums;
  double* treatment_instrument_sums;
  double* target_group_sums;

  DISALLOW_COPY_AND_ASSIGN(BTGQDumbInstrumentalSplittingRule);
};

} // namespace grf

#endif // GRF_BTGQDUMBINSTRUMENTALSPLITTINGRULE_H
