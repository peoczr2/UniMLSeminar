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

#include <map>
#include <sstream>
#include <vector>
#include <iostream>
#include "commons/globals.h"
#include "Eigen/Sparse"
#include "forest/ForestPredictors.h"
#include "forest/ForestTrainers.h"
#include "RcppUtilities.h"

using namespace grf;

// [[Rcpp::export]]
Rcpp::List btgq_dumb_causal_train(Rcpp::NumericMatrix train_matrix,
                                  Eigen::SparseMatrix<double> sparse_train_matrix,
                                  size_t outcome_index,
                                  size_t treatment_index,
                                  size_t sample_weight_index,
                                  bool use_sample_weights,
                                  std::vector<double> btgq_target_group,
                                  double btgq_lambda,
                                  unsigned int mtry,
                                  unsigned int num_trees,
                                  unsigned int min_node_size,
                                  double sample_fraction,
                                  bool honesty,
                                  double honesty_fraction,
                                  bool honesty_prune_leaves,
                                  size_t ci_group_size,
                                  double reduced_form_weight,
                                  double alpha,
                                  double imbalance_penalty,
                                  bool stabilize_splits,
                                  std::vector<size_t> clusters,
                                  unsigned int samples_per_cluster,
                                  bool compute_oob_predictions,
                                  unsigned int num_threads,
                                  unsigned int seed) {
  ForestTrainer trainer = btgq_dumb_instrumental_trainer(reduced_form_weight, stabilize_splits);

  std::unique_ptr<Data> data = RcppUtilities::convert_data(train_matrix, sparse_train_matrix);

  data->set_outcome_index(outcome_index);
  data->set_treatment_index(treatment_index);
  data->set_instrument_index(treatment_index);
  if (use_sample_weights) {
    data->set_weight_index(sample_weight_index);
  }
  data->set_btgq_target_group(btgq_target_group);
  data->set_btgq_lambda(btgq_lambda);

  ForestOptions options(num_trees, ci_group_size, sample_fraction, mtry, min_node_size, honesty,
                        honesty_fraction, honesty_prune_leaves, alpha, imbalance_penalty, num_threads, seed, clusters, samples_per_cluster);
  Forest forest = trainer.train(*data, options);

  std::vector<Prediction> predictions;
  if (compute_oob_predictions) {
    ForestPredictor predictor = instrumental_predictor(num_threads);
    predictions = predictor.predict_oob(forest, *data, false);
  }

  return RcppUtilities::create_forest_object(forest, predictions);
}
