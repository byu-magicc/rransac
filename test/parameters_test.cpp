#include <string>
#include <gtest/gtest.h>
#include "parameters.h"
#include <time.h>       /* time */
#include <stdlib.h>     /* srand, rand */

namespace rransac {



// This verifies that the default constructor 
// sets all of the variables properly.
TEST(ParametersTest, constructor) {

Parameters P;

ASSERT_GT(P.meas_time_window_, 0.0);
ASSERT_FALSE(P.transform_consensus_set_);

ASSERT_GT(P.cluster_time_threshold_, 0.0);
ASSERT_GT(P.cluster_velocity_threshold_, 0.0);
ASSERT_GT(P.cluster_position_threshold_, 0.0);
ASSERT_GT(P.cluster_min_size_requirement_, 0.0);

ASSERT_GT(P.RANSAC_max_iters_, 0.0);
ASSERT_GT(P.RANSAC_score_stopping_criteria_, 0.0);
ASSERT_GT(P.RANSAC_score_minimum_requirement_, 0.0);
ASSERT_GT(P.RANSAC_minimum_subset_, 0.0);

ASSERT_GT(P.track_good_model_threshold_, 0.0);
ASSERT_GT(P.track_max_missed_detection_time_, 0.0);
ASSERT_GT(P.track_similar_tracks_threshold_, 0.0);
ASSERT_GT(P.track_max_num_tracks_, 0.0);

ASSERT_FALSE(P.nonlinear_innov_cov_id_);
ASSERT_GT(P.nonlinear_LMLE_Ceres_threads_, 0.0);
ASSERT_GT(P.nonlinear_LMLE_Ceres_max_num_iters_, 0.0);



}
//---------------------------------------------------------------------------------------

// This verifies that the SetParameter function works properly
TEST(ParametersTest, SetParameter) {

// Use P1 to set good values to the parameters, pass them to P and 
// verify that they were copied correctly.
Parameters P;
Parameters P1;
double add = 5.0;
P1.meas_time_window_+=add;
P1.transform_consensus_set_ = true;
P1.cluster_time_threshold_+=add;
P1.cluster_velocity_threshold_+=add;
P1.cluster_position_threshold_+=add;
P1.cluster_min_size_requirement_+=add;
P1.RANSAC_max_iters_+=add;
P1.RANSAC_score_stopping_criteria_+=add;
P1.RANSAC_score_minimum_requirement_+=add;
P1.RANSAC_minimum_subset_+=add;
P1.track_good_model_threshold_+=add;
P1.track_max_missed_detection_time_+=add;
P1.track_similar_tracks_threshold_+=add;
P1.track_max_num_tracks_+=add;
P1.nonlinear_innov_cov_id_ = true;
P1.nonlinear_LMLE_Ceres_threads_+=add;
P1.nonlinear_LMLE_Ceres_max_num_iters_+=add;
P1.process_noise_covariance_ = Eigen::Matrix2d::Random();

P.SetParameters(P1);
ASSERT_EQ(P.meas_time_window_, P1.meas_time_window_);
ASSERT_EQ(P.transform_consensus_set_, P1.transform_consensus_set_);
ASSERT_EQ(P.cluster_time_threshold_, P1.cluster_time_threshold_);
ASSERT_EQ(P.cluster_velocity_threshold_, P1.cluster_velocity_threshold_);
ASSERT_EQ(P.cluster_position_threshold_, P1.cluster_position_threshold_);
ASSERT_EQ(P.cluster_min_size_requirement_, P1.cluster_min_size_requirement_);
ASSERT_EQ(P.RANSAC_max_iters_, P1.RANSAC_max_iters_);
ASSERT_EQ(P.RANSAC_score_stopping_criteria_, P1.RANSAC_score_stopping_criteria_);
ASSERT_EQ(P.RANSAC_score_minimum_requirement_, P1.RANSAC_score_minimum_requirement_);
ASSERT_EQ(P.RANSAC_minimum_subset_, P1.RANSAC_minimum_subset_);
ASSERT_EQ(P.track_good_model_threshold_, P1.track_good_model_threshold_);
ASSERT_EQ(P.track_max_missed_detection_time_, P1.track_max_missed_detection_time_);
ASSERT_EQ(P.track_similar_tracks_threshold_, P1.track_similar_tracks_threshold_);
ASSERT_EQ(P.track_max_num_tracks_, P1.track_max_num_tracks_);
ASSERT_EQ(P.nonlinear_innov_cov_id_, P1.nonlinear_innov_cov_id_);
ASSERT_EQ(P.nonlinear_LMLE_Ceres_threads_, P1.nonlinear_LMLE_Ceres_threads_);
ASSERT_EQ(P.nonlinear_LMLE_Ceres_max_num_iters_, P1.nonlinear_LMLE_Ceres_max_num_iters_);



srand (time(NULL));
// Use P1 to set bad values for the parameters that must be checked, pass
// them to P and verify the values.
P1.meas_time_window_ = -(rand() % 10 + 1);
ASSERT_ANY_THROW(P.SetParameters(P1));
P1.cluster_time_threshold_ = -(rand() % 10 + 1);
ASSERT_ANY_THROW(P.SetParameters(P1));
P1.cluster_velocity_threshold_ = -(rand() % 10 + 1);
ASSERT_ANY_THROW(P.SetParameters(P1));
P1.cluster_position_threshold_ = -(rand() % 10 + 1);
ASSERT_ANY_THROW(P.SetParameters(P1));
P1.cluster_min_size_requirement_ = -(rand() % 10 + 1);
ASSERT_ANY_THROW(P.SetParameters(P1));
P1.RANSAC_max_iters_ = -(rand() % 10 + 1);
ASSERT_ANY_THROW(P.SetParameters(P1));
P1.RANSAC_score_stopping_criteria_ = -(rand() % 10 + 1);
ASSERT_ANY_THROW(P.SetParameters(P1));
P1.RANSAC_score_minimum_requirement_ = -(rand() % 10 + 1);
ASSERT_ANY_THROW(P.SetParameters(P1));
P1.RANSAC_minimum_subset_ = -(rand() % 10 + 1);
ASSERT_ANY_THROW(P.SetParameters(P1));
P1.track_good_model_threshold_ = -(rand() % 10 + 1);
ASSERT_ANY_THROW(P.SetParameters(P1));
P1.track_max_missed_detection_time_ = -(rand() % 10 + 1);
ASSERT_ANY_THROW(P.SetParameters(P1));
P1.track_similar_tracks_threshold_ = -(rand() % 10 + 1);
ASSERT_ANY_THROW(P.SetParameters(P1));
P1.track_max_num_tracks_ = -(rand() % 10 + 1);
ASSERT_ANY_THROW(P.SetParameters(P1));
P1.nonlinear_LMLE_Ceres_threads_ = -(rand() % 10 + 1);
ASSERT_ANY_THROW(P.SetParameters(P1));
P1.nonlinear_LMLE_Ceres_max_num_iters_ = -(rand() % 10 + 1);
ASSERT_ANY_THROW(P.SetParameters(P1));


}

//---------------------------------------------------------------------------------------

// This verifies that the Copy Constructor works properly
TEST(ParametersTest, CopyConstructor) {

// Use P1 to set good values to the parameters, pass them to P and 
// verify that they were copied correctly.
Parameters P1;
double add = 5.0;
P1.meas_time_window_+=add;
P1.transform_consensus_set_ = true;
P1.cluster_time_threshold_+=add;
P1.cluster_velocity_threshold_+=add;
P1.cluster_position_threshold_+=add;
P1.cluster_min_size_requirement_+=add;
P1.RANSAC_max_iters_+=add;
P1.RANSAC_score_stopping_criteria_+=add;
P1.RANSAC_score_minimum_requirement_+=add;
P1.RANSAC_minimum_subset_+=add;
P1.track_good_model_threshold_+=add;
P1.track_max_missed_detection_time_+=add;
P1.track_similar_tracks_threshold_+=add;
P1.track_max_num_tracks_+=add;
P1.nonlinear_innov_cov_id_ = true;
P1.nonlinear_LMLE_Ceres_threads_+=add;
P1.nonlinear_LMLE_Ceres_max_num_iters_+=add;
P1.process_noise_covariance_ = Eigen::Matrix2d::Random();

Parameters P(P1);

ASSERT_EQ(P.meas_time_window_, P1.meas_time_window_);
ASSERT_EQ(P.transform_consensus_set_, P1.transform_consensus_set_);
ASSERT_EQ(P.cluster_time_threshold_, P1.cluster_time_threshold_);
ASSERT_EQ(P.cluster_velocity_threshold_, P1.cluster_velocity_threshold_);
ASSERT_EQ(P.cluster_position_threshold_, P1.cluster_position_threshold_);
ASSERT_EQ(P.cluster_min_size_requirement_, P1.cluster_min_size_requirement_);
ASSERT_EQ(P.RANSAC_max_iters_, P1.RANSAC_max_iters_);
ASSERT_EQ(P.RANSAC_score_stopping_criteria_, P1.RANSAC_score_stopping_criteria_);
ASSERT_EQ(P.RANSAC_score_minimum_requirement_, P1.RANSAC_score_minimum_requirement_);
ASSERT_EQ(P.RANSAC_minimum_subset_, P1.RANSAC_minimum_subset_);
ASSERT_EQ(P.track_good_model_threshold_, P1.track_good_model_threshold_);
ASSERT_EQ(P.track_max_missed_detection_time_, P1.track_max_missed_detection_time_);
ASSERT_EQ(P.track_similar_tracks_threshold_, P1.track_similar_tracks_threshold_);
ASSERT_EQ(P.track_max_num_tracks_, P1.track_max_num_tracks_);
ASSERT_EQ(P.nonlinear_innov_cov_id_, P1.nonlinear_innov_cov_id_);
ASSERT_EQ(P.nonlinear_LMLE_Ceres_threads_, P1.nonlinear_LMLE_Ceres_threads_);
ASSERT_EQ(P.nonlinear_LMLE_Ceres_max_num_iters_, P1.nonlinear_LMLE_Ceres_max_num_iters_);




}



}