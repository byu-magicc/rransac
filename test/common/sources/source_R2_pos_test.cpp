#include <string>
#include <gtest/gtest.h>
#include "lie_groups/liegroups/include/state.h"
#include "common/sources/source_R2_pos.h"



TEST(Source_R2_POS_TEST, ALL) {

lie_groups::State<lie_groups::Rn<2>,lie_groups::rn<2>> state;

// Construct a state
Eigen::Matrix2d g;
g << 1,2;
state.g_.data_ = g;

// Get parameters
rransac::SourceParameters params;
params.expected_num_false_meas_ = 0.1;
Eigen::Matrix2d R;
R.setRandom();
unsigned int id = 100;
params.meas_cov_ = R;
params.meas_cov_fixed_ = false;

// Initialize source
rransac::SourceR2Pos source;
source.Init(params, id);

// Construct Jacobians
Eigen::Matrix<double,2,4> H;
H << 1,0,0,0,0,1,0,0;

Eigen::Matrix2d V;
V.setIdentity();

// Test to see that parameters were set properly
ASSERT_EQ(params.expected_num_false_meas_, source.params_.expected_num_false_meas_);
ASSERT_EQ(params.meas_cov_, source.params_.meas_cov_);
ASSERT_EQ(params.meas_cov_fixed_, source.params_.meas_cov_fixed_);

// Test to see that the source ID was set properly
ASSERT_EQ(id, source.source_id_);

// Test Jacobians
ASSERT_EQ(H, source.GetLinObsMatState());
ASSERT_EQ(V, source.GetLinObsMatSensorNoise());

// Test Estimated measurement
ASSERT_EQ(g, source.GetEstMeas(state));

}