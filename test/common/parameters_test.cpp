#include <string>
#include <gtest/gtest.h>
#include "common/parameters.h"


namespace rransac {



// This verifies that the default constructor 
// sets all of the variables properly.
TEST(ParametersTest, constructor) {

Parameters P;

ASSERT_FALSE(P.fixed_time_interval_) << "Should have a default value of false.";
ASSERT_GT(P.time_interval_,0) << "Should be greater than 0.";
ASSERT_FALSE(P.transform_consensus_set_) << "Should have a default value of false.";
ASSERT_GT(P.meas_time_window_,0) <<"Should be greater than 0.";
ASSERT_GT(P.probability_of_detection_,0) <<"Should be between 0 and 1.";
ASSERT_LT(P.probability_of_detection_,1) <<"Should be between 0 and 1.";
ASSERT_GE(P.expected_num_false_meas_,0) << "Should be greater than or equal to 0.";
ASSERT_GT(P.max_RANSAC_iters_,0) << "Should be greather than 0.";
ASSERT_GT(P.RANSAC_stopping_criteria_,0) <<"Should be between 0 and 1.";
ASSERT_LT(P.RANSAC_stopping_criteria_,1) <<"Should be between 0 and 1.";

}
//---------------------------------------------------------------------------------------

// This verifies that the SetParameter function works properly
TEST(ParametersTest, SetParameter) {

// Use P1 to set good values to the parameters, pass them to P and 
// verify that they were copied correctly.
Parameters P;
Parameters P1;
P1.fixed_time_interval_ = true;
P1.time_interval_ = 2;
P1.transform_consensus_set_ = true;
P1.meas_time_window_ = 1;
P1.probability_of_detection_ = 0.5;
P1.expected_num_false_meas_ = 8;
P1.max_RANSAC_iters_ = 158;
P1.RANSAC_stopping_criteria_ = 0.72;

P.SetParameters(P1);


ASSERT_EQ(P.fixed_time_interval_,P1.fixed_time_interval_);
ASSERT_EQ(P.time_interval_,P1.time_interval_);
ASSERT_EQ(P.transform_consensus_set_,P1.transform_consensus_set_);
ASSERT_EQ(P.meas_time_window_,P1.meas_time_window_);
ASSERT_EQ(P.probability_of_detection_,P1.probability_of_detection_);
ASSERT_EQ(P.expected_num_false_meas_,P1.expected_num_false_meas_);
ASSERT_EQ(P.max_RANSAC_iters_,P1.max_RANSAC_iters_);
ASSERT_EQ(P.RANSAC_stopping_criteria_,P1.RANSAC_stopping_criteria_);

// Use P1 to set bad values for the parameters that must be checked, pass
// them to P and verify the values.
P1.time_interval_ = -1;
P1.meas_time_window_ = 0;
P1.probability_of_detection_ = 0;
P1.expected_num_false_meas_ = -1;
P1.max_RANSAC_iters_ = -1;
P1.RANSAC_stopping_criteria_ = 0;

P.SetParameters(P1);
ASSERT_GT(P.time_interval_,0) << "Should be greater than 0.";
ASSERT_GT(P.meas_time_window_,0) <<"Should be greater than 0.";
ASSERT_GT(P.probability_of_detection_,0) <<"Should be between 0 and 1.";
ASSERT_LT(P.probability_of_detection_,1) <<"Should be between 0 and 1.";
ASSERT_GE(P.expected_num_false_meas_,0) << "Should be greater than or equal to 0.";
ASSERT_GT(P.max_RANSAC_iters_,0) << "Should be greather than 0.";
ASSERT_GT(P.RANSAC_stopping_criteria_,0) <<"Should be between 0 and 1.";
ASSERT_LT(P.RANSAC_stopping_criteria_,1) <<"Should be between 0 and 1.";

// Use P1 to set bad values for the parameters that must be checked, pass
// them to P and verify the values.
P1.time_interval_ = -5;
P1.meas_time_window_ = -10;
P1.probability_of_detection_ = 1;
P1.expected_num_false_meas_ = -8;
P1.max_RANSAC_iters_ = -5;
P1.RANSAC_stopping_criteria_ = 1;

P.SetParameters(P1);
ASSERT_GT(P.time_interval_,0) << "Should be greater than 0.";
ASSERT_GT(P.meas_time_window_,0) <<"Should be greater than 0.";
ASSERT_GT(P.probability_of_detection_,0) <<"Should be between 0 and 1.";
ASSERT_LT(P.probability_of_detection_,1) <<"Should be between 0 and 1.";
ASSERT_GE(P.expected_num_false_meas_,0) << "Should be greater than or equal to 0.";
ASSERT_GT(P.max_RANSAC_iters_,0) << "Should be greather than 0.";
ASSERT_GT(P.RANSAC_stopping_criteria_,0) <<"Should be between 0 and 1.";
ASSERT_LT(P.RANSAC_stopping_criteria_,1) <<"Should be between 0 and 1.";

}

//---------------------------------------------------------------------------------------

// This verifies that the Copy Constructor works properly
TEST(ParametersTest, CopyConstructor) {

// Use P1 to set good values to the parameters, pass them to P and 
// verify that they were copied correctly.
Parameters P1;
P1.fixed_time_interval_ = true;
P1.time_interval_ = 2;
P1.transform_consensus_set_ = true;
P1.meas_time_window_ = 1;
P1.probability_of_detection_ = 0.5;
P1.expected_num_false_meas_ = 8;
P1.max_RANSAC_iters_ = 158;
P1.RANSAC_stopping_criteria_ = 0.72;

Parameters P(P1);


ASSERT_EQ(P.fixed_time_interval_,P1.fixed_time_interval_);
ASSERT_EQ(P.time_interval_,P1.time_interval_);
ASSERT_EQ(P.transform_consensus_set_,P1.transform_consensus_set_);
ASSERT_EQ(P.meas_time_window_,P1.meas_time_window_);
ASSERT_EQ(P.probability_of_detection_,P1.probability_of_detection_);
ASSERT_EQ(P.expected_num_false_meas_,P1.expected_num_false_meas_);
ASSERT_EQ(P.max_RANSAC_iters_,P1.max_RANSAC_iters_);
ASSERT_EQ(P.RANSAC_stopping_criteria_,P1.RANSAC_stopping_criteria_);


}

///--------------------------------------------------------

TEST(ParametersTest, AddSource) {

SourceParameters source_params1;
SourceParameters source_params2;
SourceParameters source_params3;
SourceParameters source_params4;
SourceParameters source_params5;
Parameters params;

// This is a valid source. Make sure we can add it
source_params1.source_id_ = 0;
source_params1.meas_cov_fixed_ = false;
source_params1.expected_num_false_meas_ = 0.1;
source_params1.type_ = SourceTypes::R2_POS;

ASSERT_TRUE(params.AddSource(source_params1));
// You shouldn't be able to add it again
ASSERT_ANY_THROW(params.AddSource(source_params1));

// This is an invalid source since the number of false measurements
// is negative.
source_params2 = source_params1;
source_params2.source_id_ = 1;
source_params2.expected_num_false_meas_ = -0.1;

ASSERT_ANY_THROW(params.AddSource(source_params2));


// This is an invalid source since the measurement covariance is fixed
// but it is not initialized
source_params3 = source_params1;
source_params3.source_id_ = 2;
source_params3.meas_cov_fixed_ = true;
// source_params3.meas_cov_ = Eigen::Matrix2d::Identity();


ASSERT_ANY_THROW(params.AddSource(source_params3));


// This is an invalid source since the type doesn't exist
source_params4.source_id_ = 3;
source_params4.meas_cov_fixed_ = false;
source_params4.expected_num_false_meas_ = 0.1;
ASSERT_ANY_THROW(params.AddSource(source_params4));

// This is an invalid source since the source id's must be in order
source_params5 = source_params1;
source_params5.source_id_ = 2;
ASSERT_ANY_THROW(params.AddSource(source_params5));

// Add valid sources
for (int i = 1; i < 20; ++i) {
    SourceParameters source_params = source_params1;
    source_params.source_id_ = i;
    params.AddSource(source_params);
}

// There should be 20 valid sources
ASSERT_EQ(20, params.sources_.size());

}

}