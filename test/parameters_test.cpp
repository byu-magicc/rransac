#include <string>
#include <gtest/gtest.h>
#include "parameters.h"


namespace rransac {



// This verifies that the default constructor 
// sets all of the variables properly.
TEST(ParametersTest, constructor) {

Parameters P;

ASSERT_FALSE(P.fixed_time_interval_) << "Should have a default value of false.";
ASSERT_GT(P.time_interval_,0) << "Should be greater than 0.";
ASSERT_FALSE(P.transform_consensus_set_) << "Should have a default value of false.";
ASSERT_GT(P.meas_time_window_,0) <<"Should be greater than 0.";
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
P1.max_RANSAC_iters_ = 158;
P1.RANSAC_stopping_criteria_ = 0.72;

P.SetParameters(P1);


ASSERT_EQ(P.fixed_time_interval_,P1.fixed_time_interval_);
ASSERT_EQ(P.time_interval_,P1.time_interval_);
ASSERT_EQ(P.transform_consensus_set_,P1.transform_consensus_set_);
ASSERT_EQ(P.meas_time_window_,P1.meas_time_window_);
ASSERT_EQ(P.max_RANSAC_iters_,P1.max_RANSAC_iters_);
ASSERT_EQ(P.RANSAC_stopping_criteria_,P1.RANSAC_stopping_criteria_);

// Use P1 to set bad values for the parameters that must be checked, pass
// them to P and verify the values.
P1.time_interval_ = -1;
ASSERT_ANY_THROW(P.SetParameters(P1));
P1 = P;
P1.meas_time_window_ = 0;
ASSERT_ANY_THROW(P.SetParameters(P1));
P1 = P;
P1.max_RANSAC_iters_ = 0;
ASSERT_ANY_THROW(P.SetParameters(P1));
P1 = P;
P1.RANSAC_stopping_criteria_ = 0;
ASSERT_ANY_THROW(P.SetParameters(P1));
P1 = P;

P1.time_interval_ = -5;
ASSERT_ANY_THROW(P.SetParameters(P1));
P1 = P;
P1.meas_time_window_ = -10;
ASSERT_ANY_THROW(P.SetParameters(P1));
P1 = P;
P1.RANSAC_stopping_criteria_ = 1;
ASSERT_ANY_THROW(P.SetParameters(P1));
P1 = P;

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
P1.max_RANSAC_iters_ = 158;
P1.RANSAC_stopping_criteria_ = 0.72;

Parameters P(P1);


ASSERT_EQ(P.fixed_time_interval_,P1.fixed_time_interval_);
ASSERT_EQ(P.time_interval_,P1.time_interval_);
ASSERT_EQ(P.transform_consensus_set_,P1.transform_consensus_set_);
ASSERT_EQ(P.meas_time_window_,P1.meas_time_window_);
ASSERT_EQ(P.max_RANSAC_iters_,P1.max_RANSAC_iters_);
ASSERT_EQ(P.RANSAC_stopping_criteria_,P1.RANSAC_stopping_criteria_);


}



}