#include <string>
#include <gtest/gtest.h>
#include <rransac.h>

namespace rransac {


//--------------------------------------------------------

TEST(RRANSACTest, AddSource) {

SourceParameters source_params1;
SourceParameters source_params2;
SourceParameters source_params3;
SourceParameters source_params4;
SourceParameters source_params5;
SourceParameters source_params6;
RRANSAC r;

// This is a valid source. Make sure we can add it
source_params1.source_id_ = 0;
source_params1.meas_cov_fixed_ = false;
source_params1.expected_num_false_meas_ = 0.1;
source_params1.type_ = SourceTypes::R2_POS;
source_params2.probability_of_detection_ = 0.89;

ASSERT_TRUE(r.AddSource(source_params1));
// You shouldn't be able to add it again
ASSERT_ANY_THROW(r.AddSource(source_params1));

// This is an invalid source since the number of false measurements
// is negative.
source_params2 = source_params1;
source_params2.source_id_ = 1;
source_params2.expected_num_false_meas_ = -0.1;

ASSERT_ANY_THROW(r.AddSource(source_params2));


// This is an invalid source since the measurement covariance is fixed
// but it is not initialized
source_params3 = source_params1;
source_params3.source_id_ = 2;
source_params3.meas_cov_fixed_ = true;
// source_params3.meas_cov_ = Eigen::Matrix2d::Identity();


ASSERT_ANY_THROW(r.AddSource(source_params3));


// This is an invalid source since the type doesn't exist
source_params4.source_id_ = 3;
source_params4.meas_cov_fixed_ = false;
source_params4.expected_num_false_meas_ = 0.1;
ASSERT_ANY_THROW(r.AddSource(source_params4));

// This is an invalid source since the source id's must be in order
source_params5 = source_params1;
source_params5.source_id_ = 2;
ASSERT_ANY_THROW(r.AddSource(source_params5));

// This is an invalid source since the probability of detection is not between 0 and 1
source_params6 = source_params1;
source_params6.probability_of_detection_ = -0.1;
ASSERT_ANY_THROW(r.AddSource(source_params6));
source_params6.probability_of_detection_ = 1.01;
ASSERT_ANY_THROW(r.AddSource(source_params6));

// Add valid sources
for (int i = 1; i < 20; ++i) {
    SourceParameters source_params = source_params1;
    source_params.source_id_ = i;
    r.AddSource(source_params);
}

// There should be 20 valid sources
ASSERT_EQ(20, params.sources_.size());

}




//-----------------------------------------------------

TEST(RRANSACTest, SetParameters) {


Parameters P;
RRANSAC r;

ASSERT_TRUE(r.SetParameters(p));

}


} // namespace rransac