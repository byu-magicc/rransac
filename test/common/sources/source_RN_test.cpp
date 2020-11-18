#include <gtest/gtest.h>
#include "common/sources/source_RN.h"

namespace rransac
{
    

TEST(Source_RN, INIT){

// This source should only have states of type RN_rn. Make sure invalid types dont work.
SourceParameters params;
params.type_ = MeasurementTypes::RN_POS;
SourceRN<lie_groups::SE2_se2> source_se2;
ASSERT_ANY_THROW(source_se2.Init(params));

SourceRN<lie_groups::SE3_se3> source_se3;
ASSERT_ANY_THROW(source_se3.Init(params));

SourceRN<lie_groups::SO2_so2> source_so2;
ASSERT_ANY_THROW(source_so2.Init(params));

SourceRN<lie_groups::SO3_so3> source_so3;
ASSERT_ANY_THROW(source_so3.Init(params));

// This source can only handle measurement types of RN_POS and RN_POS_VEL. Make sure invalid types dont work.
params.type_ = MeasurementTypes::SEN_POS;
SourceRN<lie_groups::R2_r2> source_r2;
ASSERT_ANY_THROW(source_r2.Init(params));

params.type_ = MeasurementTypes::SEN_POS_VEL;
ASSERT_ANY_THROW(source_r2.Init(params));

params.type_ = MeasurementTypes::SEN_POSE;
ASSERT_ANY_THROW(source_r2.Init(params));

params.type_ = MeasurementTypes::SEN_POSE_TWIST;
ASSERT_ANY_THROW(source_r2.Init(params));

// Give it a valid type.
params.type_ = MeasurementTypes::RN_POS;
ASSERT_NO_THROW(source_r2.Init(params));

params.type_ = MeasurementTypes::RN_POS_VEL;
ASSERT_NO_THROW(source_r2.Init(params));

// Make sure parameters were set correctly.
ASSERT_EQ(params.type_, source_r2.params_.type_);

}

//////////////////////////////////////////////////
//                Test the other functions
//////////////////////////////////////////////////


TEST(Source_RN, OTHER){

// Test the functions using R2_r2
SourceParameters params;
lie_groups::R2_r2 state = lie_groups::R2_r2::Random();
params.type_ = MeasurementTypes::RN_POS;
SourceRN<lie_groups::R2_r2> source1;
source1.Init(params);
Eigen::Matrix<double,2,4> H1;
H1 << 1,0,0,0,0,1,0,0;
Eigen::Matrix<double,2,2> V1;
V1.setIdentity();

ASSERT_EQ(source1.GetLinObsMatState(state),H1);
ASSERT_EQ(source1.GetLinObsMatSensorNoise(state),V1);
Meas m = source1.GetEstMeas(state);
ASSERT_EQ(m.pose, state.g_.data_);
ASSERT_EQ(m.twist, state.u_.data_);

// Test ToEuclidean
Meas m2;
m2.pose = Eigen::Matrix<double,2,1>::Zero();
m2.pose << 1,2;
m2.pose_euclidean = source1.ToEuclidean(m2);
ASSERT_EQ(m2.pose,m2.pose_euclidean);

// Test OMinus
Meas m3, m4;
m3.pose = Eigen::Matrix<double,2,1>::Random();
m3.twist = Eigen::Matrix<double,2,1>::Random();
m4.pose = Eigen::Matrix<double,2,1>::Random();
m4.twist = Eigen::Matrix<double,2,1>::Random();

Eigen::Matrix<double,4,1> error2;
error2.block(0,0,2,1) = m3.pose - m4.pose;
error2.block(2,0,2,1) = m3.twist - m4.twist;


SourceParameters params2;
params2.type_ = MeasurementTypes::RN_POS_VEL;
SourceRN<lie_groups::R2_r2> source2;
source2.Init(params2);

ASSERT_EQ( source1.OMinus(m3,m4), m3.pose - m4.pose);
ASSERT_EQ( source2.OMinus(m3,m4), error2);

}




} // namespace rransac
