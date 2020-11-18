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

// Test Random Measurements
const int num_rand = 10000;
double std_scalar = 0.1;
Eigen::Matrix2d std1 = Eigen::Matrix2d::Identity()*std_scalar;
Eigen::Matrix4d std2 = Eigen::Matrix4d::Identity()*std_scalar;
std::vector<Meas> rand_meas1(num_rand);
std::vector<Meas> rand_meas2(num_rand);
std::vector<Eigen::Matrix<double,2,1>> error_1(num_rand);
std::vector<Eigen::Matrix<double,4,1>> error_2(num_rand);
Eigen::Matrix<double,2,1> error_mean1;
Eigen::Matrix<double,4,1> error_mean2;
error_mean1.setZero();
error_mean2.setZero();

// std::cout << "std1: " << std::endl << std1 << std::endl;

// Generate the random measurement, calculate the error between the state and measurement, and get the mean of the error
for (unsigned long int ii = 0; ii < num_rand; ++ii) {

    rand_meas1[ii] = source1.GenerateRandomMeasurement(state,std1);
    rand_meas2[ii] = source2.GenerateRandomMeasurement(state,std2);

    error_1[ii] = source1.OMinus(rand_meas1[ii], source1.GetEstMeas(state));
    error_2[ii] = source2.OMinus(rand_meas2[ii], source2.GetEstMeas(state));
    error_mean1+=error_1[ii];
    error_mean2+=error_2[ii];

}
error_mean1 /=num_rand;
error_mean2 /=num_rand;

// Calculate the covariance
Eigen::Matrix2d cov1;
Eigen::Matrix4d cov2;
cov1.setZero();
cov2.setZero();
for (unsigned long int ii = 0; ii < num_rand; ++ii){

    cov1 += (error_mean1 - error_1[ii])*(error_mean1 - error_1[ii]).transpose();
    cov2 += (error_mean2 - error_2[ii])*(error_mean2 - error_2[ii]).transpose();

}
cov1 /= num_rand;
cov2 /= num_rand;


// std::cout << "error_mean1: " << std::endl << error_mean1 << std::endl;
// std::cout << "error_mean2: " << std::endl << error_mean2 << std::endl;
// std::cout << "cov1: " << std::endl << cov1 << std::endl;
// std::cout << "cov2: " << std::endl << cov2 << std::endl;


ASSERT_LE( (error_mean1 - Eigen::Matrix<double,2,1>::Zero()).norm(), 0.1);
ASSERT_LE( (cov1 - std1*std1).norm(), 0.1);
ASSERT_LE( (error_mean2 - Eigen::Matrix<double,4,1>::Zero()).norm(), 0.1);
ASSERT_LE( (cov2 - std2*std2).norm(), 0.1);

}

} // namespace rransac
