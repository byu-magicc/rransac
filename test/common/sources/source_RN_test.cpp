#include <gtest/gtest.h>
#include "rransac/common/sources/source_RN.h"
#include "rransac/common/measurement/measurement_base.h"
#include "rransac/common/transformations/transformation_null.h"



namespace rransac
{


typedef SourceRN<lie_groups::R2_r2, MeasurementTypes::RN_POS, TransformNULL> SourceR2Pos;
typedef SourceRN<lie_groups::R2_r2, MeasurementTypes::RN_POS_VEL, TransformNULL> SourceR2PosVel;

TEST(Source_RN, INIT){

// This source should only have states of type RN_rn. Make sure invalid types dont work.
SourceParameters params;
params.meas_cov_ = Eigen::Matrix2d::Identity();
params.spacial_density_of_false_meas_ = 0.1;
params.type_ = MeasurementTypes::RN_POS;
params.gate_probability_ = 0.8;
params.probability_of_detection_ = 0.9;

// This source can only handle measurement types of RN_POS and RN_POS_VEL. Make sure invalid types dont work.
params.type_ = MeasurementTypes::SEN_POS;
SourceR2Pos source_r2_pos;
SourceR2PosVel source_r2_pos_vel;
ASSERT_ANY_THROW(source_r2_pos.Init(params));

params.type_ = MeasurementTypes::SEN_POS_VEL;
ASSERT_ANY_THROW(source_r2_pos.Init(params));

params.type_ = MeasurementTypes::SEN_POSE;
ASSERT_ANY_THROW(source_r2_pos.Init(params));

params.type_ = MeasurementTypes::SEN_POSE_TWIST;
ASSERT_ANY_THROW(source_r2_pos.Init(params));

// Give it a valid type.
params.type_ = MeasurementTypes::RN_POS;
ASSERT_NO_THROW(source_r2_pos.Init(params));

params.type_ = MeasurementTypes::RN_POS_VEL;
params.meas_cov_ = Eigen::Matrix4d::Identity();
ASSERT_NO_THROW(source_r2_pos_vel.Init(params));

// Make sure parameters were set correctly.
ASSERT_EQ(params.type_, source_r2_pos_vel.GetParams().type_);

}

//////////////////////////////////////////////////
//                Test the other functions
//////////////////////////////////////////////////


TEST(Source_RN, OTHER){

// Test the functions using R2_r2
SourceParameters params;
params.meas_cov_ = SourceR2Pos::MatMeasCov::Identity();
params.spacial_density_of_false_meas_ = 0.1;
params.type_ = MeasurementTypes::RN_POS;
params.gate_probability_ = 0.8;
params.probability_of_detection_ = 0.9;
lie_groups::R2_r2 state = lie_groups::R2_r2::Random();


SourceR2Pos source1;
source1.Init(params);
Eigen::Matrix<double,2,4> H1;
H1 << 1,0,0,0,0,1,0,0;
Eigen::Matrix<double,2,2> V1;
V1.setIdentity();

typedef typename SourceR2Pos::Measurement Measurement;

Eigen::MatrixXd EmptyMat;
bool transform_state = false;

ASSERT_EQ(source1.GetLinObsMatState(state, transform_state, EmptyMat),H1);
ASSERT_EQ(source1.GetLinObsMatSensorNoise(state, transform_state, EmptyMat),V1);
Measurement m = source1.GetEstMeas(state, transform_state, EmptyMat);
ASSERT_EQ(m.pose, state.g_.data_);
m.pose.setZero(); 
m.twist.setZero();
m = SourceR2Pos::GetEstMeas(state, transform_state, EmptyMat);
ASSERT_EQ(m.pose, state.g_.data_);
m.pose.setZero(); 
m.twist.setZero();
m.type = MeasurementTypes::RN_POS_VEL;
m = SourceR2PosVel::GetEstMeas(state, transform_state, EmptyMat);
ASSERT_EQ(m.pose, state.g_.data_);
ASSERT_EQ(m.twist, state.u_.data_);

// Test OMinus
Measurement m3, m4;
m3.pose = Eigen::Matrix<double,2,1>::Random();
m3.twist = Eigen::Matrix<double,2,1>::Random();
m4.pose = Eigen::Matrix<double,2,1>::Random();
m4.twist = Eigen::Matrix<double,2,1>::Random();

Eigen::Matrix<double,4,1> error2;
error2.block(0,0,2,1) = m3.pose - m4.pose;
error2.block(2,0,2,1) = m3.twist - m4.twist;


SourceParameters params2;
params2.meas_cov_ = SourceR2PosVel::MatMeasCov::Identity();
params2.spacial_density_of_false_meas_ = 0.1;
params2.type_ = MeasurementTypes::RN_POS_VEL;
params2.gate_probability_ = 0.8;
params2.probability_of_detection_ = 0.9;


SourceR2PosVel source2;
source2.Init(params2);
m3.type = MeasurementTypes::RN_POS;
m4.type = MeasurementTypes::RN_POS;
ASSERT_EQ( source1.OMinus(m3,m4), m3.pose - m4.pose);
m3.type = MeasurementTypes::RN_POS_VEL;
m4.type = MeasurementTypes::RN_POS_VEL;
ASSERT_EQ( source2.OMinus(m3,m4), error2);

// Test Random Measurements
const int num_rand = 10000;
double std_scalar = 0.1;
typename SourceR2Pos::MatMeasCov std1 = SourceR2Pos::MatMeasCov::Identity()*std_scalar;
typename SourceR2PosVel::MatMeasCov std2 = SourceR2PosVel::MatMeasCov::Identity()*std_scalar;
std::vector<Measurement> rand_meas1(num_rand);
std::vector<Measurement> rand_meas2(num_rand);
std::vector<Eigen::Matrix<double,2,1>> error_1(num_rand);
std::vector<Eigen::Matrix<double,4,1>> error_2(num_rand);
Eigen::Matrix<double,2,1> error_mean1;
Eigen::Matrix<double,4,1> error_mean2;
error_mean1.setZero();
error_mean2.setZero();

// std::cout << "std1: " << std::endl << std1 << std::endl;

// Generate the random measurement, calculate the error between the state and measurement, and get the mean of the error
for (unsigned long int ii = 0; ii < num_rand; ++ii) {

    rand_meas1[ii].type = MeasurementTypes::RN_POS;
    rand_meas1[ii] = source1.GenerateRandomMeasurement(std1,state, transform_state, EmptyMat);
    rand_meas2[ii].type = MeasurementTypes::RN_POS_VEL;
    rand_meas2[ii] = source2.GenerateRandomMeasurement(std2,state, transform_state, EmptyMat);
    

    error_1[ii] = source1.OMinus(rand_meas1[ii], source1.GetEstMeas(state, transform_state, EmptyMat));
    error_2[ii] = source2.OMinus(rand_meas2[ii], source2.GetEstMeas(state, transform_state, EmptyMat));
    error_mean1+=error_1[ii];
    error_mean2+=error_2[ii];

}
error_mean1 /=num_rand;
error_mean2 /=num_rand;

// Calculate the covariance
typename SourceR2Pos::MatMeasCov cov1;
typename SourceR2PosVel::MatMeasCov cov2;
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
