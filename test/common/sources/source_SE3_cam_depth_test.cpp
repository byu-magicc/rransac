#include <gtest/gtest.h>
#include "rransac/common/sources/source_SE3_cam_depth.h"
#include "lie_groups/state.h"
#include "rransac/common/transformations/transformation_null.h"
#include "rransac/common/transformations/trans_SE3_cam_depth.h"

namespace rransac {

using namespace lie_groups;

TEST(SOURCE_SE3_CAM_DEPTH, EstimatedMeasurementAndOMinus) {
srand((unsigned int) time(0));
using State = SE3_se3;
using Source = SourceSE3CamDepth<State,MeasurementTypes::SE3_CAM_DEPTH,TransformNULL>;
typedef Source::Measurement Measurement;


SourceParameters params;
params.type_ = MeasurementTypes::SE3_CAM_DEPTH;
double noise = 0.1;
params.meas_cov_ = Eigen::Matrix<double,7,7>::Identity()*noise;
params.source_index_ = 0;

Source source;
source.Init(params);

State state;
state = state.Random();

Measurement m1,m2;
bool transform_state = false;
Eigen::MatrixXd EmptyMat;


double d = state.g_.t_.norm();
Eigen::Matrix<double,3,1> s = state.g_.t_/d;
Eigen::Matrix<double,3,1> ds = state.g_.R_*state.u_.p_/d - state.g_.t_*state.g_.t_.transpose()*state.g_.R_*state.u_.p_/pow(d,3);
m2.pose = Eigen::Matrix<double,4,1>::Zero();
m2.twist = Eigen::Matrix<double,3,1>::Zero();
m2.pose(0,0) = d;
m2.pose.block(1,0,3,1) = s;
m2.twist = ds;


// Test the Get Est Meas
m1 = source.GetEstMeas(state, transform_state, EmptyMat);
ASSERT_EQ(m1.pose, m2.pose);
ASSERT_EQ(m1.twist, m2.twist);


// std::cout << "m1 pose: " << m1.pose << std::endl;
// std::cout << "m1 twist: " << m1.twist << std::endl;



// Test the OMinus operation and the generate random measurement function
m1 = m2;
Eigen::MatrixXd error = source.OMinus(m1,m2);
Eigen::MatrixXd Rstd = Eigen::Matrix<double,7,7>::Identity()*sqrt(noise);
Eigen::Matrix<double,7,7> cov;
cov.setZero();
ASSERT_DOUBLE_EQ(error.norm(),0);

error.setZero();

m2 = source.GenerateRandomMeasurement(Rstd*0, state, transform_state, EmptyMat);
error = source.OMinus(m1,m2);
ASSERT_DOUBLE_EQ(error.norm(),0);

// Construct the covariance from the random measurements and compare it
// with the specified covariance. 
unsigned int num_iters = 1e4;
for (unsigned int ii = 0; ii < num_iters; ++ii) {
    m2 = source.GenerateRandomMeasurement(Rstd, state, transform_state, EmptyMat);
    error = source.OMinus(m1,m2);
    cov +=error*error.transpose();
}
cov/=num_iters;
ASSERT_DOUBLE_EQ(m2.pose.block(1,0,3,1).norm(),1);
ASSERT_NE(m2.pose,m1.pose);
ASSERT_NE(m2.twist,m1.twist);
ASSERT_LT( (cov-Rstd).norm(),0.9);



// Test the Jacobians
double dt = 1e-7;
Eigen::Matrix<double,7,12> Hn;
Eigen::Matrix<double,10,12> Proj;
State dstate = state;

// setup the deviations
std::vector<Eigen::Matrix<double,12,1>> dx(12);
for (int ii =0; ii < 12; ++ii) {
    dx[ii].setZero();
    dx[ii](ii) = dt;
}

// Construct the numerical approximate of the observation function Jacobian
for (int ii=0; ii < 12; ++ii) {

    dstate = state.OPlus(dx[ii]);
    m1 = source.GetEstMeas(dstate,transform_state, EmptyMat);
    m2 = source.GetEstMeas(state,transform_state, EmptyMat);

    Hn.block(0,ii,7,1) = source.OMinus(m1,m2)/dt;
}

// Construct the projection matrix
Proj.setZero();
Proj.block(0,0,7,7).setIdentity();
Proj.block(7,9,3,3).setIdentity();

// Get the analytical Jacobian
Eigen::MatrixXd Ha = source.GetLinObsMatState(state,transform_state, EmptyMat);


ASSERT_LT( (Ha-Hn*Proj.transpose()).norm(),1e-6);

Ha = source.GetLinObsMatState(state,transform_state, EmptyMat);

ASSERT_LT( (Ha-Hn*Proj.transpose()).norm(),1e-6);
}



TEST(SOURCE_SE3_CAM_DEPTH, SpatialDistanceTest) { 
srand((unsigned int) time(0));
typedef SE3_se3 State;
typedef SourceSE3CamDepth<SE3_se3,MeasurementTypes::SE3_CAM_DEPTH,TransformSE3CamDepth> Source;
typedef typename Source::Measurement Measurement;
typedef typename Source::TransformDataType TransformDataType;
typedef typename Source::Transformation Transformation;



Parameters rransac_params;

SourceParameters params;
params.type_ = Source::measurement_type_;
double noise = 0.1;
params.meas_cov_ = Source::MatMeasCov::Identity()*noise;
params.source_index_ = 0;

Source source;
source.Init(params);

Measurement m1,m2, tmp1, tmp2;
m1.source_index = 0;
m1.type = MeasurementTypes::SE3_CAM_DEPTH;
m1.transform_meas = false;
m2.source_index= 0;
m2.type = MeasurementTypes::SE3_CAM_DEPTH;
m2.transform_meas = false;

// Test without transformation
State state1_initial = State::Random(10);
State state2_initial = State::Random(10);
tmp1 = source.GetEstMeas(state1_initial,m1.transform_meas, m1.transform_data_m_t);
tmp2 = source.GetEstMeas(state2_initial,m2.transform_meas, m2.transform_data_m_t);
m1.pose = tmp1.pose;
m1.twist = tmp1.twist;
m2.pose = tmp2.pose;
m2.twist = tmp2.twist;


double d0 = source.GetSpatialDistance(m1,m2,rransac_params);

ASSERT_DOUBLE_EQ(d0, (state1_initial.g_.t_ - state2_initial.g_.t_).norm());

// Test with transformation
m1.transform_meas = true;
m1.transform_data_m_t = Transformation::GetRandomTransform(10);
tmp1 = source.GetEstMeas(state1_initial,m1.transform_meas, m1.transform_data_m_t.inverse());
m1.pose = tmp1.pose;
m1.twist = tmp1.twist;

double d1 = source.GetSpatialDistance(m1,m2,rransac_params);
ASSERT_DOUBLE_EQ(d1,d0);

double d2 = source.GetSpatialDistance(m2,m1,rransac_params);
ASSERT_DOUBLE_EQ(d2,d0);

m2.transform_meas = true;
m2.transform_data_m_t = Transformation::GetRandomTransform(10);
tmp2 = source.GetEstMeas(state2_initial,m2.transform_meas, m2.transform_data_m_t.inverse());
m2.pose = tmp2.pose;
m2.twist = tmp2.twist;


double d3 = source.GetSpatialDistance(m1,m2,rransac_params);
ASSERT_DOUBLE_EQ(d3,d0);


}


} // namespace rransac