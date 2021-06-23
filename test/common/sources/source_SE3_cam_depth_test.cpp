#include <gtest/gtest.h>
#include "rransac/common/sources/source_SE3_cam_depth.h"
#include "lie_groups/state.h"

namespace rransac {

using namespace lie_groups;

TEST(SOURCE_SE3_CAM_DEPTH, EstimatedMeasurementAndOMinus) {

using State = SE3_se3;
using Source = SourceSE3CamDepth<State>;


SourceParameters params;
params.type_ = MeasurementTypes::SE3_CAM_DEPTH;
double noise = 0.1;
params.meas_cov_ = Eigen::Matrix<double,7,7>::Identity()*noise;
params.source_index_ = 0;

Source source;
source.Init(params);

State state;
state = state.Random();

Meas<double> m1,m2;


double d = state.g_.t_.norm();
Eigen::Matrix<double,3,1> s = state.g_.t_/d;
Eigen::Matrix<double,3,1> ds = state.g_.R_*state.u_.p_/d - state.g_.t_*state.g_.t_.transpose()*state.g_.R_*state.u_.p_/pow(d,3);
m2.pose = Eigen::Matrix<double,4,1>::Zero();
m2.twist = Eigen::Matrix<double,3,1>::Zero();
m2.pose(0,0) = d;
m2.pose.block(1,0,3,1) = s;
m2.twist = ds;


// Test the Get Est Meas
m1 = source.GetEstMeas(state);
ASSERT_EQ(m1.pose, m2.pose);
ASSERT_EQ(m1.twist, m2.twist);

m1 = source.GetEstMeas(state,MeasurementTypes::SE3_CAM_DEPTH);
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

m2 = source.GenerateRandomMeasurement(state,Rstd*0);
error = source.OMinus(m1,m2);
ASSERT_DOUBLE_EQ(error.norm(),0);

// Construct the covariance from the random measurements and compare it
// with the specified covariance. 
unsigned int num_iters = 1e4;
for (unsigned int ii = 0; ii < num_iters; ++ii) {
    m2 = source.GenerateRandomMeasurement(state,Rstd);
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
    m1 = source.GetEstMeas(dstate);
    m2 = source.GetEstMeas(state);

    Hn.block(0,ii,7,1) = source.OMinus(m1,m2)/dt;
}

// Construct the projection matrix
Proj.setZero();
Proj.block(0,0,7,7).setIdentity();
Proj.block(7,9,3,3).setIdentity();

// Get the analytical Jacobian
Eigen::MatrixXd Ha = source.GetLinObsMatState(state);


ASSERT_LT( (Ha-Hn*Proj.transpose()).norm(),1e-6);

Ha = source.GetLinObsMatState(state,MeasurementTypes::SE3_CAM_DEPTH);

ASSERT_LT( (Ha-Hn*Proj.transpose()).norm(),1e-6);
}



TEST(SOURCE_SE3_CAM_DEPTH, SpatialDistanceTest) { 


using State = SE3_se3;
using Source = SourceSE3CamDepth<State>;

Parameters rransac_params;

SourceParameters params;
params.type_ = MeasurementTypes::SE3_CAM_DEPTH;
double noise = 0.1;
params.meas_cov_ = Eigen::Matrix<double,7,7>::Identity()*noise;
params.source_index_ = 0;

Source source;
source.Init(params);

State state, transform1, transform2, state_t1, state_t2;
state = state.Random();
transform1 = State::Random();
transform2 = State::Random();
state_t1 = state;
state_t1.g_.data_ = transform1.g_.data_*state_t1.g_.data_;
state_t2 = state;
state_t2.g_.data_ = transform2.g_.data_*state_t2.g_.data_;

Meas<double> m1,m2;
m1.source_index = 0;
m1.type = MeasurementTypes::SE3_CAM_DEPTH;
m2.source_index= 0;
m2.type = MeasurementTypes::SE3_CAM_DEPTH;
m1 = source.GetEstMeas(state);
m1.state_transform_data = false;
m2 = source.GetEstMeas(state);
m2.state_transform_data = false;

double d0, d1, d2, d3;

d0 = source.GetSpatialDistance(m1,m2,rransac_params);

ASSERT_EQ(d0, 0);

m1.state_transform_data = true;
m1.trans_data = transform1.g_.data_;
m2 = source.GetEstMeas(state_t1);
d1 = source.GetSpatialDistance(m1,m2,rransac_params);

ASSERT_LT(d1,1e-12);

d2 = source.GetSpatialDistance(m2,m1,rransac_params);

ASSERT_LT(d2,1e-12);

m2 = m1;

d3 = source.GetSpatialDistance(m2,m1,rransac_params);

ASSERT_LT(d3,1e-12);


}


} // namespace rransac