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
Eigen::Matrix<double,3,1> ds = state.g_.R_*state.u_.p_/d - state.g_.t_*state.g_.t_.transpose()*state.g_.R_*state.u_.p_/powf(d,3);
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

// Test the generate random measurement function






// Test the OMinus operation and the generate random measurement function
m1 = m2;
Eigen::MatrixXd error = source.OMinus(m1,m2);
ASSERT_DOUBLE_EQ(error.norm(),0);

error.setZero();

unsigned int num_iters = 1e3;
for (unsigned int ii = 0; ii < num_iters; ++ii) {
    m2 = source.GenerateRandomMeasurement(state,Eigen::Matrix<double,6,6>::Identity()*sqrt(noise));
    error += source.OMinus(m1,m2);
}

Eigen::Matrix<double,6,6> cov = error*error.transpose()/num_iters;
std::cout << "m2: " << std::endl << m2.pose << std::endl << m2.twist << std::endl;
std::cout << "cov: " << std::endl << cov << std::endl;





}


} // namespace rransac