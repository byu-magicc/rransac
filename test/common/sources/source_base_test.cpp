#include <gtest/gtest.h>
#include "common/sources/source_base.h"
#include <stdlib.h>     /* srand, rand */


namespace rransac {


///////////////////////////////////////////////////////////////////////
//                             Distance_Test
///////////////////////////////////////////////////////////////////////

TEST(SOURCE_BASE, TemporalDistance) {

/* initialize random seed: */
srand (time(NULL));
Parameters params_;
SourceBase<lie_groups::R2_r2> source;
Meas m1, m2;
m1.time_stamp = rand() % 100 -50;
m2.time_stamp = rand() % 100 -50;

ASSERT_EQ(source.GetTemporalDistance(m1,m2,params_),fabs(m1.time_stamp-m2.time_stamp));

}

TEST(SOURCE_BASE, SpatialDistance) {

/* initialize random seed: */
srand (time(NULL));
Parameters params_;
SourceBase<lie_groups::R2_r2> source;
Meas m_R2_Pose_1, m_R2_Pose_2, m_R2_Pose_Twist_1, m_R2_Pose_Twist_2;
m_R2_Pose_1.data = Eigen::Matrix<double,2,1>::Random();
m_R2_Pose_1.type = MeasurementTypes::R2_POSE;
m_R2_Pose_2.data = Eigen::Matrix<double,2,1>::Random();
m_R2_Pose_2.type = MeasurementTypes::R2_POSE;

m_R2_Pose_Twist_1.data = Eigen::Matrix<double,4,1>::Random();
m_R2_Pose_Twist_2.data = Eigen::Matrix<double,4,1>::Random();
m_R2_Pose_Twist_1.type = MeasurementTypes::R2_POSE_TWIST;
m_R2_Pose_Twist_2.type = MeasurementTypes::R2_POSE_TWIST;

ASSERT_DOUBLE_EQ(source.GetSpatialDistance( m_R2_Pose_1, m_R2_Pose_2,params_), (m_R2_Pose_1.data-m_R2_Pose_2.data).norm());
ASSERT_DOUBLE_EQ(source.GetSpatialDistance( m_R2_Pose_1, m_R2_Pose_Twist_1,params_), (m_R2_Pose_1.data-m_R2_Pose_Twist_1.data.block(0,0,2,1)).norm());
ASSERT_DOUBLE_EQ(source.GetSpatialDistance( m_R2_Pose_Twist_1, m_R2_Pose_Twist_2,params_), (m_R2_Pose_Twist_1.data.block(0,0,2,1)-m_R2_Pose_Twist_2.data.block(0,0,2,1)).norm());
ASSERT_DOUBLE_EQ(source.GetSpatialDistance( m_R2_Pose_Twist_1, m_R2_Pose_1,params_), (m_R2_Pose_Twist_1.data.block(0,0,2,1)-m_R2_Pose_1.data).norm());

}

///////////////////////////////////////////////////////////////////////
//                             R2_POS
///////////////////////////////////////////////////////////////////////


// TEST(SOURCE_BASE, R2_POS) {

// // Construct a state
// lie_groups::R2_r2 state;
// Eigen::Matrix<double,2,1> g;
// g << 1,2;
// state.g_.data_ = g;

// // Get parameters
// SourceParameters params;
// params.source_index_ = 100;
// params.type_ = MeasurementTypes::R2_POSE;
// params.expected_num_false_meas_ = 0.1;
// Eigen::Matrix2d R;
// R.setRandom();
// unsigned int id = 100;
// params.meas_cov_ = R;
// params.meas_cov_fixed_ = false;

// // Initialize source
// SourceBase source;
// const MeasurementTypes type = params.type_;
// source.Init<type>(params);

// // Construct Jacobians
// Eigen::Matrix<double,2,4> H;
// H << 1,0,0,0,0,1,0,0;

// Eigen::Matrix2d V;
// V.setIdentity();

// // Test to see that parameters were set properly
// ASSERT_EQ(params.expected_num_false_meas_, source.params_.expected_num_false_meas_);
// ASSERT_EQ(params.meas_cov_, source.params_.meas_cov_);
// ASSERT_EQ(params.meas_cov_fixed_, source.params_.meas_cov_fixed_);
// ASSERT_EQ(params.source_index_, source.params_.source_index_);
// ASSERT_EQ(params.type_, source.params_.type_);

// Test Jacobia.
// ASSERT_EQ(H, source.GetLinObsMatState<source.params_.type_>(state));
// ASSERT_EQ(V, source.GetLinObsMatSensorNoise<source.params_.type_>(state));

// source->GetEstMeas<rransac::SourceTypes::R2_POS>(state)

// Test Estimated measurement
// ASSERT_EQ(g, source.GetEstMeas<rransac::SourceTypes::R2_POS>(state));


// }


///////////////////////////////////////////////////////////////////////
//                             R2_POS_VEL
///////////////////////////////////////////////////////////////////////


// TEST(SOURCE_BASE, R2_POS_VEL) {

// // Construct a state
// lie_groups::R2_r2 state;
// Eigen::Matrix<double,2,1> g;
// Eigen::Matrix<double,2,1> u;
// g.setRandom();
// u.setRandom();
// state.g_.data_ = g;
// state.u_.data_ = u;

// // Get parameters
// rransac::SourceParameters params;
// params.source_id_ = 100;
// params.type_ = SourceTypes::R2_POS_VEL;
// params.expected_num_false_meas_ = 0.1;
// Eigen::Matrix4d R;
// R.setRandom();
// unsigned int id = 100;
// params.meas_cov_ = R;
// params.meas_cov_fixed_ = false;

// // Initialize source
// SourceBase source;
// source.Init<SourceTypes::R2_POS_VEL>(params);

// // Construct Jacobians
// Eigen::Matrix<double,4,4> H;
// H.setIdentity();

// Eigen::Matrix4d V;
// V.setIdentity();

// // Construct expected measurement
// Eigen::Matrix<double,4,1> meas;
// meas << g(0), g(1), u(0), u(1);


// // Test to see that parameters were set properly
// ASSERT_EQ(params.expected_num_false_meas_, source.params_.expected_num_false_meas_);
// ASSERT_EQ(params.meas_cov_, source.params_.meas_cov_);
// ASSERT_EQ(params.meas_cov_fixed_, source.params_.meas_cov_fixed_);
// ASSERT_EQ(id, source.params_.source_id_);
// ASSERT_EQ(SourceTypes::R2_POS_VEL, source.params_.type_);



// // Test Jacobians
// ASSERT_EQ(H, source.GetLinObsMatState<SourceTypes::R2_POS_VEL>(state));
// ASSERT_EQ(V, source.GetLinObsMatSensorNoise<SourceTypes::R2_POS_VEL>(state));

// // Test Estimated measurement
// ASSERT_EQ(meas, source.GetEstMeas<SourceTypes::R2_POS_VEL>(state));

// // #if CMAKE_BUILD_TYPE==Debug
// // std::cout << "debug" << std::endl;
// // #endif


// }


}