#include <gtest/gtest.h>
#include <chrono> 
#include <iostream>

#include "rransac/common/transformations/trans_SE3_cam_depth.h"
#include "rransac/common/measurement/measurement_base.h"
#include "rransac/common/sources/source_SE3_cam_depth.h"

namespace rransac
{

using namespace lie_groups;


TEST(TransformationSE3CamDepthTest, AllFunctions) {

// Setup source
typedef typename lie_groups::SE3_se3 State;
typedef SourceSE3CamDepth<State> Source;

SourceParameters params;
params.type_ = MeasurementTypes::SE3_CAM_DEPTH;
double noise = 0.1;
params.meas_cov_ = Eigen::Matrix<double,7,7>::Identity()*noise;
params.source_index_ = 0;
Source source;
source.Init(params);

// Setup target
Eigen::Matrix<double,10,10> cov, cov_transformed;
cov.setRandom();
cov_transformed = cov;
State target, transform_data, target_transformed, target_c;
target = State::Random();
target.g_.t_*=10;
target_transformed = target;
transform_data = State::Random();
target_c = target;

// Construct measurements
Meas<double> m1, m2, m3;
m1.source_index = 0;
m2.source_index = 0;
m1.state_transform_data = false;
m2.state_transform_data = false;
m1.type = MeasurementTypes::SE3_CAM_DEPTH;
m2.type = MeasurementTypes::SE3_CAM_DEPTH;

m1 = source.GetEstMeas(target);
target_transformed.g_.data_ = transform_data.g_.data_*target_transformed.g_.data_;
m2 = m1;
m2.trans_data = transform_data.g_.data_;
m2.state_transform_data = true;
m3 = m2;
m3.trans_data = transform_data.g_.data_*m3.trans_data;


// Setup transform
TransformSE3CamDepth<State> trans; 
trans.SetData(transform_data.g_.data_);


// Test the transform target
trans.TransformTrack(target,cov);
ASSERT_EQ(cov,cov_transformed);
ASSERT_EQ(target.g_.data_, target_transformed.g_.data_);
ASSERT_EQ(target.u_.data_, target_transformed.u_.data_);

trans.SetData(Eigen::Matrix4d::Identity());
trans.TransformTrack(target_c,cov,transform_data.g_.data_);
ASSERT_EQ(cov,cov_transformed);
ASSERT_EQ(target_c.g_.data_, target_transformed.g_.data_);
ASSERT_EQ(target_c.u_.data_, target_transformed.u_.data_);
ASSERT_EQ(trans.GetData(), Eigen::Matrix4d::Identity());
trans.SetData(transform_data.g_.data_);

// Test the transform measurements
trans.TransformMeasurement(m1);
ASSERT_EQ(m1.pose, m2.pose);
ASSERT_EQ(m1.twist, m2.twist);
ASSERT_EQ(m1.trans_data, m2.trans_data);

trans.TransformMeasurement(m1);
ASSERT_EQ(m1.pose, m3.pose);
ASSERT_EQ(m1.twist, m3.twist);
ASSERT_EQ(m1.trans_data, m3.trans_data);







}

    
} // namespace rransac
