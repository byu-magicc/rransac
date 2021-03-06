#include <gtest/gtest.h>
#include <chrono> 
#include <iostream>

#include "rransac/common/transformations/trans_SE3_cam_depth.h"
#include "rransac/common/measurement/measurement_base.h"
#include "rransac/common/sources/source_SE3_cam_depth.h"
#include "lie_groups/state.h"

namespace rransac
{

using namespace lie_groups;


TEST(TransformationSE3CamDepthTest, AllFunctions) {

// Setup source
typedef SE3_se3 State;
typedef SourceSE3CamDepth<State,MeasurementTypes::SE3_CAM_DEPTH,TransformSE3CamDepth> Source;
typedef typename Source::Measurement Measurement;
typedef typename TransformSE3CamDepth<State>::TransformDataType TransformDataType;

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
target = State::Random(10);
target.g_.t_*=10;
target_transformed = target;
transform_data = State::Random();
target_c = target;

// Construct measurements
Measurement m_init, m1, m2, m3, m_final;
m1.source_index = 0;
m2.source_index = 0;
m1.transform_state = false;
m2.transform_state = false;
m1.type = MeasurementTypes::SE3_CAM_DEPTH;
m2.type = MeasurementTypes::SE3_CAM_DEPTH;

TransformDataType EmptyMat;

m1 = source.GetEstMeas(target,false,EmptyMat);
m_init = m1;
target_transformed.g_.data_ = transform_data.g_.data_*target_transformed.g_.data_;
m_final = source.GetEstMeas(target_transformed,false,EmptyMat);
m2 = m1;
m2.transform_data_m_t = transform_data.g_.data_;
m2.transform_data_t_m = transform_data.g_.data_.inverse();
m2.transform_state = true;
m2.transform_meas = true;
m3 = m2;
m3.transform_data_m_t = transform_data.g_.data_*m3.transform_data_m_t;
m3.transform_data_t_m = m3.transform_data_t_m * transform_data.g_.data_.inverse();


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
ASSERT_TRUE(m1.transform_meas);
ASSERT_TRUE(m1.transform_state);
ASSERT_EQ(m1.pose, m2.pose);
ASSERT_EQ(m1.twist, m2.twist);
ASSERT_EQ(m1.transform_data_m_t, m2.transform_data_m_t);
ASSERT_EQ(m1.transform_data_t_m, m2.transform_data_t_m);

trans.TransformMeasurement(m1);
ASSERT_EQ(m1.pose, m3.pose);
ASSERT_EQ(m1.twist, m3.twist);
ASSERT_EQ(m1.transform_data_m_t, m3.transform_data_m_t);
ASSERT_EQ(m1.transform_data_t_m, m3.transform_data_t_m);

// Verify measurement transformation when a transformation is given.
Measurement m_transformed = trans.TransformMeasurement(m1,transform_data.g_.data_);
ASSERT_LT( (m_transformed.pose - m_final.pose).norm(), 1e-9);


// Verify transformation data
TransformDataType transformation_data;
transformation_data = TransformDataType::Identity();

ASSERT_TRUE(trans.IsAcceptableTransformData(transformation_data)); 

transformation_data(3,1) = 1;
ASSERT_FALSE(trans.IsAcceptableTransformData(transformation_data));

transformation_data(3,1) = 0;
transformation_data(1,0) = 1;
ASSERT_FALSE(trans.IsAcceptableTransformData(transformation_data)); // Rotation matrix not correct.

transformation_data(1,0) = 0;
transformation_data(1,1) = 2;
ASSERT_FALSE(trans.IsAcceptableTransformData(transformation_data)); // Determininat isn't one.

// Verify random transform
transformation_data = trans.GetRandomTransform(10);
ASSERT_TRUE(trans.IsAcceptableTransformData(transformation_data)); 

}

    
} // namespace rransac
