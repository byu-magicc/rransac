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
State target, transform_data, target_transformed;
target = State::Random();
target.g_.t_*=10;
target_transformed = target;
transform_data = State::Random();

// Construct measurements
Meas<double> m1, m2, m1_transformed;
m1.source_index = 0;
m2.source_index = 0;
m1.type = MeasurementTypes::SE3_CAM_DEPTH;
m2.type = MeasurementTypes::SE3_CAM_DEPTH;

m1 = source.GetEstMeas(target);
m1_transformed = m1;
target_transformed.g_.data_ = transform_data.g_.data_*target_transformed.g_.data_;
m2 = source.GetEstMeas(target_transformed);

// Setup transform
TransformSE3CamDepth<State> trans; 
trans.SetData(transform_data.g_.data_);

std::cout << "t" << std::endl << target.g_.t_ << std::endl;
std::cout << "t" << std::endl << target_transformed.g_.t_ << std::endl;
std::cout << "Rp" << std::endl << target.g_.R_*target.u_.p_ << std::endl;

// Test the transform target
trans.TransformTrack(target,cov);
ASSERT_EQ(cov,cov_transformed);
ASSERT_EQ(target.g_.data_, target_transformed.g_.data_);
ASSERT_EQ(target.u_.data_, target_transformed.u_.data_);



// Test the transform measurements
trans.TransformMeasurement(m1_transformed);
ASSERT_EQ(m1_transformed.pose, m2.pose);
ASSERT_EQ(m1_transformed.twist, m2.twist);








}

    
} // namespace rransac
