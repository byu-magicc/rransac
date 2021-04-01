#include <string>
#include <gtest/gtest.h>

#include "lie_groups/state.h"

#include "rransac/rransac.h"
#include "rransac/common/models/model_RN.h"
#include "rransac/common/transformations/transformation_null.h"
#include "rransac/common/transformations/trans_homography.h"
#include "rransac/common/data_association/model_policies/model_pdf_policy.h"
#include "rransac/common/data_association/cluster_data_tree_policies/data_tree_cluster_association_policy.h"
#include "rransac/track_initialization/seed_policies/null_seed_policy.h"
#include "rransac/track_initialization/lmle_policies/linear_lmle_policy.h"



using namespace rransac;
using namespace lie_groups;

template <typename State>
struct CallbackClass {

bool func(const State& state) {
    if (state.g_.data_.norm() < 5)
        return true;
    else 
        return false;
}
};

typedef ModelRN<R2_r2,TransformHomography,SourceRN> Model;
typedef RRANSACTemplateParameters<R2_r2,SourceRN,TransformHomography,ModelRN,NULLSeedPolicy,LinearLMLEPolicy,ModelPDFPolicy, DataTreeClusterAssociationPolicy> RANSACParams;
typedef Eigen::Matrix<double,2,1> MatPose;

TEST(RRANSACTest, AddSource_AddMeasurement) {

RRANSAC<RANSACParams> rransac;
const System<RANSACParams::Model>* sys = rransac.GetSystemInformation();
CallbackClass<RANSACParams::State> call;

//
// Set System Parameters test
//

Parameters params;
params.cluster_time_threshold_ = 10;
ASSERT_ANY_THROW(rransac.SetSystemParameters(params));
params.process_noise_covariance_ = Eigen::Matrix4d::Identity();
ASSERT_NO_THROW(rransac.SetSystemParameters(params));
ASSERT_EQ(sys->params_.process_noise_covariance_,params.process_noise_covariance_);

// Init sources
SourceParameters source_params1, source_params2, source_params3, source_params4;

source_params1.meas_cov_ = Eigen::Matrix2d::Identity();
source_params1.type_ = MeasurementTypes::RN_POS;
source_params2.meas_cov_ = Eigen::Matrix2d::Identity();
source_params2.type_ = MeasurementTypes::RN_POS;
source_params3.meas_cov_ = Eigen::Matrix4d::Identity();
source_params3.type_ = MeasurementTypes::RN_POS_VEL;
source_params4.meas_cov_ = Eigen::Matrix4d::Identity();
source_params4.type_ = MeasurementTypes::RN_POS_VEL;

// Test invalid source index
source_params1.source_index_ = -1;
ASSERT_ANY_THROW(rransac.AddSource(source_params1));
ASSERT_ANY_THROW(rransac.AddSource(source_params1,std::bind(&CallbackClass<RANSACParams::State>::func, call, std::placeholders::_1)));
source_params1.source_index_ = 1;
ASSERT_ANY_THROW(rransac.AddSource(source_params1));
ASSERT_ANY_THROW(rransac.AddSource(source_params1,std::bind(&CallbackClass<RANSACParams::State>::func, call, std::placeholders::_1)));

// Add valid source index
source_params1.source_index_ = 0;
ASSERT_NO_THROW(rransac.AddSource(source_params1));

// Test invalid source index
ASSERT_ANY_THROW(rransac.AddSource(source_params1,std::bind(&CallbackClass<RANSACParams::State>::func, call, std::placeholders::_1)));

// Test valid source index
source_params2.source_index_ = 1;
ASSERT_NO_THROW(rransac.AddSource(source_params2,std::bind(&CallbackClass<RANSACParams::State>::func, call, std::placeholders::_1)));

// Test invalid source index
source_params3.source_index_ = 0;
ASSERT_ANY_THROW(rransac.AddSource(source_params3));

// Test valid source index
source_params3.source_index_ = 2;
ASSERT_NO_THROW(rransac.AddSource(source_params3));

// Test invalid source index
source_params4.source_index_ = 4;
ASSERT_ANY_THROW(rransac.AddSource(source_params4));

// Test valid source index
source_params4.source_index_ = 3;
ASSERT_NO_THROW(rransac.AddSource(source_params4));



// Make sure the size is right.
ASSERT_EQ(sys->sources_.size(),4);

typename RANSACParams::State state;
state.g_.data_ << 10,1;

// Test that the callback was set properly.
ASSERT_FALSE(sys->sources_[1].StateInsideSurveillanceRegion(state));

// Test change parameters
source_params1.source_index_ = -1;
ASSERT_ANY_THROW(rransac.ChangeSourceParameters(source_params1));
source_params1.source_index_ = sys->sources_.size();
ASSERT_ANY_THROW(rransac.ChangeSourceParameters(source_params1));
source_params1.source_index_ = 0;
source_params1.type_ = MeasurementTypes::RN_POS_VEL;
ASSERT_ANY_THROW(rransac.ChangeSourceParameters(source_params1));

source_params1.type_ = MeasurementTypes::RN_POS;
source_params1.gate_probability_ = 0.1234;
ASSERT_NO_THROW(rransac.ChangeSourceParameters(source_params1));
ASSERT_EQ(sys->sources_[0].params_.gate_probability_, source_params1.gate_probability_);

//
// Add Measurement test
//
std::list<Meas<double>> new_measurements, empty_measurements;
double time = 1;
Meas<double> m1,m2,m3,m4,m5;
m1.type = source_params1.type_;
m1.source_index = source_params1.source_index_;
m1.time_stamp = time;
m2.type = source_params2.type_;
m2.source_index = source_params2.source_index_;
m2.time_stamp = time;
m3.type = source_params3.type_;
m3.source_index = source_params3.source_index_;
m3.time_stamp = time;
m4.type = source_params4.type_;
m4.source_index = source_params4.source_index_;
m4.time_stamp = time+1;

m1.pose = MatPose::Random();
m2.pose = MatPose::Random();
m3.pose = MatPose::Random();
m4.pose = MatPose::Random();
m3.twist = MatPose::Random();
m4.twist = MatPose::Random();

new_measurements.push_back(m1);
new_measurements.push_back(m2);
new_measurements.push_back(m3);
new_measurements.push_back(m4);

// Add invalid measurements
#ifdef DEBUG_BUILD
ASSERT_ANY_THROW(rransac.AddMeasurements(new_measurements,time));
ASSERT_EQ(sys->new_meas_.size(), 0);
ASSERT_FALSE(sys->time_set_);
#endif

new_measurements.back().time_stamp = time;
ASSERT_NO_THROW(rransac.AddMeasurements(new_measurements,time));
ASSERT_EQ(sys->new_meas_.size(), 0);
ASSERT_EQ(sys->data_tree_.Size(), new_measurements.size());
ASSERT_TRUE(sys->time_set_);
ASSERT_EQ(sys->current_time_,m1.time_stamp);

// Add invalid measurements

m5.time_stamp = 0;
m5.source_index = 2;
m5.type = MeasurementTypes::RN_POS_VEL;
m5.pose = MatPose::Random();
m5.twist = MatPose::Random();
new_measurements.push_back(m5);
#ifdef DEBUG_BUILD
ASSERT_ANY_THROW(rransac.AddMeasurements(new_measurements,time));

for (auto iter = new_measurements.begin(); iter != new_measurements.end(); ++iter) {
    iter->time_stamp = 0;
}

ASSERT_ANY_THROW(rransac.AddMeasurements(new_measurements,time));

for (auto iter = new_measurements.begin(); iter != new_measurements.end(); ++iter) {
    iter->time_stamp = 2;
}

new_measurements.back().source_index = -1;
ASSERT_ANY_THROW(rransac.AddMeasurements(new_measurements,time)); // Source index out of scope
ASSERT_EQ(sys->data_tree_.Size(), new_measurements.size()-1);

new_measurements.back().source_index = sys->sources_.size();
ASSERT_ANY_THROW(rransac.AddMeasurements(new_measurements,time)); // Source index out of scope
ASSERT_EQ(sys->data_tree_.Size(), new_measurements.size()-1);

new_measurements.back().source_index = 0;
ASSERT_ANY_THROW(rransac.AddMeasurements(new_measurements,time)); // source index and measurement type do not match.
new_measurements.back().source_index = 2;
new_measurements.back().pose = Eigen::Matrix3d::Identity();
ASSERT_ANY_THROW(rransac.AddMeasurements(new_measurements,time)); // Pose not proper measurement type.
new_measurements.back().pose = MatPose::Identity();
new_measurements.back().twist = Eigen::Matrix3d::Identity();
ASSERT_ANY_THROW(rransac.AddMeasurements(new_measurements,time)); // Twist not proper measurement type.
new_measurements.back().twist = MatPose::Identity();
#endif

for (auto iter = new_measurements.begin(); iter != new_measurements.end(); ++iter) {
    iter->time_stamp = 2;
}

ASSERT_NO_THROW(rransac.AddMeasurements(new_measurements,time));

Eigen::Matrix3d transform_data; // Should transform all measurements to the zero vector
transform_data << 0, 0, 0, 0, 0, 0, 0, 0, 1;
rransac.AddMeasurements(empty_measurements,time,transform_data);
ASSERT_EQ(sys->new_meas_.size(), 0);
ASSERT_EQ(sys->data_tree_.Size(), 9);
ASSERT_TRUE(sys->time_set_);
ASSERT_EQ(sys->current_time_,time);

// std::list<Cluster<double>>

for (auto cluster_iter = sys->data_tree_.data_.begin(); cluster_iter != sys->data_tree_.data_.end(); ++ cluster_iter) {
    for (auto outer_iter = cluster_iter->data_.begin(); outer_iter != cluster_iter->data_.end(); ++ outer_iter) {
        for (auto inner_iter = outer_iter->begin(); inner_iter != outer_iter->end(); ++ inner_iter) {
            ASSERT_DOUBLE_EQ(inner_iter->pose.norm(),0);
        }
    }
}


for (auto iter = new_measurements.begin(); iter != new_measurements.end(); ++iter) {
    iter->time_stamp = 3;
}

ASSERT_NO_THROW(rransac.AddMeasurements(new_measurements,new_measurements.back().time_stamp));
ASSERT_EQ(sys->new_meas_.size(), 0);
ASSERT_EQ(sys->data_tree_.Size(), 14);
ASSERT_TRUE(sys->time_set_);
ASSERT_EQ(sys->current_time_,new_measurements.back().time_stamp);


for (auto iter = new_measurements.begin(); iter != new_measurements.end(); ++iter) {
    iter->time_stamp = time + sys->params_.meas_time_window_ + new_measurements.back().time_stamp;
}

ASSERT_NO_THROW(rransac.AddMeasurements(new_measurements,new_measurements.back().time_stamp));
ASSERT_EQ(sys->new_meas_.size(), 0);
ASSERT_EQ(sys->data_tree_.Size(), 5);
ASSERT_EQ(sys->current_time_,new_measurements.back().time_stamp);

}







