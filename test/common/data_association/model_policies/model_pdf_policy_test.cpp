#include <gtest/gtest.h>
#include "system.h"
#include "common/models/model_RN.h"
#include "common/sources/source_RN.h"
#include "common/transformations/transformation_null.h"
#include "common/data_association/data_association_host.h"
#include "common/data_association/model_policies/model_pdf_policy.h"


using namespace rransac;
using namespace lie_groups;





//----------------------------------------------------------------------------------


bool InSurveillanceRegionSource2(const R2_r2& state ) {

Eigen::Matrix<double,2,1> center;
center << 5,5;

if( (state.g_.data_ - center ).norm() < 50)
    return true;
else
{
        return false;
}

}

//----------------------------------------------------------------------------------

template<typename tModel>
class DummyClusterDataTreeAssocPolicy {
public:

    static void  PolicyDataAssociationClusterDataTree(System<tModel>& sys) {
        // Store the size as the current time so we can verify it.
        sys.current_time_ = sys.new_meas_.size();
        sys.new_meas_.clear();
    }

};

//----------------------------------------------------------------------------------

class MoelPDFPolicyTestObject : public ::testing::Test {

protected:

typedef R2_r2 State;
typedef ModelRN<State,TransformNULL> Model;
typedef Model::Source Source;
typedef System<Model> Sys;

void SetUp() override {

// None of these paramters matter for these tests. They just need to be set.
Parameters params;
params.meas_time_window_ = 10;
params.fixed_time_interval_ = true;
params.time_interval_ = 0.1;
params.transform_consensus_set_ = false;
params.cluster_position_threshold_ = 1;
params.cluster_time_threshold_ = 1;
params.cluster_position_threshold_ = 1;
params.RANSAC_minimum_subset_ = 2;
params.RANSAC_max_iters_ = 100;
params.RANSAC_stopping_criteria_ = 20;
params.process_noise_covariance_ = Eigen::Matrix4d::Identity();

// Setup the sources
SourceParameters source_params1,source_params2, source_params3;
Source source1, source2, source3;

source_params1.meas_cov_fixed_ = true;
source_params1.meas_cov_ = Eigen::Matrix2d::Identity();
source_params1.expected_num_false_meas_ = 0.1;
source_params1.type_ = MeasurementTypes::RN_POS;
source_params1.probability_of_detection_ = 0.9;
source_params1.gate_probability_ = 0.393469340287367;  // select it so that the gate threshold is 1
source_params1.source_index_ = 0;

source_params2.meas_cov_fixed_ = false;
source_params2.expected_num_false_meas_ = 0.2;
source_params2.type_ = MeasurementTypes::RN_POS_VEL;
source_params2.probability_of_detection_ = 0.9;
source_params2.gate_probability_ = 0.090204010431050;  // select it so that the gate threshold is 1
source_params2.source_index_ = 1;

source_params3 = source_params1;
source_params3.source_index_ = 2;

// All models will be in the surveillance region of source1, model 4 will not be in the surveillance region of source 2
source1.Init(source_params1);
source2.Init(source_params2, &InSurveillanceRegionSource2);
source3.Init(source_params3);  // This source will not produce any measurements

sys.sources_.push_back(source1);
sys.sources_.push_back(source2);
sys.sources_.push_back(source3);



model1.Init(sys.sources_, params);
model1.state_.g_.data_ << 0,0;
model1.state_.u_.data_ << 0,0;
model1.err_cov_ = Eigen::Matrix4d::Identity();

model2.Init(sys.sources_, params);
model2.state_.g_.data_ << -2,-0.5;
model2.state_.u_.data_ << 1,1;
model2.err_cov_ = Eigen::Matrix4d::Identity();

model3.Init(sys.sources_, params);
model3.state_.g_.data_ << 10,-10;
model3.state_.u_.data_ << -2,2;
model3.err_cov_ = Eigen::Matrix4d::Identity();

model4.Init(sys.sources_, params);
model4.state_.g_.data_ << 100, 100;
model4.state_.u_.data_ << -10, -40;
model4.err_cov_ = Eigen::Matrix4d::Identity();

sys.models_.push_back(model1);
sys.models_.push_back(model2);
sys.models_.push_back(model3);
sys.models_.push_back(model4);



// Will be in the validation region of model1 and model2
m1.time_stamp = 0;
m1.type = MeasurementTypes::RN_POS;
m1.source_index = 0;
m1.pose = Eigen::Matrix<double,2,1>::Zero();
m1.pose << -1.3, 0;

// Will be in the validation region of model1
m2.time_stamp = 0;
m2.type = MeasurementTypes::RN_POS;
m2.source_index = 0;
m2.pose = Eigen::Matrix<double,2,1>::Zero();
m2.pose << -0.1, 0.9;

// Will be in the validation region of model 4
m3.time_stamp = 0;
m3.type = MeasurementTypes::RN_POS;
m3.source_index = 0;
m3.pose = Eigen::Matrix<double,2,1>::Zero();
m3.pose << 99.5, 100.5;

// Will be in the validation region of model 1
m4.time_stamp = 0;
m4.type = MeasurementTypes::RN_POS_VEL;
m4.source_index = 1;
m4.pose = Eigen::Matrix<double,2,1>::Zero();
m4.pose << 0.5, 0.5;
m4.twist = Eigen::Matrix<double,2,1>::Zero();
m4.twist << -0.5, 0.5;
m4.meas_cov = Eigen::Matrix4d::Identity()*3;

// Will be in the validation region of model 2
m5.time_stamp = 0;
m5.type = MeasurementTypes::RN_POS_VEL;
m5.source_index = 1;
m5.pose = Eigen::Matrix<double,2,1>::Zero();
m5.pose << -2.3, 0;
m5.twist = Eigen::Matrix<double,2,1>::Zero();
m5.twist << 0.7, 1.3;
m5.meas_cov = Eigen::Matrix4d::Identity()*3;


// Will be in the validation region of model 3
m6.time_stamp = 0;
m6.type = MeasurementTypes::RN_POS_VEL;
m6.source_index = 1;
m6.pose = Eigen::Matrix<double,2,1>::Zero();
m6.pose << 10.1, -9.6;
m6.twist = Eigen::Matrix<double,2,1>::Zero();
m6.twist << -2.3, 1.8;
m6.meas_cov = Eigen::Matrix4d::Identity()*3;


// Noisy Measurement
m7.time_stamp = 0;
m7.type = MeasurementTypes::RN_POS_VEL;
m7.source_index = 1;
m7.pose = Eigen::Matrix<double,2,1>::Zero();
m7.pose << -100, -100;
m7.twist = Eigen::Matrix<double,2,1>::Zero();
m7.twist << -2.3, 1.8;
m7.meas_cov = Eigen::Matrix4d::Identity()*3;


// Noisy Measurement
m8.time_stamp = 0;
m8.type = MeasurementTypes::RN_POS;
m8.source_index = 0;
m8.pose = Eigen::Matrix<double,2,1>::Zero();
m8.pose << 100, -100;
m8.twist = Eigen::Matrix<double,2,1>::Zero();
m8.twist << -2.3, 1.8;

// Noisy Measurement
m9.time_stamp = 0;
m9.type = MeasurementTypes::RN_POS;
m9.source_index = 0;
m9.pose = Eigen::Matrix<double,2,1>::Zero();
m9.pose << 100, -100;
m9.twist = Eigen::Matrix<double,2,1>::Zero();
m9.twist << -50, 1.8;

sys.new_meas_.push_back(m9);
sys.new_meas_.push_back(m1);
sys.new_meas_.push_back(m2);
sys.new_meas_.push_back(m3);
sys.new_meas_.push_back(m7);
sys.new_meas_.push_back(m4);
sys.new_meas_.push_back(m5);
sys.new_meas_.push_back(m6);
sys.new_meas_.push_back(m8);

}

bool HasMeasurement(const Model& model, const Meas& meas) {

    for(auto outer_iter = model.new_assoc_meas_.begin(); outer_iter != model.new_assoc_meas_.end(); ++ outer_iter) {
        for(auto inner_iter = outer_iter->begin(); inner_iter != outer_iter->end(); ++inner_iter) {
            if(meas.pose == inner_iter->pose)
                return true;
        }
    }
    return false;
}

Meas m1, m2, m3, m4, m5, m6, m7, m8, m9;

Model model1, model2, model3, model4;
Sys sys;

DataAssociationHost< Model, ModelPDFPolicy ,DummyClusterDataTreeAssocPolicy> data_associate;

};

TEST_F(MoelPDFPolicyTestObject, PolicyDataAssociationModel) {


data_associate.AssociateNewMeasurements(sys);

auto model_iter = sys.models_.begin();

// Model 1
ASSERT_EQ(model_iter->new_assoc_meas_.size(),2);
ASSERT_TRUE(HasMeasurement(*model_iter, m1));
ASSERT_TRUE(HasMeasurement(*model_iter, m2));
ASSERT_FALSE(HasMeasurement(*model_iter, m3));
ASSERT_TRUE(HasMeasurement(*model_iter, m4));
ASSERT_FALSE(HasMeasurement(*model_iter, m5));
ASSERT_FALSE(HasMeasurement(*model_iter, m6));
ASSERT_FALSE(HasMeasurement(*model_iter, m7));
ASSERT_FALSE(HasMeasurement(*model_iter, m8));
ASSERT_FALSE(HasMeasurement(*model_iter, m9));

++model_iter;

// Model 2
ASSERT_TRUE(HasMeasurement(*model_iter, m1));
ASSERT_FALSE(HasMeasurement(*model_iter, m2));
ASSERT_FALSE(HasMeasurement(*model_iter, m3));
ASSERT_FALSE(HasMeasurement(*model_iter, m4));
ASSERT_TRUE(HasMeasurement(*model_iter, m5));
ASSERT_FALSE(HasMeasurement(*model_iter, m6));
ASSERT_FALSE(HasMeasurement(*model_iter, m7));
ASSERT_FALSE(HasMeasurement(*model_iter, m8));
ASSERT_FALSE(HasMeasurement(*model_iter, m9));

++model_iter;

// Model 3
ASSERT_FALSE(HasMeasurement(*model_iter, m1));
ASSERT_FALSE(HasMeasurement(*model_iter, m2));
ASSERT_FALSE(HasMeasurement(*model_iter, m3));
ASSERT_FALSE(HasMeasurement(*model_iter, m4));
ASSERT_FALSE(HasMeasurement(*model_iter, m5));
ASSERT_TRUE(HasMeasurement(*model_iter, m6));
ASSERT_FALSE(HasMeasurement(*model_iter, m7));
ASSERT_FALSE(HasMeasurement(*model_iter, m8));
ASSERT_FALSE(HasMeasurement(*model_iter, m9));

++model_iter;

// Model 4
ASSERT_FALSE(HasMeasurement(*model_iter, m1));
ASSERT_FALSE(HasMeasurement(*model_iter, m2));
ASSERT_TRUE(HasMeasurement(*model_iter, m3));
ASSERT_FALSE(HasMeasurement(*model_iter, m4));
ASSERT_FALSE(HasMeasurement(*model_iter, m5));
ASSERT_FALSE(HasMeasurement(*model_iter, m6));
ASSERT_FALSE(HasMeasurement(*model_iter, m7));
ASSERT_FALSE(HasMeasurement(*model_iter, m8));
ASSERT_FALSE(HasMeasurement(*model_iter, m9));











}