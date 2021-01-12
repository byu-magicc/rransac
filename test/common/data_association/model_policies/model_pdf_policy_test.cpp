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



model1.Init(params);
model1.state_.g_.data_ << 0,0;
model1.state_.u_.data_ << 0,0;
model1.err_cov_ = Eigen::Matrix4d::Identity();

model2.Init(params);
model2.state_.g_.data_ << -2,-0.5;
model2.state_.u_.data_ << 1,1;
model2.err_cov_ = Eigen::Matrix4d::Identity();

model3.Init(params);
model3.state_.g_.data_ << 10,-10;
model3.state_.u_.data_ << -2,2;
model3.err_cov_ = Eigen::Matrix4d::Identity();

model4.Init(params);
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

//--------------------------------------------------------------------------------------------

bool HasMeasurement(const Model& model, const Meas& meas) {

    for(auto outer_iter = model.new_assoc_meas_.begin(); outer_iter != model.new_assoc_meas_.end(); ++ outer_iter) {
        for(auto inner_iter = outer_iter->begin(); inner_iter != outer_iter->end(); ++inner_iter) {
            if(meas.pose == inner_iter->pose)
                return true;
        }
    }
    return false;
}

bool HasMeasurementWeights(const Model& model, const Meas& meas){

    for(auto outer_iter = model.new_assoc_meas_.begin(); outer_iter != model.new_assoc_meas_.end(); ++ outer_iter) {
        for(auto inner_iter = outer_iter->begin(); inner_iter != outer_iter->end(); ++inner_iter) {
            if(meas.pose == inner_iter->pose && fabs(meas.likelihood - inner_iter->likelihood) <1e-10 && fabs(meas.weight - inner_iter->weight)< 1e-7 )
                return true;
        }
    }
    return false;
}

//---------------------------------------------------------------------------------------------

bool HasUpdateInfo(const Model& model, const ModelLikelihoodUpdateInfo& info) {
    for(auto iter = model.model_likelihood_update_info_.begin(); iter != model.model_likelihood_update_info_.end(); ++ iter) {

        if(iter->num_assoc_meas > 0) {
            if(iter->in_local_surveillance_region == info.in_local_surveillance_region && iter->num_assoc_meas == info.num_assoc_meas && iter->source_index == info.source_index && (iter->volume - info.volume) < 1e-6) {
                    return true;
                }
        } else {
            if(iter->in_local_surveillance_region == info.in_local_surveillance_region && iter->num_assoc_meas == info.num_assoc_meas && iter->source_index == info.source_index) {
                    return true;
                }
        }
  


    }
    return false;
}

//-----------------------------------------------------------------------------------------------------

void CalculateLikelihood(Meas& meas, const Model& model) {
    if (meas.source_index == 0) {
        Eigen::Matrix2d S = Eigen::Matrix2d::Identity()*2.0;
        double det = 4;
        Eigen::Matrix<double,2,1> err = meas.pose - model.state_.g_.data_;
        meas.likelihood = exp(  - (err.transpose()*S.inverse()*err)(0,0)/2.0)/sqrt(det*pow(2.0*M_PI,2));
    } else {
        Eigen::Matrix4d S = Eigen::Matrix4d::Identity()*4.0;
        double det = S.determinant();
        Eigen::Matrix<double,4,1> err;
        err.block(0,0,2,1) = meas.pose - model.state_.g_.data_;
        err.block(2,0,2,1) = meas.twist - model.state_.u_.data_;
        meas.likelihood = exp(  - (err.transpose()*S.inverse()*err)(0,0)/2.0)/sqrt(det*pow(2.0*M_PI,4));
    }
}

//-----------------------------------------------------------------------------------------------------


Meas m1, m2, m3, m4, m5, m6, m7, m8, m9;

Model model1, model2, model3, model4;
Sys sys;

DataAssociationHost< Model, ModelPDFPolicy ,DummyClusterDataTreeAssocPolicy> data_associate;

};

TEST_F(MoelPDFPolicyTestObject, PolicyDataAssociationModel) {


data_associate.AssociateNewMeasurements(sys);

ModelLikelihoodUpdateInfo info_model1_source1;
ModelLikelihoodUpdateInfo info_model1_source2;
ModelLikelihoodUpdateInfo info_model2_source1;
ModelLikelihoodUpdateInfo info_model2_source2;
ModelLikelihoodUpdateInfo info_model3_source1;
ModelLikelihoodUpdateInfo info_model3_source2;
ModelLikelihoodUpdateInfo info_model4_source1;
ModelLikelihoodUpdateInfo info_model4_source2;

info_model1_source1.in_local_surveillance_region = true;
info_model1_source1.num_assoc_meas = 2;
info_model1_source1.source_index = 0;
info_model1_source1.volume = M_PI *2.0;

info_model1_source2.in_local_surveillance_region = true;
info_model1_source2.num_assoc_meas = 1;
info_model1_source2.source_index = 1;
info_model1_source2.volume = M_PI*M_PI*8.0;

info_model2_source1.in_local_surveillance_region = true;
info_model2_source1.num_assoc_meas = 1;
info_model2_source1.source_index = 0;
info_model2_source1.volume = M_PI *2.0;

info_model2_source2.in_local_surveillance_region = true;
info_model2_source2.num_assoc_meas = 1;
info_model2_source2.source_index = 1;
info_model2_source2.volume = M_PI*M_PI*8.0;

info_model3_source1.in_local_surveillance_region = true;
info_model3_source1.num_assoc_meas = 0;
info_model3_source1.source_index = 0;
info_model3_source1.volume = M_PI *2.0;

info_model3_source2.in_local_surveillance_region = true;
info_model3_source2.num_assoc_meas = 1;
info_model3_source2.source_index = 1;
info_model3_source2.volume = M_PI*M_PI*8.0;

info_model4_source1.in_local_surveillance_region = true;
info_model4_source1.num_assoc_meas = 1;
info_model4_source1.source_index = 0;
info_model4_source1.volume = M_PI *2.0;

info_model4_source2.in_local_surveillance_region = false;
info_model4_source2.num_assoc_meas = 0;
info_model4_source2.source_index = 1;
info_model4_source2.volume = M_PI*M_PI*8.0;



Source& source1 = sys.sources_[0];
Source& source2 = sys.sources_[1];

auto model_iter = sys.models_.begin();

CalculateLikelihood(m1, *model_iter);
CalculateLikelihood(m2, *model_iter);
CalculateLikelihood(m4, *model_iter);

double PDG1 = source1.params_.probability_of_detection_/source1.params_.expected_num_false_meas_;
double PDT1 = source1.params_.probability_of_detection_*source1.params_.gate_probability_;
double PDG2 = source2.params_.probability_of_detection_/source2.params_.expected_num_false_meas_;
double PDT2 = source2.params_.probability_of_detection_*source2.params_.gate_probability_;
double den = (m1.likelihood + m2.likelihood)*PDG1 + 1 - PDT1;
m1.weight = m1.likelihood*PDG1/den;
m2.weight = m2.likelihood*PDG1/den;

den = (m4.likelihood)*PDG2 + 1 - PDT2;
m4.weight = m4.likelihood*PDG2/den;


// Model 1
ASSERT_EQ(model_iter->new_assoc_meas_.size(),2);
ASSERT_TRUE(HasMeasurementWeights(*model_iter, m1));
ASSERT_TRUE(HasMeasurementWeights(*model_iter, m2));
ASSERT_FALSE(HasMeasurement(*model_iter, m3));
ASSERT_TRUE(HasMeasurementWeights(*model_iter, m4));
ASSERT_FALSE(HasMeasurement(*model_iter, m5));
ASSERT_FALSE(HasMeasurement(*model_iter, m6));
ASSERT_FALSE(HasMeasurement(*model_iter, m7));
ASSERT_FALSE(HasMeasurement(*model_iter, m8));
ASSERT_FALSE(HasMeasurement(*model_iter, m9));
ASSERT_EQ(model_iter->model_likelihood_update_info_.size(), 2);
ASSERT_TRUE(HasUpdateInfo(*model_iter,info_model1_source1));
ASSERT_TRUE(HasUpdateInfo(*model_iter,info_model1_source2));


++model_iter;

// Model 2
CalculateLikelihood(m1, *model_iter);
CalculateLikelihood(m5, *model_iter);

den = m1.likelihood*PDG1 + 1 - PDT1;
m1.weight = m1.likelihood*PDG1/den;

den = (m5.likelihood)*PDG2 + 1 - PDT2;
m5.weight = m5.likelihood*PDG2/den;

ASSERT_TRUE(HasMeasurementWeights(*model_iter, m1));
ASSERT_FALSE(HasMeasurement(*model_iter, m2));
ASSERT_FALSE(HasMeasurement(*model_iter, m3));
ASSERT_FALSE(HasMeasurement(*model_iter, m4));
ASSERT_TRUE(HasMeasurementWeights(*model_iter, m5));
ASSERT_FALSE(HasMeasurement(*model_iter, m6));
ASSERT_FALSE(HasMeasurement(*model_iter, m7));
ASSERT_FALSE(HasMeasurement(*model_iter, m8));
ASSERT_FALSE(HasMeasurement(*model_iter, m9));
ASSERT_TRUE(HasUpdateInfo(*model_iter,info_model2_source1));
ASSERT_TRUE(HasUpdateInfo(*model_iter,info_model2_source2));




++model_iter;

// Model 3
CalculateLikelihood(m6, *model_iter);

den = (m6.likelihood)*PDG2 + 1 - PDT2;
m6.weight = m6.likelihood*PDG2/den;
ASSERT_FALSE(HasMeasurement(*model_iter, m1));
ASSERT_FALSE(HasMeasurement(*model_iter, m2));
ASSERT_FALSE(HasMeasurement(*model_iter, m3));
ASSERT_FALSE(HasMeasurement(*model_iter, m4));
ASSERT_FALSE(HasMeasurement(*model_iter, m5));
ASSERT_TRUE(HasMeasurementWeights(*model_iter, m6));
ASSERT_FALSE(HasMeasurement(*model_iter, m7));
ASSERT_FALSE(HasMeasurement(*model_iter, m8));
ASSERT_FALSE(HasMeasurement(*model_iter, m9));
ASSERT_TRUE(HasUpdateInfo(*model_iter,info_model3_source1));
ASSERT_TRUE(HasUpdateInfo(*model_iter,info_model3_source2));

++model_iter;

// Model 4
CalculateLikelihood(m3, *model_iter);
den = (m3.likelihood)*PDG1 + 1 - PDT1;
m3.weight = m3.likelihood*PDG1/den;
ASSERT_FALSE(HasMeasurement(*model_iter, m1));
ASSERT_FALSE(HasMeasurement(*model_iter, m2));
ASSERT_TRUE(HasMeasurementWeights(*model_iter, m3));
ASSERT_FALSE(HasMeasurement(*model_iter, m4));
ASSERT_FALSE(HasMeasurement(*model_iter, m5));
ASSERT_FALSE(HasMeasurement(*model_iter, m6));
ASSERT_FALSE(HasMeasurement(*model_iter, m7));
ASSERT_FALSE(HasMeasurement(*model_iter, m8));
ASSERT_FALSE(HasMeasurement(*model_iter, m9));
ASSERT_TRUE(HasUpdateInfo(*model_iter,info_model4_source1));
ASSERT_TRUE(HasUpdateInfo(*model_iter,info_model4_source2));
ASSERT_FALSE(HasUpdateInfo(*model_iter,info_model1_source1));










}