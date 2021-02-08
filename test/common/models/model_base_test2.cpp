#include <gtest/gtest.h>
#include <typeinfo>
#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>


#include "rransac/common/models/model_base.h"
#include "rransac/common/models/model_RN.h"
#include "rransac/common/models/model_SEN_pos_vel.h"
#include "rransac/common/models/model_SEN_pose_twist.h"
#include "rransac/common/sources/source_base.h"
#include "rransac/common/sources/source_RN.h"
#include "rransac/parameters.h"
#include "rransac/common/transformations/transformation_base.h"
#include "rransac/common/transformations/trans_homography.h"
#include "rransac/common/transformations/transformation_null.h"

namespace rransac
{

using namespace lie_groups;


typedef ModelRN<R3_r3,TransformNULL,SourceRN> Model2;


class ModelBaseTest2 : public ::testing::Test {
public: 

void SetUp() override {

// Setup the sources
source_params1.meas_cov_ = Eigen::Matrix3d::Identity();
source_params1.spacial_density_of_false_meas_ = 0.5;
source_params1.type_ = MeasurementTypes::RN_POS;
source_params1.gate_probability_ = 0.75;
source_params1.probability_of_detection_ = 0.95;
source_params1.source_index_ = 0;

source_params2.meas_cov_ = Eigen::Matrix<double,6,6>::Identity();
source_params2.spacial_density_of_false_meas_ = 0.2;
source_params2.type_ = MeasurementTypes::RN_POS_VEL;
source_params2.gate_probability_ = 0.8;
source_params2.probability_of_detection_ = 0.94;
source_params2.source_index_ = 1;

source_params3.meas_cov_ = Eigen::Matrix<double,6,6>::Identity();
source_params3.spacial_density_of_false_meas_ = 0.45;
source_params3.type_ = MeasurementTypes::RN_POS_VEL;
source_params3.gate_probability_ = 0.85;
source_params3.probability_of_detection_ = 0.65;
source_params3.source_index_ = 2;

source1.Init(source_params1);
source2.Init(source_params2);
source3.Init(source_params3);
sources.push_back(source1);
sources.push_back(source2);
sources.push_back(source3);

// Setup info
ModelLikelihoodUpdateInfo info1, info2, info3;

info1.in_local_surveillance_region = true;
info1.num_assoc_meas = 2;
info1.source_index = 0;
info1.volume = 2;

info2.in_local_surveillance_region = true;
info2.num_assoc_meas = 1;
info2.source_index = 1;
info2.volume = 3;

info3.in_local_surveillance_region = true;
info3.num_assoc_meas = 0;
info3.source_index = 2;
info3.volume = 1;

update_info.push_back(info1);
update_info.push_back(info2);
update_info.push_back(info3);

// Setup the parameters
params.process_noise_covariance_ = Eigen::Matrix<double,6,6>::Identity();

// Setup the model
model.Init(params);
model.model_likelihood_update_info_ = update_info;




}

Model2 model;
Model2::Source source1, source2, source3;
SourceParameters source_params1, source_params2, source_params3;
Parameters params;
std::vector<Model2::Source> sources;
std::vector<ModelLikelihoodUpdateInfo> update_info;

};


TEST_F(ModelBaseTest2, UpdateLikelihood_TEST) {
    ASSERT_EQ(this->model.model_likelihood_,0);
    this->model.UpdateModelLikelihood(this->sources);

    double tmp1 = std::log(1+this->source1.params_.gate_probability_*this->source1.params_.probability_of_detection_*(this->update_info[0].num_assoc_meas/(source1.params_.spacial_density_of_false_meas_*update_info[0].volume) -1));   
    double tmp2 = std::log(1+this->source2.params_.gate_probability_*this->source2.params_.probability_of_detection_*(this->update_info[1].num_assoc_meas/(source2.params_.spacial_density_of_false_meas_*update_info[1].volume) -1));   
    double tmp3 = std::log(1+this->source3.params_.gate_probability_*this->source3.params_.probability_of_detection_*(this->update_info[2].num_assoc_meas/(source3.params_.spacial_density_of_false_meas_*update_info[2].volume) -1));   


    ASSERT_DOUBLE_EQ(this->model.model_likelihood_,tmp1+tmp2+tmp3);

    update_info[1].num_assoc_meas = 3;
    update_info[2].num_assoc_meas = 2;
    model.model_likelihood_ = 0;
    model.model_likelihood_update_info_ = update_info;

    // Update it with associated measurements so that the likelihood increases
    for (int ii = 0; ii < 100; ++ii) {
        model.model_likelihood_update_info_ = update_info;
        this->model.UpdateModelLikelihood(this->sources);
    }
    // std::cout << this->model.model_likelihood_ << std::endl;
    
    ASSERT_GE(model.model_likelihood_, 50); // It should be at about 299

    update_info[0].num_assoc_meas = 0;
    update_info[1].num_assoc_meas = 0;
    update_info[2].num_assoc_meas = 0;
    model.model_likelihood_ = 0;
    model.model_likelihood_update_info_ = update_info;

    // Update it with no associated measurements so that the likelihood decreases
    for (int ii = 0; ii < 100; ++ii) {
        model.model_likelihood_update_info_ = update_info;
        this->model.UpdateModelLikelihood(this->sources);
    }

    ASSERT_LE(model.model_likelihood_, -50); // It should be at about -344


    // std::cout << this->model.model_likelihood_ << std::endl;

}

// ----------------------------------------------------------------------------------

TEST(ModelBaseTest2_, OMinus){

ModelSENPosVel<SE2_se2,TransformNULL,SourceSENPosVel> model1, model2;

typename SE2_se2::Mat_SC cartesian;
cartesian << 1,2,0.1,3,4,0.2;

model1.state_ = model1.state_.Random();
model2.state_ = model1.state_.OPlus(cartesian);

Eigen::MatrixXd diff = model1.OMinus(model2,model1);

ASSERT_DOUBLE_EQ(diff(0,0), cartesian(0,0));
ASSERT_DOUBLE_EQ(diff(1,0), cartesian(1,0));
ASSERT_DOUBLE_EQ(diff(2,0), cartesian(2,0));
ASSERT_DOUBLE_EQ(diff(3,0), cartesian(3,0));
ASSERT_DOUBLE_EQ(diff(4,0), cartesian(5,0));

ModelSENPoseTwist<SE2_se2,TransformNULL,SourceSENPoseTwist> model3, model4;
model3.state_ = model3.state_.Random();
model4.state_ = model3.state_.OPlus(cartesian);

Eigen::MatrixXd diff2 = model3.OMinus(model4,model3);
ASSERT_LE( ( diff2 - cartesian).norm(), 1e-10);

ModelRN<R3_r3,TransformNULL,SourceRN> model5, model6;
model5.state_ = model5.state_.Random();
model6.state_ = model5.state_.OPlus(cartesian);

Eigen::MatrixXd diff3 = model5.OMinus(model6,model5);
ASSERT_LE( ( diff3 - cartesian).norm(), 1e-10);


}



} // namespace rransac
