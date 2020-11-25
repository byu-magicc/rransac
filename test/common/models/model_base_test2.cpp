#include <gtest/gtest.h>
#include <typeinfo>
#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>


#include "common/models/model_base.h"
#include "common/models/model_RN.h"
#include "common/models/model_SEN_pos_vel.h"
#include "common/models/model_SEN_pose_twist.h"
#include "common/sources/source_base.h"
#include "common/sources/source_RN.h"
#include "parameters.h"
#include "common/transformations/transformation_base.h"
#include "common/transformations/trans_homography.h"
#include "common/transformations/transformation_null.h"

namespace rransac
{

using namespace lie_groups;

// typedef ModelBase<SourceR2, TransformNULL<R2_r2>, 4, ModelRN<R2_r2, TransformNULL<R2_r2>>> Model1;
typedef ModelBase<SourceR3, TransformNULL<R3_r3>, 6, ModelRN<R3_r3, TransformNULL<R3_r3>>> Model2;
// typedef ModelBase<SourceSE2PosVel, TransformNULL<SE2_se2>, 5, ModelSENPosVel<SE2_se2, TransformNULL<SE2_se2>>> Model3;
// typedef ModelBase<SourceSE2PoseTwist, TransformNULL<SE2_se2>, 6, ModelSENPoseTwist<SE2_se2, TransformNULL<SE2_se2>>> Model4;
// typedef ModelBase<SourceSE3PosVel, TransformNULL<SE3_se3>, 10, ModelSENPosVel<SE3_se3, TransformNULL<SE3_se3>>> Model5;
// typedef ModelBase<SourceSE3PoseTwist, TransformNULL<SE3_se3>,12, ModelSENPoseTwist<SE3_se3, TransformNULL<SE3_se3>>> Model6;


class ModelBaseTest2 : public ::testing::Test {
public: 

void SetUp() override {

// Setup the sources
source_params1.meas_cov_fixed_ = false;
source_params1.expected_num_false_meas_ = 0.5;
source_params1.type_ = MeasurementTypes::RN_POS;
source_params1.gate_probability_ = 0.75;
source_params1.probability_of_detection_ = 0.95;
source_params1.source_index_ = 0;

source_params2.meas_cov_fixed_ = false;
source_params2.expected_num_false_meas_ = 0.2;
source_params2.type_ = MeasurementTypes::RN_POS_VEL;
source_params2.gate_probability_ = 0.8;
source_params2.probability_of_detection_ = 0.94;
source_params2.source_index_ = 1;

source_params3.meas_cov_fixed_ = false;
source_params3.expected_num_false_meas_ = 0.45;
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
model.Init(sources, params);
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
    this->model.UpdateModelLikelihood();

    double tmp1 = std::log(1+this->source1.params_.gate_probability_*this->source1.params_.probability_of_detection_*(this->update_info[0].num_assoc_meas/(source1.params_.expected_num_false_meas_*update_info[0].volume) -1));   
    double tmp2 = std::log(1+this->source2.params_.gate_probability_*this->source2.params_.probability_of_detection_*(this->update_info[1].num_assoc_meas/(source2.params_.expected_num_false_meas_*update_info[1].volume) -1));   
    double tmp3 = std::log(1+this->source3.params_.gate_probability_*this->source3.params_.probability_of_detection_*(this->update_info[2].num_assoc_meas/(source3.params_.expected_num_false_meas_*update_info[2].volume) -1));   


    ASSERT_DOUBLE_EQ(this->model.model_likelihood_,tmp1+tmp2+tmp3);

    update_info[1].num_assoc_meas = 3;
    update_info[2].num_assoc_meas = 2;
    model.model_likelihood_ = 0;
    model.model_likelihood_update_info_ = update_info;

    // Update it with associated measurements so that the likelihood increases
    for (int ii = 0; ii < 100; ++ii) {
        this->model.UpdateModelLikelihood();
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
        this->model.UpdateModelLikelihood();
    }

    ASSERT_LE(model.model_likelihood_, -50); // It should be at about -344


    // std::cout << this->model.model_likelihood_ << std::endl;

}


} // namespace rransac
