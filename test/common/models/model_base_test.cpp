#include <gtest/gtest.h>
#include <typeinfo>
#include <Eigen/Core>


#include "common/models/model_base.h"
#include "common/sources/source_base.h"
#include "common/sources/source_RN.h"
#include "common/sources/source_SEN_pos_vel.h"
#include "common/sources/source_SEN_pose_twist.h"
#include "parameters.h"
#include "common/transformations/transformation_base.h"
#include "common/transformations/trans_homography.h"
#include "common/transformations/transformation_null.h"


namespace rransac
{

using namespace lie_groups;

typedef ModelBase<R2_r2, R2_r2::StateType, SourceRN<R2_r2>, TransformNULL<R2_r2>> Model1;
typedef ModelBase<R3_r3, R3_r3::StateType, SourceRN<R3_r3>, TransformNULL<R3_r3>> Model2;
typedef ModelBase<SE2_se2, SE2_se2::StateType, SourceSENPosVel<SE2_se2>,    TransformNULL<SE2_se2>> Model3;
typedef ModelBase<SE2_se2, SE2_se2::StateType, SourceSENPoseTwist<SE2_se2>, TransformNULL<SE2_se2>> Model4;
typedef ModelBase<SE3_se3, SE3_se3::StateType, SourceSENPosVel<SE3_se3>,    TransformNULL<SE3_se3>> Model5;
typedef ModelBase<SE3_se3, SE3_se3::StateType, SourceSENPoseTwist<SE3_se3>, TransformNULL<SE3_se3>> Model6;

template<class M, MeasurementTypes MT1, MeasurementTypes MT2>
struct ModelHelper {
    typedef M Model;
    MeasurementTypes MeasType1 = MT1;
    MeasurementTypes MeasType2 = MT2;
};

typedef ModelHelper<Model1, MeasurementTypes::RN_POS, MeasurementTypes::RN_POS_VEL> ModelHelper1;
typedef ModelHelper<Model2, MeasurementTypes::RN_POS, MeasurementTypes::RN_POS_VEL> ModelHelper2;
typedef ModelHelper<Model3, MeasurementTypes::SEN_POS, MeasurementTypes::SEN_POS_VEL> ModelHelper3;
typedef ModelHelper<Model4, MeasurementTypes::SEN_POS, MeasurementTypes::SEN_POS_VEL> ModelHelper4;
typedef ModelHelper<Model5, MeasurementTypes::SEN_POSE, MeasurementTypes::SEN_POSE_TWIST> ModelHelper5;
typedef ModelHelper<Model6, MeasurementTypes::SEN_POSE, MeasurementTypes::SEN_POSE_TWIST> ModelHelper6;


using MyTypes = ::testing::Types<ModelHelper1, ModelHelper2, ModelHelper3, ModelHelper4, ModelHelper5, ModelHelper6 >;

template<class Model>
class ModelTest : public ::testing::Test {

public:

static constexpr unsigned int meas_dim = Model::Model::Source::Meas_dim;
static constexpr unsigned int state_dim = Model::Model::State::g_type_::dim_*2;

void SetUp() override {

// Setup the parameters
double meas_cov_scale = 0.01;
double system_cov_scale = 0.01;
source_params1.meas_cov_fixed_ = true;
source_params1.meas_cov_ = Eigen::Matrix<double,meas_dim,meas_dim>::Identity()*meas_cov_scale;
source_params1.type_ = Model::MeasType1;
source_params1.source_index_ = 0;

source_params2.meas_cov_fixed_ = false;
source_params2.type_ = Model::MeasType2;
source_params2.source_index_ = 1;


Model::Model::Source source1;
Model::Model::Source source2;
source1.Init(source_params1);
source2.Init(source_params2);

sources.push_back(source1);
sources.push_back(source2);

params.process_noise_covariance_ = Eigen::Matrix<double,state_dim,state_dim>::Identity()*system_cov_scale;

}

SourceParameters source_params1;
SourceParameters source_params2;

Parameters params;

Model::Model::State state;

std::vector<Model::Model::Source> sources;

std::vector<std::vector<Meas>> new_meas;




};


TYPED_TEST_SUITE(ModelTest, MyTypes);



TYPED_TEST(ModelTest, AllFunctions) {
int i = 0;
}

    
} // namespace rransac
