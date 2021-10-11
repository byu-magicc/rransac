#include <gtest/gtest.h>
#include <typeinfo>
#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>
#include <time.h>

#include "rransac/common/models/centralized_measurement_fusion.h"
#include "lie_groups/state.h"
#include "rransac/common/sources/source_RN.h"
#include "rransac/common/measurement/measurement_base.h"
#include "rransac/common/sources/source_SEN_pose_twist.h"
#include "rransac/common/sources/source_SEN_pos_vel.h"
#include "rransac/common/sources/source_container.h"
#include "rransac/common/transformations/transformation_null.h"
#include "rransac/common/models/model_RN.h"
#include "rransac/common/models/model_SEN_pose_twist.h"
#include "rransac/common/models/model_SEN_pos_vel.h"
#include "rransac/parameters.h"
#include "rransac/common/transformations/trans_homography.h"

using namespace rransac;
using namespace lie_groups;




struct TestRn {


typedef SourceRN<R2_r2, MeasurementTypes::RN_POS, TransformNULL> SourceR2Pos;
typedef SourceRN<R2_r2, MeasurementTypes::RN_POS_VEL, TransformNULL> SourceR2PosVel;
typedef SourceContainer<SourceR2Pos,SourceR2Pos,SourceR2PosVel> SC;
typedef ModelRN<SC> Model;
typedef typename Model::MatModelCov MatModelCov;
typedef typename Model::Transformation Transformation;
typedef typename Model::State State;
typedef typename Model::TransformDataType TransformDataType;

static constexpr bool transform_state_ = false;
TransformDataType transform_data_t_m_;


};

//--------------------------------------------------------------------

struct TestSE2Pose {

typedef SourceSENPoseTwist<SE2_se2, MeasurementTypes::SEN_POSE, TransformNULL> SourceSE2Pose;
typedef SourceSENPoseTwist<SE2_se2, MeasurementTypes::SEN_POSE_TWIST, TransformNULL> SourceSE2PoseTwist;
typedef SourceContainer<SourceSE2Pose,SourceSE2Pose,SourceSE2PoseTwist> SC;
typedef ModelSENPoseTwist<SC> Model;
typedef typename Model::TransformDataType TransformDataType;
static constexpr bool transform_state_ = false;
TransformDataType transform_data_t_m_;



};

//--------------------------------------------------------------------

struct TestSE2Pos {

typedef SourceSENPosVel<SE2_se2, MeasurementTypes::SEN_POS, TransformNULL> SourceSE2Pos;
typedef SourceSENPosVel<SE2_se2, MeasurementTypes::SEN_POS_VEL, TransformNULL> SourceSE2PosVel;
typedef SourceContainer<SourceSE2Pos,SourceSE2Pos,SourceSE2PosVel> SC;
typedef ModelSENPosVel<SC> Model;
typedef typename Model::TransformDataType TransformDataType;
static constexpr bool transform_state_ = false;
TransformDataType transform_data_t_m_;



};

//--------------------------------------------------------------------

struct TestSE3Pos {

typedef SourceSENPosVel<SE3_se3, MeasurementTypes::SEN_POS, TransformNULL> SourceSE3Pos;
typedef SourceSENPosVel<SE3_se3, MeasurementTypes::SEN_POS_VEL, TransformNULL> SourceSE3PosVel;
typedef SourceContainer<SourceSE3Pos,SourceSE3Pos,SourceSE3PosVel> SC;
typedef ModelSENPosVel<SC> Model;
typedef typename Model::TransformDataType TransformDataType;
static constexpr bool transform_state_ = false;
TransformDataType transform_data_t_m_;



};

template<typename T>
class CMFTest : public testing::Test {

public:
// typedef ModelRN<SC> Model;
typedef typename T::Model Model;
typedef typename Model::MatModelCov MatModelCov;
typedef typename Model::State State;
typedef typename T::SC SC;
typedef typename SC::Source0 Source0;
typedef typename SC::Source1 Source1;
typedef typename SC::Source2 Source2;
typedef typename Model::Transformation Transformation;
typedef typename Model::Measurement Measurement;
typedef typename Model::TransformDataType TransformDataType;

void SetUp() override {

// Setup sources
source_params.resize(SC::num_sources_);

source_params[0].type_ = Source0::measurement_type_;
source_params[0].source_index_ = 0;
source_params[0].meas_cov_ = Source0::MatMeasCov::Identity()*meas_noise_[0];

source_params[1].type_ = Source1::measurement_type_;
source_params[1].source_index_ = 1;
source_params[1].meas_cov_ = Source1::MatMeasCov::Identity()*meas_noise_[1];

source_params[2].type_ = Source2::measurement_type_;
source_params[2].source_index_ = 2;
source_params[2].meas_cov_ = Source2::MatMeasCov::Identity()*meas_noise_[2];

for (size_t ii = 0; ii < SC::num_sources_; ii++)
{
    source_container_.AddSource(source_params[ii]);
}

// Setup system parameters
system_params_.process_noise_covariance_ = Model::MatModelCov::Identity()*cov_noise_;


// Setup model
true_state_ = State::Random(10);
model_.Init(system_params_);
model_.err_cov_ = Model::MatModelCov::Identity()*cov_noise_;
model_.model_likelihood_ = 0.5;
model_.newest_measurement_time_stamp = 0;
model_.state_ = true_state_;
model_.OPlusEQ( model_.err_cov_.sqrt()* rransac::utilities::GaussianRandomGenerator(Model::cov_dim_) );
for (size_t ii = 0; ii < 10; ++ii) {
    model_.PropagateModel(0.1);
}

// Setup measurements
bool at_least_one_source = false;
Measurement tmp;
std::vector<Measurement> meas(SC::num_sources_);
event_weights_.resize(SC::num_sources_);
for(size_t meas_index =0; meas_index < SC::num_sources_; ++meas_index) {
    meas[meas_index].time_stamp = 0;
    meas[meas_index].source_index = meas_index;
    meas[meas_index].weight = (1.0-weight_all_measurements_false_)/(num_meas_[meas_index]);
    meas[meas_index].transform_state = test_.transform_state_;
    meas[meas_index].transform_data_t_m = test_.transform_data_t_m_;
    event_weights_[meas_index] = (1.0-weight_all_measurements_false_)/(num_meas_[meas_index]);
}

meas[0].type = Source0::measurement_type_;
meas[1].type = Source1::measurement_type_;
meas[2].type = Source2::measurement_type_;

srand(time(NULL));
int rand_num = 0;
int prob_source_produced_meas = 7;
sources_produced_measurements_.resize(SC::num_sources_,false);
typename State::Vec_SC state_error = State::OMinus(true_state_, model_.state_);
double pose_error = state_error.block(0,0,State::Group::dim_,1).norm();
double twist_error = state_error.block(State::Group::dim_,0,State::U::total_num_dim_,1).norm();
bool good_twist = true;
bool good_pose = false;
for (size_t source_index=0; source_index < SC::num_sources_; ++source_index) {

    rand_num = rand() % 10 +1;
    if (rand_num <= prob_source_produced_meas || (source_index ==2 && !at_least_one_source) ) {
        at_least_one_source = true;
        sources_produced_measurements_[source_index] = true;
        num_sources_produced_meas_++;
        num_joint_association_events_ *= ( num_meas_[source_index] + 1);
        for(size_t meas_index = 0; meas_index < num_meas_[source_index]; ++meas_index) {

            do
            {
                tmp = source_container_.GenerateRandomMeasurement(source_index,source_container_.GetParams(source_index).meas_cov_.sqrt()*0.1,true_state_,test_.transform_state_,test_.transform_data_t_m_);
            
                good_twist = true;
                good_pose =  Source0::OMinus(tmp, Source0::GetEstMeas(true_state_,test_.transform_state_,test_.transform_data_t_m_)).norm() < pose_error;
                if (source_index == 2) {
                   Measurement tmp2 = Source2::GetEstMeas(true_state_,test_.transform_state_,test_.transform_data_t_m_);
                   good_twist =  (tmp.twist - tmp2.twist).norm() < twist_error;
                }

            } while ( !good_pose && !good_twist);
            
            

            
            meas[source_index].pose = tmp.pose;
            // std::cout << tmp.pose << std::endl << std::endl;
            meas[source_index].twist = tmp.twist;
            model_.AddNewMeasurement(meas[source_index]);
        }
    }

}

MatModelCov err_cov_inv = model_.err_cov_.inverse();
std::vector<MatModelCov> HT_VRVTinv_H;
MatModelCov sum;

for (size_t source_index =0; source_index < SC::num_sources_; ++source_index ) {

    typename Model::MatH H = source_container_.GetLinObsMatState(source_index, model_.state_, test_.transform_state_, test_.transform_data_t_m_);      // Jacobian of observation function w.r.t. state
    typename Model::MatV V = source_container_.GetLinObsMatSensorNoise(source_index, model_.state_, test_.transform_state_, test_.transform_data_t_m_);                             // Jacobian of observation function w.r.t. noise
    
    HT_VRVTinv_H.push_back(H.transpose()*(V*source_container_.GetParams(source_index).meas_cov_*V.transpose()).inverse()*H);

}

sum.setZero();
association_cov_.clear();

for (int ii = 0; ii <2; ++ ii) {

    for(int jj = 0; jj<2; ++jj) {

        for (int kk=0; kk<2; ++kk) {

            if( ii != 0 && ! sources_produced_measurements_[2]) {
                continue;
            }
            if( jj != 0 && ! sources_produced_measurements_[1]) {
                continue;
            }
            if( kk != 0 && ! sources_produced_measurements_[0]) {
                continue;
            }


            if(ii > 0 && sources_produced_measurements_[2]) {
                sum += HT_VRVTinv_H[2];
            }
            if(jj > 0 && sources_produced_measurements_[1]) {
                sum += HT_VRVTinv_H[1];
            }
            if(kk > 0 && sources_produced_measurements_[0]) {
                sum += HT_VRVTinv_H[0];
            }
            association_cov_.push_back( (sum + err_cov_inv).inverse()    );
            sum.setZero();


        }
    }

}


}

std::vector<SourceParameters> source_params;
Parameters system_params_;
std::vector<double> meas_noise_ = {0.01,0.1,0.05};
double cov_noise_ = 1.5;
std::vector<unsigned int> num_meas_{5,9,7};
std::vector<double> event_weights_;
std::vector<bool> sources_produced_measurements_;
unsigned int num_sources_produced_meas_ = 0;
unsigned int num_joint_association_events_ = 1;
double weight_all_measurements_false_ = 0.05;
std::vector<MatModelCov> association_cov_;

SC source_container_;
Model model_;
State true_state_;
T test_;


};


//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------

// using MyTypes = ::testing::Types< TestRn,TestSE2Pose,TestSE2Pos,TestSE3Pos>;
using MyTypes = ::testing::Types< TestRn,TestSE2Pose>;
TYPED_TEST_SUITE(CMFTest, MyTypes);

TYPED_TEST(CMFTest, SequentialCentralizedMeasurementFusion) {

    typedef typename TypeParam::Model Model;
    typedef typename TypeParam::Model::State State;
    typedef typename TypeParam::Model::MatModelCov MatModelCov;

    MatModelCov err_cov_initial = this->model_.err_cov_;
    State state_initial = this->model_.state_;

    // std::cout << "sources produced measurements" << std::endl;
    // for (auto source : this->sources_produced_measurements_) {
    //     std::cout << source << std::endl;
    // }

    CentralizedMeasurementFusion<typename Model::ModelBase> cmf;

    bool sequential_else_parallel_fusion = true; 

    cmf.PerformCentralizedMeasurementFusion(this->source_container_, sequential_else_parallel_fusion, this->model_);


    double initial_error = State::OMinus(state_initial,this->true_state_).norm();
    double final_error = State::OMinus(this->model_.state_,this->true_state_).norm();



    // std::cout << "true state: " << std::endl << this->true_state_.g_.data_ << std::endl << this->true_state_.u_.data_ << std::endl;
    // std::cout << "initial state: " << std::endl << state_initial.g_.data_ << std::endl << state_initial.u_.data_ << std::endl;
    // std::cout << "final state: " << std::endl << this->model_.state_.g_.data_ << std::endl << this->model_.state_.u_.data_ << std::endl;

    // std::cout << "initial err cov: " << std::endl << err_cov_initial << std::endl;
    // std::cout << "final err cov: " << std::endl << this->model_.err_cov_ << std::endl;

    ASSERT_LT(final_error, initial_error);

    ASSERT_LT( this->model_.err_cov_.determinant(), err_cov_initial.determinant());


}


TYPED_TEST(CMFTest, ParallelCentralizedMeasurementFusion) {

    typedef typename TypeParam::Model Model;
    typedef typename TypeParam::Model::State State;
    typedef typename TypeParam::Model::MatModelCov MatModelCov;

    MatModelCov err_cov_initial = this->model_.err_cov_;
    State state_initial = this->model_.state_;

    CentralizedMeasurementFusion<typename Model::ModelBase> cmf;

    bool sequential_else_parallel_fusion = false; 
    cmf.PerformCentralizedMeasurementFusion(this->source_container_, sequential_else_parallel_fusion, this->model_);


    const std::vector<unsigned int>& sources_produced_measurements_ = cmf.GetSourcesProducedMeasurements();
    unsigned int num_sources_produced_measurements_ = cmf.GetNumSourcesProducedMeasurements();
    const std::vector<unsigned int>& num_single_association_events = cmf.GetNumSingleAssociationEvents();
    unsigned int num_joint_association_events = cmf.GetNumJointAssociationEvents();
    const std::vector<std::vector<double>>& single_association_events_weights = cmf.GetSingleAssociationEventsWeights();

    const std::vector<MatModelCov>& joint_association_cov = cmf.GetJointAssociationCov();





    // std::cout << "num sources produced measurements " << this->num_sources_produced_meas_ << std::endl;
    // for(auto& cov: joint_association_cov) {
    //     std::cout <<std::endl << cov << std::endl;
    // }

    // std::cout << "ture cov:" << std::endl;

    // for(auto& cov: this->association_cov_) {
    //     std::cout <<std::endl << cov << std::endl;
    // }



    for( size_t source_index = 0; source_index < Model::num_sources_; ++source_index) {

        if(this->sources_produced_measurements_[source_index]) {

            bool found = std::find(sources_produced_measurements_.begin(), sources_produced_measurements_.end(), source_index) != sources_produced_measurements_.end();
            ASSERT_TRUE(found) << "Source was not identified to having produced measurements " << std::endl;
            ASSERT_EQ(this->num_meas_[source_index] +1, num_single_association_events[source_index]);

            ASSERT_NEAR(single_association_events_weights[source_index][0], this->weight_all_measurements_false_, 1e-9);

            for(size_t association_index =1; association_index < single_association_events_weights[source_index].size(); ++ association_index) {
                ASSERT_DOUBLE_EQ(single_association_events_weights[source_index][association_index], this->event_weights_[source_index]   );
            }

        }

    }

    ASSERT_EQ(this->num_sources_produced_meas_, num_sources_produced_measurements_);
    ASSERT_EQ(this->num_joint_association_events_, num_joint_association_events);

    ASSERT_EQ(this->association_cov_.size(), joint_association_cov.size());

    for( size_t ii = 0; ii < this->association_cov_.size(); ++ii) {

        ASSERT_LT( (this->association_cov_[ii] - joint_association_cov[ii]).norm(), 1e-9);

    }

    double initial_error = State::OMinus(state_initial,this->true_state_).norm();
    double final_error = State::OMinus(this->model_.state_,this->true_state_).norm();

    // std::cout << "true state: " << std::endl << this->true_state_.g_.data_ << std::endl << this->true_state_.u_.data_ << std::endl;
    // std::cout << "initial state: " << std::endl << state_initial.g_.data_ << std::endl << state_initial.u_.data_ << std::endl;
    // std::cout << "final state: " << std::endl << this->model_.state_.g_.data_ << std::endl << this->model_.state_.u_.data_ << std::endl;

    // std::cout << "initial err cov: " << std::endl << err_cov_initial << std::endl;
    // std::cout << "final err cov: " << std::endl << this->model_.err_cov_ << std::endl;

    ASSERT_LT(final_error, initial_error);

    ASSERT_LT( this->model_.err_cov_.determinant(), err_cov_initial.determinant());


}