#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <math.h> 
#include <time.h>
#include <stdlib.h>

#include "lie_groups/state.h"
#include "rransac/system.h"
#include "rransac/common/models/model_RN.h"
#include "rransac/common/sources/source_R2_R3_radar.h"
#include "rransac/common/transformations/trans_radar_R2_R3_with_SE2_SE3.h"
#include "rransac/track_initialization/seed_policies/radar_R2_R3_seed_policy.h"

#include "rransac/data_containers/cluster.h"

using namespace rransac;
using namespace lie_groups;

typedef SourceRadarR2_R3<R2_r2,MeasurementTypes::R2_R3_RADAR,TransRadarR2R3WithSE2SE3> SourceR2_Radar;
typedef SourceRadarR2_R3<R2_r2,MeasurementTypes::R2_R3_RADAR_DEPTH_DERIV,TransRadarR2R3WithSE2SE3> SourceR2_RadarDeriv;

typedef SourceRadarR2_R3<R3_r3,MeasurementTypes::R2_R3_RADAR,TransRadarR2R3WithSE2SE3> SourceR3_Radar;
typedef SourceRadarR2_R3<R3_r3,MeasurementTypes::R2_R3_RADAR_DEPTH_DERIV,TransRadarR2R3WithSE2SE3> SourceR3_RadarDeriv;

typedef SourceContainer<SourceR2_Radar,SourceR2_RadarDeriv> SourceContainerR2;
typedef SourceContainer<SourceR3_Radar,SourceR3_RadarDeriv> SourceContainerR3;


template<typename _SourceContainer>
class RadarR2R3SeedPolicyTest :  public testing::Test {

public:

typedef _SourceContainer SourceContainer;
typedef ModelRN<SourceContainer> Model;
typedef typename Model::Base::Transformation Transformation;
typedef typename Model::Base::TransformDataType TransformDataType;
typedef typename Model::Base::Measurement Measurement;
typedef typename SourceContainer::Source0 Source0;
typedef typename SourceContainer::Source1 Source1;
typedef System<Model> Sys;
typedef typename Sys::ClusterT ClusterT;
typedef RadarR2R3SeedPolicy<Model> SeedPolicy;


void SetUp() override {
srand((unsigned int) time(0));
// Setup sources
source_params_.resize(SourceContainer::num_sources_);
source_params_[0].source_index_ = 0;
source_params_[0].type_ = Source0::measurement_type_;
source_params_[0].meas_cov_ = Source0::MatMeasCov::Identity()*noise_;

source_params_[1].source_index_ = 1;
source_params_[1].type_ = Source1::measurement_type_;
source_params_[1].meas_cov_ = Source1::MatMeasCov::Identity()*noise_;

// setup measurements
meas_.resize(SourceContainer::num_sources_);
for (int ii = 0; ii < SourceContainer::num_sources_; ++ii) {
    meas_[ii].type = source_params_[ii].type_;
    meas_[ii].source_index = source_params_[ii].source_index_;
    meas_[ii].transform_state = true;
    meas_[ii].transform_meas = true;
    meas_[ii].transform_data_t_m = Transformation::GetRandomTransform(10);
    meas_[ii].transform_data_m_t = meas_[ii].transform_data_t_m.inverse();
}

// Setup system
sys_params_.process_noise_covariance_ = Model::MatModelCov::Identity();
sys_.params_ = sys_params_;
sys_.source_container_.AddSource(source_params_[0]);
sys_.source_container_.AddSource(source_params_[1]);

// Setup the track
track_.state_ = Model::GetRandomState(10);

// Generate the measurements

Measurement tmp;
std::list<Measurement> meas_same_time_step;
typename ClusterT::IteratorPair iter_pair;

for (double ii = start_time_; ii < end_time_; ii += dt_) {

    if (ii != start_time_) {
        track_.PropagateModel(dt_);
    }
    meas_same_time_step.clear();

    for (int jj = 0; jj < SourceContainer::num_sources_; ++jj) {
        tmp = sys_.source_container_.GenerateRandomMeasurement(jj,sys_.source_container_.GetParams(jj).meas_cov_*0,track_.state_,meas_[jj].transform_state,meas_[jj].transform_data_t_m);
        meas_[jj].pose = tmp.pose;
        meas_[jj].twist = tmp.twist;
        meas_[jj].time_stamp = ii;
        meas_same_time_step.push_back(meas_[jj]);
    }

    measurements_.push_back(meas_same_time_step);
    sys_.current_time_ = ii;

}

for (auto outer_iter = measurements_.begin(); outer_iter != measurements_.end(); ++outer_iter) {
    for (auto inner_iter = outer_iter->begin(); inner_iter != outer_iter->end(); ++ inner_iter) {
        iter_pair.inner_it = inner_iter;
        iter_pair.outer_it = outer_iter;
        meas_subset_.push_back(iter_pair);
    }
}

std::random_shuffle(meas_subset_.begin(), meas_subset_.end());


}

std::vector<SourceParameters> source_params_;
double noise_ = 1e-2;
std::vector<Measurement> meas_;
Parameters sys_params_;
Sys sys_;
Model track_;

std::list<std::list<Measurement>> measurements_;
std::vector<typename ClusterT::IteratorPair> meas_subset_;

double start_time_ = 0;
double end_time_ = 1;
double dt_ = 0.1;

};

//---------------------------------------------------------------------------------------------

using MyTypes = ::testing::Types<SourceContainerR2,SourceContainerR3>;
TYPED_TEST_SUITE(RadarR2R3SeedPolicyTest, MyTypes);

//--------------------------------------------------------------------------------------------

TYPED_TEST(RadarR2R3SeedPolicyTest, FullTest) {

typedef RadarR2R3SeedPolicyTest<TypeParam> TestType;

double x[TestType::Model::cov_dim_];
typename TestType::Model::VecCov x_vec;
typename TestType::Model::State est_state;
typename TestType::SeedPolicy seed_policy;

seed_policy.GenerateSeedPolicy(this->meas_subset_, this->sys_, x, TestType::Model::cov_dim_);

for(int ii = 0; ii < TestType::Model::cov_dim_; ++ii) {
    x_vec(ii) = x[ii];
}

est_state.OPlusEQ(x_vec);

// std::cout << "est track pose: " << std::endl << est_state.g_.data_ << std::endl;
// std::cout << "est track twist: " << std::endl << est_state.u_.data_ << std::endl;

// std::cout << "true track pose: " << std::endl << this->track_.state_.g_.data_ << std::endl;
// std::cout << "true track twist: " << std::endl << this->track_.state_.u_.data_ << std::endl;

ASSERT_LT( (est_state.g_.data_ - this->track_.state_.g_.data_).norm(), 1e-9);
ASSERT_LT( (est_state.u_.data_ - this->track_.state_.u_.data_).norm(), 1e-9);

}