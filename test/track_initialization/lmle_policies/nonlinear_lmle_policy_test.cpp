#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <math.h> 


#include "state.h"
#include "system.h"
#include "common/transformations/transformation_null.h"
#include "common/models/model_SEN_pose_twist.h"
#include "common/models/model_SEN_pos_vel.h"
// #include "common/models/model_SEN_pos_vel.h"
// #include "common/sources/source_SEN_pose_twist.h"
#include "common/sources/source_base.h"
#include "data_containers/cluster.h"
#include "track_initialization/lmle_policies/nonlinear_lmle_policy.h"
#include "track_initialization/seed_policies/null_policy.h"
#include "track_initialization/seed_policies/SE2_pos_policy.h"

using namespace rransac;
using namespace lie_groups;


TEST(NONLinearPolicyTest, SE3PoseTest) {

typedef ModelSENPoseTwist<SE3_se3,TransformNULL> Model;

// Setup sources
double noise = 1e-3;
typename Model::Source source1, source2;
SourceParameters source_params1, source_params2;
source_params1.meas_cov_ = Eigen::Matrix<double,6,6>::Identity()*noise;
source_params1.source_index_ = 0;
source_params1.type_ = MeasurementTypes::SEN_POSE;
source_params2.meas_cov_ = Eigen::Matrix<double,12,12>::Identity()*noise;
source_params2.source_index_ = 1;
source_params2.type_ = MeasurementTypes::SEN_POSE_TWIST;
source1.Init(source_params1);
source2.Init(source_params2);

// Setup measurements
Meas<double> m1, m2;
m1.source_index = 0;
m1.type = MeasurementTypes::SEN_POSE;
m2.source_index = 1;
m2.type = MeasurementTypes::SEN_POSE_TWIST;
// Setup system
Parameters params;
params.process_noise_covariance_ = Eigen::Matrix<double,12,12>::Identity()*noise;
// params.nonlinear_innov_cov_id_ = true;
System<Model> sys;
sys.params_ = params;
sys.sources_.push_back(source1);
sys.sources_.push_back(source2);

// Get minimum subset;
Model track;
track.state_ = Model::State::Random();
double start_time = 0;
double end_time = 0.5;
double dt = 0.1;
Meas<double> tmp1, tmp2;
std::list<std::list<Meas<double>>> measurements;
std::list<Meas<double>> meas_time;
std::vector<Cluster<double>::IteratorPair> meas_subset;
Cluster<double>::IteratorPair iter_pair;



for (double ii = start_time; ii < end_time; ii += dt) {

    if (ii != start_time) {
        track.PropagateModel(dt);
    }
    meas_time.clear();

    tmp1 = sys.sources_[m1.source_index].GenerateRandomMeasurement(track.state_,Eigen::Matrix<double,6,6>::Identity()*sqrt(noise));
    // tmp2 = sys.sources_[m2.source_index].GenerateRandomMeasurement(track.state_,Eigen::Matrix<double,6,6>::Identity()*sqrt(noise));
    tmp2 = sys.sources_[m2.source_index].GenerateRandomMeasurement(track.state_,Eigen::Matrix<double,12,12>::Identity()*sqrt(noise));
    m1.pose = tmp1.pose;
    m1.time_stamp = ii;
    m2.pose = tmp2.pose;
    m2.twist = tmp2.twist;
    m2.time_stamp = ii;
    meas_time.push_back(m1);
    meas_time.push_back(m2);
    measurements.push_back(meas_time);
    sys.current_time_ = ii;

}

for (auto outer_iter = measurements.begin(); outer_iter != measurements.end(); ++outer_iter) {
    for (auto inner_iter = outer_iter->begin(); inner_iter != outer_iter->end(); ++ inner_iter) {
        iter_pair.inner_it = inner_iter;
        iter_pair.outer_it = outer_iter;
        meas_subset.push_back(iter_pair);
    }
}

// iter_pair.outer_it = measurements.begin();
// iter_pair.inner_it = measurements.begin()->begin();
// meas_subset.push_back(iter_pair);
// iter_pair.outer_it = std::prev(measurements.end());
// iter_pair.inner_it = std::prev(std::prev(measurements.end())->end());
// meas_subset.push_back(iter_pair);

std::random_shuffle(meas_subset.begin(), meas_subset.end());

// Generate Current state estimate
NonLinearLMLEPolicy<Model,NULLSeedPolicy> policy;
bool success;
typename Model::State state = policy.GenerateHypotheticalStateEstimatePolicy(meas_subset,sys,success);

// std::cout << " track g: " << std::endl << track.state_.g_.data_ << std::endl;
// std::cout << " track u: " << std::endl << track.state_.u_.data_ << std::endl;
// std::cout << " est g: " << std::endl << state.g_.data_ << std::endl;
// std::cout << " est u: " << std::endl << state.u_.data_ << std::endl;

ASSERT_LT( (track.state_.g_.data_-state.g_.data_).norm(), 0.1  );
ASSERT_LT( (track.state_.u_.data_-state.u_.data_).norm(), 0.1  );



}


//---------------------------------------------------------------------------------------


TEST(NONLinearPolicyTest, SE2PosTest) {

srand((unsigned int) time(0));

typedef SE2_se2 State;
typedef SourceSENPosVel<State> Source;
typedef ModelSENPosVel<State,TransformNULL> Model;

System<Model> sys;


double noise = 1e-3;
SourceParameters source_params1, source_params2;
source_params1.source_index_ = 0;
source_params1.type_ = MeasurementTypes::SEN_POS;
source_params1.meas_cov_ = Eigen::Matrix2d::Identity()*noise;

source_params2.source_index_ = 1;
source_params2.type_ = MeasurementTypes::SEN_POS_VEL;
source_params2.meas_cov_ = Eigen::Matrix4d::Identity()*noise;

Source source1, source2;
source1.Init(source_params1);
source2.Init(source_params2);

Parameters params;
params.process_noise_covariance_ = Eigen::Matrix<double,5,5>::Identity()*noise;
params.nonlinear_innov_cov_id_ = false;
params.nonlinear_LMLE_Ceres_max_num_iters_ = 10;
sys.params_ = params;
sys.sources_.push_back(source1);
sys.sources_.push_back(source2);

Meas<double> m1,m2;
m1.source_index = 0;
m1.type = MeasurementTypes::SEN_POS;
m2.source_index = 1;
m2.type = MeasurementTypes::SEN_POS_VEL;

Model track;
track.state_ = State::Random();
track.state_.g_.t_*=10;
track.state_.u_.data_(1) = 0;
track.state_.u_.data_(0) = fabs(track.state_.u_.data_(0));

double start_time = 0;
double end_time = 1;
double dt = 0.1;
Meas<double> tmp1, tmp2;
std::list<std::list<Meas<double>>> measurements;
std::list<Meas<double>> meas_time;
std::vector<Cluster<double>::IteratorPair> meas_subset;
Cluster<double>::IteratorPair iter_pair;

for (double ii = start_time; ii < end_time; ii += dt) {

    if (ii != start_time) {
        track.PropagateModel(dt);
    }
    meas_time.clear();

    tmp1 = sys.sources_[m1.source_index].GenerateRandomMeasurement(track.state_,Eigen::Matrix<double,2,2>::Identity()*sqrt(noise));
    // tmp2 = sys.sources_[m2.source_index].GenerateRandomMeasurement(track.state_,Eigen::Matrix<double,6,6>::Identity()*sqrt(noise));
    tmp2 = sys.sources_[m2.source_index].GenerateRandomMeasurement(track.state_,Eigen::Matrix<double,4,4>::Identity()*sqrt(noise));
    m1.pose = tmp1.pose;
    m1.time_stamp = ii;
    m2.pose = tmp2.pose;
    m2.twist = tmp2.twist;
    m2.time_stamp = ii;
    meas_time.push_back(m1);
    meas_time.push_back(m2);
    measurements.push_back(meas_time);
    sys.current_time_ = ii;

}

// for (auto outer_iter = measurements.begin(); outer_iter != measurements.end(); ++outer_iter) {
//     for (auto inner_iter = outer_iter->begin(); inner_iter != outer_iter->end(); ++ inner_iter) {
//         iter_pair.inner_it = inner_iter;
//         iter_pair.outer_it = outer_iter;
//         meas_subset.push_back(iter_pair);
//     }
// }

iter_pair.outer_it = measurements.begin();
iter_pair.inner_it = measurements.begin()->begin();
meas_subset.push_back(iter_pair);
iter_pair.outer_it = std::prev(measurements.end());
iter_pair.inner_it = std::prev(std::prev(measurements.end())->end());
meas_subset.push_back(iter_pair);
iter_pair.outer_it = std::prev(measurements.end(),5);
iter_pair.inner_it = std::prev(std::prev(measurements.end(),5)->end());
meas_subset.push_back(iter_pair);

// NonLinearLMLEPolicy<Model,NULLSeedPolicy> policy;
NonLinearLMLEPolicy<Model,SE2PosSeedPolicy> policy;
bool success;
typename Model::State state = policy.GenerateHypotheticalStateEstimatePolicy(meas_subset,sys,success);

ASSERT_LT( (track.state_.g_.data_-state.g_.data_).norm(), 0.2  );
ASSERT_LT( (track.state_.u_.data_-state.u_.data_).norm(), 0.3  );

// std::cout << " track g: " << std::endl << track.state_.g_.data_ << std::endl;
// std::cout << " track u: " << std::endl << track.state_.u_.data_ << std::endl;
// std::cout << " est g: " << std::endl << state.g_.data_ << std::endl;
// std::cout << " est u: " << std::endl << state.u_.data_ << std::endl;
// std::cout << "success: " << success << std::endl;


}