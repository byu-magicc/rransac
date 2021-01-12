#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <math.h> 

#include "system.h"
#include "common/sources/source_base.h"
#include "common/sources/source_RN.h"
#include "common/models/model_base.h"
#include "common/models/model_RN.h"
#include "common/transformations/transformation_null.h"
#include "data_containers/cluster.h"
#include "track_initialization/lmle_policies/linear_lmle_policy.h"

using namespace rransac;
using namespace lie_groups;


TEST(LINEAR_LMLE_POLICY_TEST, MainTest){

typedef ModelRN<R3_r3, TransformNULL> Model;

// Setup sources
SourceParameters source_params1;
SourceParameters source_params2;

double noise = 1e-3;

source_params1.meas_cov_fixed_ = true;
source_params1.meas_cov_ = Eigen::Matrix3d::Identity()*noise;
source_params1.type_ = MeasurementTypes::RN_POS;
source_params1.source_index_ = 0;
source_params1.expected_num_false_meas_ = 0.8;
source_params1.probability_of_detection_ = 0.8;
source_params1.gate_probability_ = 0.8;

source_params2.meas_cov_fixed_ = false;
source_params2.type_ = MeasurementTypes::RN_POS_VEL;
source_params2.source_index_ = 1;
source_params2.expected_num_false_meas_ = 0.8;
source_params2.probability_of_detection_ = 0.8;
source_params2.gate_probability_ = 0.8;

SourceR3 source1;
SourceR3 source2;

source1.Init(source_params1);
source2.Init(source_params2);

// Setup system parameters
Parameters params;
params.process_noise_covariance_ = Eigen::Matrix<double,6,6>::Identity()*noise;

// Setup system
System<Model> sys;
sys.sources_.push_back(source1);
sys.sources_.push_back(source2);
sys.params_ = params;

// Setup Measurements
Meas m1, m2;
m1.source_index = 0;
m1.type = MeasurementTypes::RN_POS;


m2.source_index = 1;
m2.type = MeasurementTypes::RN_POS_VEL;
m2.meas_cov = Eigen::Matrix<double,6,6>::Identity()*noise;

std::list<std::list<Meas>> measurements;
std::list<Meas> meas_time;
std::vector<Cluster::ConstIteratorPair> meas_subset;
Cluster::ConstIteratorPair iter_pair;

// Generate Measurements
Model x;
x.state_ = Model::State::Random();

int steps = 3;
double dt = 0.1;
double start_time = 0;

for (double ii = start_time; ii < steps*dt; ii += dt ) {
    x.PropagateModel(dt);
    meas_time.clear();
    Meas tmp1 = sys.sources_[m1.source_index].GenerateRandomMeasurement(x.state_,Eigen::Matrix3d::Identity()*sqrt(noise));
    Meas tmp2 = sys.sources_[m2.source_index].GenerateRandomMeasurement(x.state_,Eigen::Matrix<double,6,6>::Identity()*sqrt(noise));
    m1.time_stamp = ii + dt;
    m1.pose = tmp1.pose;
    m2.time_stamp = ii +dt;
    m2.pose = tmp2.pose;
    m2.twist = tmp2.twist;
    meas_time.push_back(m1);
    meas_time.push_back(m2);
    measurements.push_back(meas_time);
}

for (auto outer_iter = measurements.begin(); outer_iter != measurements.end(); ++outer_iter) {
    for (auto inner_iter = outer_iter->begin(); inner_iter != outer_iter->end(); ++ inner_iter) {
        iter_pair.inner_it = inner_iter;
        iter_pair.outer_it = outer_iter;
        meas_subset.push_back(iter_pair);
    }


}

sys.current_time_ = start_time + steps*dt;

// Generate Current state estimate
LinearLMLEPolicy<Model> policy;

typename Model::State state = policy.GenerateStateEstimatePolicy(meas_subset,sys);

ASSERT_LT( state.OMinus(x.state_).norm(), 1e-1 );

}