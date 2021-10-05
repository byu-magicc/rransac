#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <math.h> 

#include "rransac/system.h"
#include "rransac/common/sources/source_base.h"
#include "rransac/common/sources/source_RN.h"
#include "rransac/common/sources/source_container.h"
#include "rransac/common/models/model_base.h"
#include "rransac/common/models/model_RN.h"
#include "rransac/common/transformations/transformation_null.h"
#include "rransac/data_containers/cluster.h"
#include "rransac/track_initialization/lmle_policies/linear_lmle_policy.h"
#include "rransac/track_initialization/seed_policies/null_seed_policy.h"
#include "rransac/common/transformations/trans_radar_R2_R3_with_SE2_SE3.h"
#include "rransac/common/sources/source_R2_R3_radar.h"

using namespace rransac;
using namespace lie_groups;

typedef State<Rn,double,2,2> StateR2_2;
typedef State<Rn,double,2,3> StateR2_3;
typedef State<Rn,double,3,2> StateR3_2;
typedef State<Rn,double,3,3> StateR3_3;

typedef SourceRN<R3_r3,MeasurementTypes::RN_POS,TransformNULL> SourceR3_1Pos;
typedef SourceRN<R3_r3,MeasurementTypes::RN_POS_VEL,TransformNULL> SourceR3_1PosVel;

typedef SourceRN<StateR2_2,MeasurementTypes::RN_POS,TransformNULL> SourceR2_2Pos;
typedef SourceRN<StateR2_2,MeasurementTypes::RN_POS_VEL,TransformNULL> SourceR2_2PosVel;

typedef SourceRN<StateR2_3,MeasurementTypes::RN_POS,TransformNULL> SourceR2_3Pos;
typedef SourceRN<StateR2_3,MeasurementTypes::RN_POS_VEL,TransformNULL> SourceR2_3PosVel;

typedef SourceRN<StateR3_2,MeasurementTypes::RN_POS,TransformNULL> SourceR3_2Pos;
typedef SourceRN<StateR3_2,MeasurementTypes::RN_POS_VEL,TransformNULL> SourceR3_2PosVel;

typedef SourceRN<StateR3_3,MeasurementTypes::RN_POS,TransformNULL> SourceR3_3Pos;
typedef SourceRN<StateR3_3,MeasurementTypes::RN_POS_VEL,TransformNULL> SourceR3_3PosVel;

typedef SourceContainer<SourceR3_1Pos,SourceR3_1PosVel> SourceContainerR3_1;
typedef SourceContainer<SourceR2_2Pos,SourceR2_2PosVel> SourceContainerR2_2;
typedef SourceContainer<SourceR2_3Pos,SourceR2_3PosVel> SourceContainerR2_3;
typedef SourceContainer<SourceR3_2Pos,SourceR3_2PosVel> SourceContainerR3_2;
typedef SourceContainer<SourceR3_3Pos,SourceR3_3PosVel> SourceContainerR3_3;

using MyTypes = ::testing::Types<SourceContainerR3_1,SourceContainerR2_2,SourceContainerR2_3,SourceContainerR3_2,SourceContainerR3_3>;
// using MyTypes = ::testing::Types<SourceContainerR2_3>;


template<typename _SourceContainer>
class LinearLmlePolicyTest : public ::testing::Test {

public:

typedef _SourceContainer SC;
typedef typename SC::Source0 Source0;
typedef typename SC::Source1 Source1;
typedef ModelRN<SC> Model;
typedef typename Source0::Measurement Measurement;
typedef typename Source0::Base::DataType DataType;
typedef typename Source0::Base::TransformDataType TransformDataType;
typedef Cluster<DataType,TransformDataType> ClusterT;

protected:

void SetUp() override {

}

};

TYPED_TEST_SUITE(LinearLmlePolicyTest, MyTypes);

TYPED_TEST(LinearLmlePolicyTest, MainTest){
srand(time(NULL));
typedef LinearLmlePolicyTest<TypeParam> TestType;
typedef typename TestType::Measurement Measurement;

// Setup sources
SourceParameters source_params0;
SourceParameters source_params1;

double noise = 1e-4;

source_params0.meas_cov_ = TestType::Source0::MatMeasCov::Identity()*noise;
source_params0.type_ = TestType::Source0::measurement_type_;
source_params0.source_index_ = 0;
source_params0.spacial_density_of_false_meas_ = 0.8;
source_params0.probability_of_detection_ = 0.8;
source_params0.gate_probability_ = 0.8;

source_params1.type_ = TestType::Source1::measurement_type_;
source_params1.meas_cov_ = TestType::Source1::MatMeasCov::Identity()*noise;
source_params1.source_index_ = 1;
source_params1.spacial_density_of_false_meas_ = 0.8;
source_params1.probability_of_detection_ = 0.8;
source_params1.gate_probability_ = 0.8;



// Setup system parameters
Parameters params;
params.process_noise_covariance_ = TestType::Model::MatModelCov::Identity()*noise;

// Setup system
System<typename TestType::Model> sys;
sys.source_container_.AddSource(source_params0);
sys.source_container_.AddSource(source_params1);
sys.params_ = params;

// Setup Measurements
typename TestType::Measurement m1, m2;
m1.source_index = 0;
m1.type = TestType::Source0::measurement_type_;


m2.source_index = 1;
m2.type = TestType::Source1::measurement_type_;


std::list<std::list<Measurement>> measurements;
std::list<Measurement> meas_time;
std::vector<typename TestType::ClusterT::IteratorPair> meas_subset;
typename TestType::ClusterT::IteratorPair iter_pair;

// Generate Measurements
typename TestType::Model x;
x.state_ = TestType::Model::State::Random();

int steps = 5;
double dt = 0.1;
double start_time = 0;

bool transform_state = false;
typename TestType::Source0::TransformDataType EmptyMat;

for (double ii = start_time; ii < steps*dt; ii += dt ) {
    x.PropagateModel(dt);
    meas_time.clear();
    Measurement tmp1 = sys.source_container_.GenerateRandomMeasurement(m1.source_index, TestType::Source0::MatMeasCov::Identity()*sqrt(noise)*0, x.state_,transform_state,EmptyMat);
    Measurement tmp2 = sys.source_container_.GenerateRandomMeasurement(m2.source_index, TestType::Source1::MatMeasCov::Identity()*sqrt(noise)*0,x.state_,transform_state,EmptyMat);
    m1.time_stamp = ii + dt;
    m1.pose = tmp1.pose;
    m2.time_stamp = ii +dt;
    m2.pose = tmp2.pose;
    m2.twist = tmp2.twist;
    meas_time.push_back(m1);
    meas_time.push_back(m2);
    measurements.push_back(meas_time);

    // std::cout << "m1 pose: " << std::endl << m1.pose << std::endl;
    // std::cout << "m2 pose: " << std::endl << m2.pose << std::endl;
    // std::cout << "m2 twist: " << std::endl << m2.twist << std::endl;
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
LinearLMLEPolicy<typename TestType::Model, NULLSeedPolicy> policy;
bool success = false;
typename TestType::Model::State state = policy.GenerateHypotheticalStateEstimatePolicy(meas_subset,sys,success);
ASSERT_LT( state.OMinus(x.state_).norm(), 1e-10 );
// std::cout << "x " << x.state_.g_.data_ << std::endl;
// std::cout << "state" << state.g_.data_ << std::endl;

}