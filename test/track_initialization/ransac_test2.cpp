#include <gtest/gtest.h>
#include <stdlib.h>
#include <time.h>
#include <Eigen/Dense>
#include <string.h>
#include <chrono>
#include <vector>

#include "lie_groups/state.h"

#include "rransac/system.h"
#include "rransac/common/models/model_RN.h"
#include "rransac/common/models/model_SEN_pos_vel.h"
#include "rransac/common/models/model_SEN_pose_twist.h"
#include "rransac/common/sources/source_RN.h"
#include "rransac/common/sources/source_SEN_pos_vel.h"
#include "rransac/common/sources/source_SEN_pose_twist.h"
#include "rransac/common/measurement/measurement_base.h"
#include "rransac/parameters.h"
#include "rransac/track_initialization/seed_policies/null_seed_policy.h"
#include "rransac/track_initialization/seed_policies/SE2_pos_seed_policy.h"
#include "rransac/track_initialization/lmle_policies/linear_lmle_policy.h"
#include "rransac/track_initialization/lmle_policies/nonlinear_lmle_policy.h"
#include "rransac/common/transformations/transformation_null.h"
#include "rransac/track_initialization/ransac.h"
#include "rransac/common/data_association/validation_region_policies/validation_region_innov_policy.h"
#include "rransac/common/data_association/track_likelihood_info_policies/tli_ipdaf_policy.h"
#include "rransac/common/data_association/measurement_weight_policies/mw_ipdaf_policy.h"

using namespace lie_groups;
using namespace rransac;

struct Test1 {
    public:
    typedef ModelRN<R2_r2, TransformNULL,SourceRN> Model;
    typedef typename Model::State State;
    typedef typename State::Algebra Algebra;
    typedef typename Model::Source Source;
    typedef Ransac<Model, NULLSeedPolicy, LinearLMLEPolicy, ValidationRegionInnovPolicy, TLI_IPDAFPolicy, MW_IPDAFPolicy> RANSAC;
    typedef Eigen::Matrix<double,2,2> MatR;
    typedef Eigen::Matrix<double,4,4> MatR2;
    static constexpr MeasurementTypes MeasurementType1= MeasurementTypes::RN_POS;
    static constexpr MeasurementTypes MeasurementType2= MeasurementTypes::RN_POS_VEL;
    typedef Eigen::Matrix<double,4,4> ProcessNoiseCov;
    std::vector<State> states;
    typedef Eigen::Matrix<double,2,1> VecU;
    std::string test_name = "R2 Test";

    static State GenerateRandomState(const double fov) {
        State rand_state;
        rand_state.g_.data_ = State::Algebra::Exp(Eigen::Matrix<double,State::Group::dim_,1>::Random()*fov);
        rand_state.u_.data_ = VecU::Random()*2;
        return rand_state;
    }

    Test1() {
        double pos = 5;
        double vel = 0.5;
        states.resize(4);
        states[0].g_.data_ << pos,pos;
        states[0].u_.data_ << 0, -vel;
        states[1].g_.data_ << pos, -pos;
        states[1].u_.data_ << -vel,0;
        states[2].g_.data_ << -pos, -pos;
        states[2].u_.data_ << 0, vel;
        states[3].g_.data_ << -pos, pos;
        states[3].u_.data_ << vel,0;
    }

  
};

//---------------------------------------------------------------------------------------------------------

struct Test2 {
    public:
    typedef ModelRN<R3_r3, TransformNULL,SourceRN> Model;
    typedef typename Model::State State;
    typedef typename State::Algebra Algebra;
    typedef typename Model::Source Source;
    typedef Ransac<Model, NULLSeedPolicy, LinearLMLEPolicy, ValidationRegionInnovPolicy, TLI_IPDAFPolicy, MW_IPDAFPolicy> RANSAC;
    typedef Eigen::Matrix<double,3,3> MatR;
    typedef Eigen::Matrix<double,6,6> MatR2;
    static constexpr MeasurementTypes MeasurementType1= MeasurementTypes::RN_POS;
    static constexpr MeasurementTypes MeasurementType2= MeasurementTypes::RN_POS_VEL;
    typedef Eigen::Matrix<double,6,6> ProcessNoiseCov;
    std::vector<State> states;
    typedef Eigen::Matrix<double,3,1> VecU;
    std::string test_name = "R3 Test";

    static State GenerateRandomState(const double fov) {
        State rand_state;
        rand_state.g_.data_ = State::Algebra::Exp(Eigen::Matrix<double,State::Group::dim_,1>::Random()*fov);
        rand_state.u_.data_ = VecU::Random();
        return rand_state;
    }

    Test2() {
        double pos = 5;
        double vel = 0.5;
        states.resize(4);
        states[0].g_.data_ << pos,pos, pos;
        states[0].u_.data_ << 0, -vel, 0;
        states[1].g_.data_ << pos, -pos, -pos;
        states[1].u_.data_ << -vel,0,-vel;
        states[2].g_.data_ << -pos, -pos, -pos;
        states[2].u_.data_ << 0, vel, vel;
        states[3].g_.data_ << -pos, pos, pos;
        states[3].u_.data_ << vel,0,0;
    }

  
};

//---------------------------------------------------------------------------------------------------------

struct Test3 {
    public:
    typedef ModelSENPoseTwist<SE2_se2, TransformNULL,SourceSENPoseTwist> Model;
    typedef typename Model::State State;
    typedef typename State::Algebra Algebra;
    typedef typename Model::Source Source;
    typedef Ransac<Model, NULLSeedPolicy, NonLinearLMLEPolicy, ValidationRegionInnovPolicy, TLI_IPDAFPolicy, MW_IPDAFPolicy> RANSAC;
    typedef Eigen::Matrix<double,3,3> MatR;
    typedef Eigen::Matrix<double,6,6> MatR2;
    static constexpr MeasurementTypes MeasurementType1= MeasurementTypes::SEN_POSE;
    static constexpr MeasurementTypes MeasurementType2= MeasurementTypes::SEN_POSE_TWIST;
    typedef Eigen::Matrix<double,6,6> ProcessNoiseCov;
    std::vector<State> states;
    typedef Eigen::Matrix<double,3,1> VecU;
    std::string test_name = "SE2 Pose Test";

    static State GenerateRandomState(const double fov) {

        State rand_state;
        rand_state.g_.R_ = so2<double>::Exp(Eigen::Matrix<double,1,1>::Random()*3);
        rand_state.g_.t_ = Eigen::Matrix<double,2,1>::Random()*fov;
        rand_state.u_.data_ = VecU::Random();
        return rand_state;
    }

    Test3() {
        double pos = 5;
        double rot = 0.1;
        double t_vel = 0.5;
        double a_vel = 0.1;
        Eigen::Matrix<double,3,1> pose;
        states.resize(4);
        pose << pos, pos, 0;
        states[0].g_.data_ = State::Algebra::Exp(pose);
        states[0].u_.data_ << t_vel, t_vel,0;
        pose << pos, -pos, 0;
        states[1].g_.data_ = State::Algebra::Exp(pose);
        states[1].u_.data_ << -t_vel,-t_vel,0;
        pose << -pos, -pos, rot;
        states[2].g_.data_ = State::Algebra::Exp(pose);
        states[2].u_.data_ << t_vel, 0, a_vel;
        pose << -pos, pos, -rot;
        states[3].g_.data_ = State::Algebra::Exp(pose);
        states[3].u_.data_ << t_vel,0,-a_vel;
    }

  
};

//---------------------------------------------------------------------------------------------------------

struct Test4 {
    public:
    typedef ModelSENPoseTwist<SE3_se3, TransformNULL,SourceSENPoseTwist> Model;
    typedef typename Model::State State;
    typedef typename State::Algebra Algebra;
    typedef typename Model::Source Source;
    typedef Ransac<Model, NULLSeedPolicy, NonLinearLMLEPolicy, ValidationRegionInnovPolicy, TLI_IPDAFPolicy, MW_IPDAFPolicy> RANSAC;
    typedef Eigen::Matrix<double,6,6> MatR;
    typedef Eigen::Matrix<double,12,12> MatR2;
    static constexpr MeasurementTypes MeasurementType1= MeasurementTypes::SEN_POSE;
    static constexpr MeasurementTypes MeasurementType2= MeasurementTypes::SEN_POSE_TWIST;
    typedef Eigen::Matrix<double,12,12> ProcessNoiseCov;
    std::vector<State> states;
    typedef Eigen::Matrix<double,6,1> VecU;
    std::string test_name = "SE3 Pose Test";

    static State GenerateRandomState(const double fov) {

        State rand_state;
        rand_state.g_.R_ = so3<double>::Exp(Eigen::Matrix<double,3,1>::Random()*3);
        rand_state.g_.t_ = Eigen::Matrix<double,3,1>::Random()*fov;
        rand_state.u_.data_ = VecU::Random()*2;
        return rand_state;
    }

    Test4() {
        double pos = 10;
        double rot1 = 0.2;
        double rot2 = -0.2;
        double rot3 = 0;
        double t_vel = 0.5;
        double a_vel = 0.1;
        Eigen::Matrix<double,6,1> pose;
        states.resize(4);
        pose << pos, pos, pos, rot1, rot2, rot3;
        states[0].g_.data_ = State::Algebra::Exp(pose);
        states[0].u_.data_ << t_vel, t_vel, t_vel, 0, 0, 0;
        pose << -pos, pos, -pos, rot1, -rot2, rot3;
        states[1].g_.data_ = State::Algebra::Exp(pose);
        states[1].u_.data_ << -t_vel, t_vel, t_vel, a_vel, -a_vel, a_vel;
        pose << pos, -pos, pos, -rot1, rot2, -rot3;
        states[2].g_.data_ = State::Algebra::Exp(pose);
        states[2].u_.data_ << t_vel, -t_vel, t_vel, -a_vel, a_vel, -a_vel;
        pose << -pos, -pos, -pos, -rot1, -rot2, -rot3;
        states[3].g_.data_ = State::Algebra::Exp(pose);
        states[3].u_.data_ << -t_vel, -t_vel, -t_vel, 0, 0, 0;
    }

  
};

//---------------------------------------------------------------------------------------------------------

struct Test5 {
    public:
    typedef ModelSENPosVel<SE2_se2, TransformNULL,SourceSENPosVel> Model;
    typedef typename Model::State State;
    typedef typename State::Algebra Algebra;
    typedef typename Model::Source Source;
    typedef Ransac<Model, SE2PosSeedPolicy, NonLinearLMLEPolicy, ValidationRegionInnovPolicy, TLI_IPDAFPolicy, MW_IPDAFPolicy> RANSAC;
    typedef Eigen::Matrix<double,2,2> MatR;
    typedef Eigen::Matrix<double,4,4> MatR2;
    static constexpr MeasurementTypes MeasurementType1= MeasurementTypes::SEN_POS;
    static constexpr MeasurementTypes MeasurementType2= MeasurementTypes::SEN_POS_VEL;
    typedef Eigen::Matrix<double,5,5> ProcessNoiseCov;
    std::vector<State> states;
    typedef Eigen::Matrix<double,3,1> VecU;
    std::string test_name = "SE2 Pos Test";

    static State GenerateRandomState(const double fov) {

        State rand_state;
        rand_state.g_.R_ = so2<double>::Exp(Eigen::Matrix<double,1,1>::Random()*3);
        rand_state.g_.t_ = Eigen::Matrix<double,2,1>::Random()*fov;
        rand_state.u_.data_ = VecU::Random()*2;
        return rand_state;
    }

    Test5() {
        double pos = 5;
        double rot = 0.1;
        double t_vel = 0.5;
        double a_vel = 0.1;
        Eigen::Matrix<double,3,1> pose;
        states.resize(4);
        pose << pos, pos, 0;
        states[0].g_.data_ = State::Algebra::Exp(pose);
        states[0].u_.data_ << t_vel, 0,0;
        pose << pos, -pos, rot*2;
        states[1].g_.data_ = State::Algebra::Exp(pose);
        states[1].u_.data_ << t_vel,0,0;
        pose << -pos, -pos, rot;
        states[2].g_.data_ = State::Algebra::Exp(pose);
        states[2].u_.data_ << t_vel, 0, a_vel;
        pose << -pos, pos, -rot;
        states[3].g_.data_ = State::Algebra::Exp(pose);
        states[3].u_.data_ << t_vel,0,-a_vel;
    }

  
};

//---------------------------------------------------------------------------------------------------------

// This test needs a good seed policy to work 

// struct Test6 {
//     public:
//     typedef ModelSENPosVel<SE3_se3, TransformNULL> Model;
//     typedef typename Model::State State;
//     typedef typename State::Algebra Algebra;
//     typedef typename Model::Source Source;
//     typedef Ransac<Model, NULLSeedPolicy, NonLinearLMLEPolicy, ModelPDFPolicy> RANSAC;
//     typedef Eigen::Matrix<double,3,3> MatR;
//     typedef Eigen::Matrix<double,6,6> MatR2;
//     static constexpr MeasurementTypes MeasurementType1= MeasurementTypes::SEN_POS;
//     static constexpr MeasurementTypes MeasurementType2= MeasurementTypes::SEN_POS_VEL;
//     typedef Eigen::Matrix<double,10,10> ProcessNoiseCov;
//     std::vector<State> states;
//     typedef Eigen::Matrix<double,6,1> VecU;
//     std::string test_name = "SE3 Pos Test";

//     Test6() {
//         double pos = 10;
//         double rot1 = 0;
//         double rot2 = 0.2;
//         double rot3 = 0;
//         double t_vel = 0.5;
//         double a_vel = 0.1;
//         Eigen::Matrix<double,6,1> pose;
//         states.resize(4);
//         pose << pos, pos, pos, rot1, rot2, rot3;
//         states[0].g_.data_ = State::Algebra::Exp(pose);
//         states[0].u_.data_ << t_vel, 0, 0, 0, 0, 0;
//         pose << -pos, pos, -pos, rot1, -rot2, rot3;
//         states[1].g_.data_ = State::Algebra::Exp(pose);
//         states[1].u_.data_ << t_vel, 0, 0, 0, 0, 0;
//         pose << pos, -pos, pos, -rot1, rot2, -rot3;
//         states[2].g_.data_ = State::Algebra::Exp(pose);
//         states[2].u_.data_ << t_vel, 0, 0, 0, -a_vel, 0;
//         pose << -pos, -pos, -pos, -rot1, -rot2, -rot3;
//         states[3].g_.data_ = State::Algebra::Exp(pose);
//         states[3].u_.data_ << t_vel, 0, 0, 0, 0, -a_vel;
//     }

  
// };

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------

template<typename T> 
class RANSACTest : public testing::Test {

public:

typedef typename T::Model Model;
typedef typename T::State State;
typedef typename T::Source Source;
typedef typename T::RANSAC RANSAC;

void SetUp() override {

SourceParameters source_params1, source_params2, source_params3;
source_params1.type_ = T::MeasurementType1;
source_params1.source_index_ = 0;
source_params1.meas_cov_ = T::MatR::Identity()*noise;
source_params1.gate_probability_ = 0.9;
source_params1.probability_of_detection_ = 0.85;

source_params2.type_ = T::MeasurementType2;
source_params2.source_index_ = 1;
source_params2.meas_cov_ = T::MatR2::Identity()*noise;
source_params2.gate_probability_ = 0.9;
source_params2.probability_of_detection_ = 0.85;


source_params3.type_ = T::MeasurementType2;
source_params3.source_index_ = 2;
source_params3.meas_cov_ = T::MatR2::Identity()*noise;
source_params3.gate_probability_ = 0.9;
source_params3.probability_of_detection_ = 0.85;


Source source1,source2, source3;
source1.Init(source_params1);
source2.Init(source_params2);
source3.Init(source_params3);

// Setup system
Parameters params;
params.process_noise_covariance_ = T::ProcessNoiseCov::Identity()*noise;
params.RANSAC_max_iters_ = 100;
params.RANSAC_minimum_subset_ = 3;
params.RANSAC_score_stopping_criteria_ = 10;
params.RANSAC_score_minimum_requirement_ = 6;
params.meas_time_window_ = 5;                   // 5 seconds
params.cluster_time_threshold_ = 2;
params.cluster_velocity_threshold_ = 1;
params.cluster_position_threshold_ = 0.5;
params.track_max_num_tracks_ = 5;
// params.nonlinear_innov_cov_id_ = true;

sys.sources_.push_back(source1);
sys.sources_.push_back(source2);
sys.sources_.push_back(source3);
sys.params_ = params;

// Setup Measurements
Meas<double> m1, m2, m3, m4;
m1.source_index = 0;
m1.type = source_params1.type_ ;


m2.source_index = 1;
m2.type = source_params2.type_;

// This measurement is noise
m3.source_index = 2;
m3.type = source_params3.type_;

m4.source_index = 0;
m4.type =source_params1.type_ ;

// Setup tracks
tracks.resize(test_data.states.size());
for (int ii = 0; ii < test_data.states.size(); ++ii) {
    tracks[ii].Init(sys.params_,sys.sources_.size());
    tracks[ii].state_ = test_data.states[ii];
}

// Create simulation data
double dt = 0.1;
double end_time = 5; // seconds;
double start_time = 0; // seconds;
double fov = 50;  // The surveillance region is a square centered at zero with side length 20
Meas<double> tmp1, tmp2, tmp3, tmp4;

for (double ii =start_time; ii < end_time; ii += dt) {

    for (auto& track: tracks) {

        if (ii !=start_time) {
            track.PropagateModel(dt);
        }

        State rand_state = T::GenerateRandomState(fov);

        tmp1 = sys.sources_[m1.source_index].GenerateRandomMeasurement(track.state_,T::MatR::Identity()*sqrt(noise)*0);
        tmp2 = sys.sources_[m2.source_index].GenerateRandomMeasurement(track.state_,T::MatR2::Identity()*sqrt(noise)*0);
        tmp3 = sys.sources_[m3.source_index].GenerateRandomMeasurement(rand_state,T::MatR2::Identity()*sqrt(noise)*0);
        tmp4 = sys.sources_[m4.source_index].GenerateRandomMeasurement(track.state_,T::MatR::Identity()*sqrt(noise)*0);

        m1.time_stamp = ii;
        m1.pose = tmp1.pose;
        m2.time_stamp = ii;
        m2.pose = tmp2.pose;
        m2.twist = tmp2.twist;
        m3.time_stamp = ii;
        m3.pose = tmp3.pose;
        m3.twist = tmp3.twist;
        m4.time_stamp = ii;
        m4.pose = tmp4.pose;
        sys.current_time_ = ii;

        sys.data_tree_.AddMeasurement(sys, m1);
        sys.data_tree_.AddMeasurement(sys, m2);
        sys.data_tree_.AddMeasurement(sys, m3);
        sys.data_tree_.AddMeasurement(sys, m4);
    }
}

sys.data_tree_.ConstructClusters(sys);
sys.dt_ = dt;


}

double noise = 1e-2;
System<Model> sys;
T test_data;
std::vector<Model> tracks;
RANSAC ransac;



};

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------

// using MyTypes = ::testing::Types<Test1,Test2,Test3,Test4,Test5>;
using MyTypes = ::testing::Types<Test1>;
TYPED_TEST_SUITE(RANSACTest, MyTypes);

TYPED_TEST(RANSACTest, FullTest) {

    auto start = std::chrono::system_clock::now();
    this->ransac.Run(this->sys);
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()/1e9;
    std::cout << "RANSAC test benchmark for " << this->test_data.test_name << " : " <<  elapsed << '\n';


// for (auto& created_track: this->sys.models_) {
//     std::cout << "created_track g: " << std::endl << created_track.state_.g_.data_ << std::endl;
//     std::cout << "created_track u: " << std::endl << created_track.state_.u_.data_ << std::endl << std::endl;

// }

// for (auto& sim_track: this->tracks) {
//     std::cout << "sim_track g: " << std::endl << sim_track.state_.g_.data_ << std::endl;
//     std::cout << "sim_track u: " << std::endl << sim_track.state_.u_.data_ << std::endl << std::endl;

// }


// make sure that the tracks were created
ASSERT_GE(this->sys.models_.size(), 4 );
for (auto& created_track: this->sys.models_) {

    bool found = false;
    for (auto& sim_track: this->tracks) {

        if (sim_track.state_.OMinus(created_track.state_).norm() < 5e-1) {
            found = true;
            ASSERT_LT(created_track.err_cov_.norm(), 1);
        }
        
    }
    ASSERT_TRUE(found);

}

// Make sure that the consensus set is not empty and that the model likelihood is greater than 0.5.
for (auto& track : this->sys.models_) {
    ASSERT_GE(track.cs_.Size(), this->sys.params_.RANSAC_score_minimum_requirement_);
    ASSERT_GE(track.model_likelihood_, 0.5);
}


// make sure that all of the measurement are removed
for (auto track_iter = this->sys.models_.begin(); track_iter != this->sys.models_.end(); ++track_iter) {
for (auto track_meas_outer = track_iter->cs_.consensus_set_.begin(); track_meas_outer != track_iter->cs_.consensus_set_.end(); ++track_meas_outer) {
for (auto track_meas_inner = track_meas_outer->begin(); track_meas_inner != track_meas_outer->end(); ++ track_meas_inner) {
for (auto cluster_iter = this->sys.data_tree_.data_.begin(); cluster_iter != this->sys.data_tree_.data_.end(); ++cluster_iter) {
for (auto cluster_meas_outer = cluster_iter->data_.begin(); cluster_meas_outer != cluster_iter->data_.end(); ++ cluster_meas_outer) {
for (auto cluster_meas_inner = cluster_meas_outer->begin(); cluster_meas_inner != cluster_meas_outer->end(); ++ cluster_meas_inner) {

    
    if (cluster_meas_inner->source_index == track_meas_inner->source_index && cluster_meas_inner->time_stamp == track_meas_inner->time_stamp) {
        
        ASSERT_GT(this->sys.sources_[cluster_meas_inner->source_index].OMinus(*cluster_meas_inner, *track_meas_inner).norm(), 1e-8 );
        
    }


}}}}}}


}