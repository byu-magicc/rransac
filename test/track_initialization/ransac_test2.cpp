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
#include "rransac/common/sources/source_R2_R3_radar.h"
#include "rransac/common/measurement/measurement_base.h"
#include "rransac/parameters.h"
#include "rransac/track_initialization/seed_policies/null_seed_policy.h"
#include "rransac/track_initialization/seed_policies/SE2_pos_seed_policy.h"
#include "rransac/track_initialization/lmle_policies/linear_lmle_policy.h"
#include "rransac/track_initialization/lmle_policies/nonlinear_lmle_policy.h"
#include "rransac/common/transformations/transformation_null.h"
#include "rransac/common/transformations/trans_radar_R2_R3_with_SE2_SE3.h"
#include "rransac/track_initialization/ransac.h"
#include "rransac/common/data_association/validation_region_policies/validation_region_innov_policy.h"
#include "rransac/common/data_association/track_likelihood_info_policies/tli_ipdaf_policy.h"
#include "rransac/common/data_association/measurement_weight_policies/mw_ipdaf_policy.h"

using namespace lie_groups;
using namespace rransac;

struct Test1 {
    public:

    typedef SourceRN<R2_r2,MeasurementTypes::RN_POS,TransformNULL> SourceR2Pos;
    typedef SourceRN<R2_r2,MeasurementTypes::RN_POS_VEL,TransformNULL> SourceR2PosVEL;
    typedef SourceContainer<SourceR2Pos,SourceR2PosVEL,SourceR2PosVEL> SC;
    typedef ModelRN<SC> Model;
    typedef typename Model::State State;
    typedef typename State::Algebra Algebra;  
    typedef Ransac<Model, NULLSeedPolicy, LinearLMLEPolicy, ValidationRegionInnovPolicy, TLI_IPDAFPolicy, MW_IPDAFPolicy> RANSAC;
    typedef typename Model::Base::Measurement Measurement;
    typedef typename Model::Base::Transformation Transformation;
    typedef typename Model::Base::TransformDataType TransformDataType;

    bool transform_state0_ = false;
    bool transform_state1_ = false;
    bool transform_state2_ = false;
    bool transform_measurement0_ = false;
    bool transform_measurement1_ = false;
    bool transform_measurement2_ = false;
    TransformDataType transform_data_t_m_0_;
    TransformDataType transform_data_t_m_1_;
    TransformDataType transform_data_t_m_2_;
    TransformDataType transform_data_m_t_0_;
    TransformDataType transform_data_m_t_1_;
    TransformDataType transform_data_m_t_2_;


    std::vector<State> states;
    std::string test_name = "R2 Test";

  

    Test1() {
        double pos = 10;
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

    typedef SourceRN<R3_r3,MeasurementTypes::RN_POS,TransformNULL> SourceRnPos;
    typedef SourceRN<R3_r3,MeasurementTypes::RN_POS_VEL,TransformNULL> SourceRnPosVEL;
    typedef SourceContainer<SourceRnPos,SourceRnPosVEL,SourceRnPosVEL> SC;

    typedef ModelRN<SC> Model;
    typedef typename Model::State State;
    typedef typename State::Algebra Algebra;
    typedef Ransac<Model, NULLSeedPolicy, LinearLMLEPolicy, ValidationRegionInnovPolicy, TLI_IPDAFPolicy, MW_IPDAFPolicy> RANSAC;

    typedef typename Model::Base::Measurement Measurement;
    typedef typename Model::Base::Transformation Transformation;
    typedef typename Model::Base::TransformDataType TransformDataType;

    bool transform_state0_ = false;
    bool transform_state1_ = false;
    bool transform_state2_ = false;
    bool transform_measurement0_ = false;
    bool transform_measurement1_ = false;
    bool transform_measurement2_ = false;
    TransformDataType transform_data_t_m_0_;
    TransformDataType transform_data_t_m_1_;
    TransformDataType transform_data_t_m_2_;
    TransformDataType transform_data_m_t_0_;
    TransformDataType transform_data_m_t_1_;
    TransformDataType transform_data_m_t_2_;

    std::vector<State> states;
    std::string test_name = "R3 Test";





    Test2() {
        double pos = 10;
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

    typedef SourceSENPoseTwist<SE2_se2,MeasurementTypes::SEN_POSE,TransformNULL> SourceSE2Pose;
    typedef SourceSENPoseTwist<SE2_se2,MeasurementTypes::SEN_POSE_TWIST,TransformNULL> SourceSE2PoseTwist;
    typedef SourceContainer<SourceSE2Pose,SourceSE2PoseTwist,SourceSE2PoseTwist> SC;

    typedef ModelSENPoseTwist<SC> Model;
    typedef typename Model::State State;
    typedef typename State::Algebra Algebra;
    typedef Ransac<Model, NULLSeedPolicy, NonLinearLMLEPolicy, ValidationRegionInnovPolicy, TLI_IPDAFPolicy, MW_IPDAFPolicy> RANSAC;
    typedef typename Model::Base::Measurement Measurement;
    typedef typename Model::Base::Transformation Transformation;
    typedef typename Model::Base::TransformDataType TransformDataType;

    bool transform_state0_ = false;
    bool transform_state1_ = false;
    bool transform_state2_ = false;
    bool transform_measurement0_ = false;
    bool transform_measurement1_ = false;
    bool transform_measurement2_ = false;
    TransformDataType transform_data_t_m_0_;
    TransformDataType transform_data_t_m_1_;
    TransformDataType transform_data_t_m_2_;
    TransformDataType transform_data_m_t_0_;
    TransformDataType transform_data_m_t_1_;
    TransformDataType transform_data_m_t_2_;
    std::vector<State> states;
    std::string test_name = "SE2 Pose Test";


    Test3() {
        double pos = 10;
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
    typedef SourceSENPoseTwist<SE3_se3,MeasurementTypes::SEN_POSE,TransformNULL> SourceSE3Pose;
    typedef SourceSENPoseTwist<SE3_se3,MeasurementTypes::SEN_POSE_TWIST,TransformNULL> SourceSE3PoseTwist;
    typedef SourceContainer<SourceSE3Pose,SourceSE3PoseTwist,SourceSE3PoseTwist> SC;

    typedef ModelSENPoseTwist<SC> Model;
    typedef typename Model::State State;
    typedef typename State::Algebra Algebra;
    typedef Ransac<Model, NULLSeedPolicy, NonLinearLMLEPolicy, ValidationRegionInnovPolicy, TLI_IPDAFPolicy, MW_IPDAFPolicy> RANSAC;
    typedef typename Model::Base::Measurement Measurement;
    typedef typename Model::Base::Transformation Transformation;
    typedef typename Model::Base::TransformDataType TransformDataType;

    bool transform_state0_ = false;
    bool transform_state1_ = false;
    bool transform_state2_ = false;
    bool transform_measurement0_ = false;
    bool transform_measurement1_ = false;
    bool transform_measurement2_ = false;
    TransformDataType transform_data_t_m_0_;
    TransformDataType transform_data_t_m_1_;
    TransformDataType transform_data_t_m_2_;
    TransformDataType transform_data_m_t_0_;
    TransformDataType transform_data_m_t_1_;
    TransformDataType transform_data_m_t_2_;
    std::vector<State> states;
    std::string test_name = "SE3 Pose Test";


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
    typedef SourceSENPosVel<SE2_se2,MeasurementTypes::SEN_POS,TransformNULL> SourceSE2Pos;
    typedef SourceSENPosVel<SE2_se2,MeasurementTypes::SEN_POS_VEL,TransformNULL> SourceSE2PosVel;
    typedef SourceContainer<SourceSE2Pos,SourceSE2PosVel,SourceSE2PosVel> SC;

    typedef ModelSENPosVel<SC> Model;
    typedef typename Model::State State;
    typedef typename State::Algebra Algebra;
    typedef Ransac<Model, SE2PosSeedPolicy, NonLinearLMLEPolicy, ValidationRegionInnovPolicy, TLI_IPDAFPolicy, MW_IPDAFPolicy> RANSAC;
    typedef typename Model::Base::Measurement Measurement;
    typedef typename Model::Base::Transformation Transformation;
    typedef typename Model::Base::TransformDataType TransformDataType;

    bool transform_state0_ = false;
    bool transform_state1_ = false;
    bool transform_state2_ = false;
    bool transform_measurement0_ = false;
    bool transform_measurement1_ = false;
    bool transform_measurement2_ = false;
    TransformDataType transform_data_t_m_0_;
    TransformDataType transform_data_t_m_1_;
    TransformDataType transform_data_t_m_2_;
    TransformDataType transform_data_m_t_0_;
    TransformDataType transform_data_m_t_1_;
    TransformDataType transform_data_m_t_2_;
    std::vector<State> states;
    std::string test_name = "SE2 Pos Test";



    Test5() {
        double pos = 10;
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



//---------------------------------------------------------------------------------------------------------

struct Test7 {
    public:
    typedef lie_groups::State<Rn,double, 3,2> StateR3_2;
    typedef SourceRadarR2_R3<StateR3_2,MeasurementTypes::R2_R3_RADAR,TransRadarR2R3WithSE2SE3> SourceRadar;
    typedef SourceRadarR2_R3<StateR3_2,MeasurementTypes::R2_R3_RADAR_DEPTH_DERIV,TransRadarR2R3WithSE2SE3> SourceRadarDepth;
    typedef SourceContainer<SourceRadar,SourceRadarDepth,SourceRadarDepth> SC;

    typedef ModelRN<SC> Model;
    typedef typename Model::State State;
    typedef typename State::Algebra Algebra;
    typedef Ransac<Model, NULLSeedPolicy, NonLinearLMLEPolicy, ValidationRegionInnovPolicy, TLI_IPDAFPolicy, MW_IPDAFPolicy> RANSAC;

    typedef typename Model::Base::Measurement Measurement;
    typedef typename Model::Base::Transformation Transformation;
    typedef typename Model::Base::TransformDataType TransformDataType;

    bool transform_state0_ = true;
    bool transform_state1_ = true;
    bool transform_state2_ = false;
    bool transform_measurement0_ = true;
    bool transform_measurement1_ = true;
    bool transform_measurement2_ = false;
    TransformDataType transform_data_t_m_0_;
    TransformDataType transform_data_t_m_1_;
    TransformDataType transform_data_t_m_2_;
    TransformDataType transform_data_m_t_0_;
    TransformDataType transform_data_m_t_1_;
    TransformDataType transform_data_m_t_2_;

    std::vector<State> states;
    std::string test_name = "R3 Test";





    Test7() {
        transform_data_t_m_0_ = Transformation::GetRandomTransform(10);
        transform_data_t_m_1_ = Transformation::GetRandomTransform(10);
        transform_data_t_m_2_.setIdentity();
        transform_data_m_t_0_ = transform_data_t_m_0_.inverse();
        transform_data_m_t_1_ = transform_data_t_m_1_.inverse();
        transform_data_m_t_2_ = transform_data_t_m_2_.inverse();
        double pos = 20;
        double vel = 0.5;
        double accel = 0;
        states.resize(4);
        states[0].g_.data_ << pos,pos, pos;
        states[0].u_.data_ << 0, -vel, 0, 0, 0, accel;
        states[1].g_.data_ << 0, -pos, -pos;
        states[1].u_.data_ << -vel,0,-vel, 0, 0, -accel;
        states[2].g_.data_ << -pos, -pos, -pos*2;
        states[2].u_.data_ << 0, vel, vel, accel, 0, 0;
        states[3].g_.data_ << -pos, pos, pos;
        states[3].u_.data_ << vel,0,0, -accel, 0, 0;

    }

  
};

//---------------------------------------------------------------------------------------------------------

struct Test8 {
    public:
    typedef lie_groups::State<Rn,double, 2,2> StateR2_2;
    typedef SourceRadarR2_R3<StateR2_2,MeasurementTypes::R2_R3_RADAR,TransRadarR2R3WithSE2SE3> SourceRadar;
    typedef SourceRadarR2_R3<StateR2_2,MeasurementTypes::R2_R3_RADAR_DEPTH_DERIV,TransRadarR2R3WithSE2SE3> SourceRadarDepth;
    typedef SourceContainer<SourceRadar,SourceRadarDepth,SourceRadarDepth> SC;

    typedef ModelRN<SC> Model;
    typedef typename Model::State State;
    typedef typename State::Algebra Algebra;
    typedef Ransac<Model, NULLSeedPolicy, NonLinearLMLEPolicy, ValidationRegionInnovPolicy, TLI_IPDAFPolicy, MW_IPDAFPolicy> RANSAC;

    typedef typename Model::Base::Measurement Measurement;
    typedef typename Model::Base::Transformation Transformation;
    typedef typename Model::Base::TransformDataType TransformDataType;

    bool transform_state0_ = true;
    bool transform_state1_ = true;
    bool transform_state2_ = false;
    bool transform_measurement0_ = true;
    bool transform_measurement1_ = true;
    bool transform_measurement2_ = false;
    TransformDataType transform_data_t_m_0_;
    TransformDataType transform_data_t_m_1_;
    TransformDataType transform_data_t_m_2_;
    TransformDataType transform_data_m_t_0_;
    TransformDataType transform_data_m_t_1_;
    TransformDataType transform_data_m_t_2_;

    std::vector<State> states;
    std::string test_name = "R2 radar Test";





    Test8() {
        transform_data_t_m_0_ = Transformation::GetRandomTransform(10);
        transform_data_t_m_1_ = Transformation::GetRandomTransform(10);
        // transform_data_t_m_0_.setIdentity();
        // transform_data_t_m_1_.setIdentity();
        transform_data_t_m_2_.setIdentity();
        transform_data_m_t_0_ = transform_data_t_m_0_.inverse();
        transform_data_m_t_1_ = transform_data_t_m_1_.inverse();
        transform_data_m_t_2_ = transform_data_t_m_2_.inverse();
        double pos = 20;
        double vel = 0.5;
        double accel = 0;
        states.resize(4);
        states[0].g_.data_ << pos,0;
        states[0].u_.data_ << 0, -vel, accel, 0;
        states[1].g_.data_ << pos, -pos/2;
        states[1].u_.data_ << -vel,0, 0, -accel;
        states[2].g_.data_ << -pos*2, -pos;
        states[2].u_.data_ << 0, vel, 0, accel;
        states[3].g_.data_ << -pos, pos;
        states[3].u_.data_ << vel,0, accel,accel;
    }

  
};

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------

template<typename T> 
class RANSACTest : public testing::Test {

public:

typedef typename T::Model Model;
typedef typename T::State State;
typedef typename T::RANSAC RANSAC;
typedef typename T::Measurement Measurement;
typedef typename T::Transformation Transformation;
typedef typename T::TransformDataType TransformDataType;
typedef typename T::SC::Source0 Source0;
typedef typename T::SC::Source1 Source1;
typedef typename T::SC::Source2 Source2;




void SetUp() override {

SourceParameters source_params0, source_params1, source_params2;
source_params0.type_ = Source0::measurement_type_;
source_params0.source_index_ = 0;
source_params0.meas_cov_ = Source0::MatMeasCov::Identity()*noise;
source_params0.gate_probability_ = 0.9;
source_params0.probability_of_detection_ = 0.85;

source_params1.type_ = Source1::measurement_type_;
source_params1.source_index_ = 1;
source_params1.meas_cov_ = Source1::MatMeasCov::Identity()*noise;
source_params1.gate_probability_ = 0.9;
source_params1.probability_of_detection_ = 0.85;


source_params2.type_ = Source2::measurement_type_;
source_params2.source_index_ = 2;
source_params2.meas_cov_ = Source2::MatMeasCov::Identity()*noise;
source_params2.gate_probability_ = 0.9;
source_params2.probability_of_detection_ = 0.85;



// Setup system
Parameters params;
params.process_noise_covariance_ = Model::Base::MatModelCov::Identity()*noise;
params.RANSAC_max_iters_ = 100;
params.RANSAC_minimum_subset_ = 5;
params.RANSAC_score_stopping_criteria_ = 10;
params.RANSAC_score_minimum_requirement_ = 6;
params.meas_time_window_ = 5;                   // 5 seconds
params.cluster_time_threshold_ = 1;
params.cluster_velocity_threshold_ = 1.5;
params.cluster_position_threshold_ = 0.6;
params.track_max_num_tracks_ = 5;
// params.nonlinear_innov_cov_id_ = true;

sys.source_container_.AddSource(source_params0);
sys.source_container_.AddSource(source_params1);
sys.source_container_.AddSource(source_params2);
sys.params_ = params;

// Setup Measurements
Measurement m0, m1, m2, m3;
m0.source_index = 0;
m0.type = source_params0.type_ ;


m1.source_index = 1;
m1.type = source_params1.type_;

// This measurement is noise
m2.source_index = 2;
m2.type = source_params2.type_;

m3.source_index = 0;
m3.type =source_params0.type_ ;

// Setup tracks
tracks.resize(test_data.states.size());
for (int ii = 0; ii < test_data.states.size(); ++ii) {
    tracks[ii].Init(sys.params_);
    tracks[ii].state_ = test_data.states[ii];
}

// Create simulation data

double dt = 0.1;
double end_time = 5; // seconds;
double start_time = 0; // seconds;
double fov = 50;  // The surveillance region is a square centered at zero with side length 20
Measurement tmp0, tmp1, tmp2, tmp3;

for (double ii =start_time; ii < end_time; ii += dt) {

    for (auto& track: tracks) {

        if (ii !=start_time) {
            track.PropagateModel(dt);
        }

        State rand_state = Model::GetRandomState(fov);

        tmp0 = sys.source_container_.GenerateRandomMeasurement(m0.source_index,Source0::MatMeasCov::Identity()*sqrt(noise)*0, track.state_,test_data.transform_state0_, test_data.transform_data_t_m_0_);
        tmp1 = sys.source_container_.GenerateRandomMeasurement(m1.source_index,Source1::MatMeasCov::Identity()*sqrt(noise)*0, track.state_,test_data.transform_state1_, test_data.transform_data_t_m_1_);
        tmp2 = sys.source_container_.GenerateRandomMeasurement(m2.source_index,Source2::MatMeasCov::Identity()*sqrt(noise)*0, rand_state,  test_data.transform_state2_, test_data.transform_data_t_m_2_);
        tmp3 = sys.source_container_.GenerateRandomMeasurement(m3.source_index,Source0::MatMeasCov::Identity()*sqrt(noise)*0, track.state_,test_data.transform_state0_, test_data.transform_data_t_m_0_);

        m0.time_stamp = ii;
        m0.pose = tmp0.pose;
        m0.transform_state = test_data.transform_state0_;
        m0.transform_data_t_m = test_data.transform_data_t_m_0_;
        m0.transform_meas = test_data.transform_measurement0_;
        m0.transform_data_m_t = test_data.transform_data_m_t_0_;
        m1.time_stamp = ii;
        m1.pose = tmp1.pose;
        m1.twist = tmp1.twist;
        m1.transform_state = test_data.transform_state1_;
        m1.transform_data_t_m = test_data.transform_data_t_m_1_;
        m1.transform_meas = test_data.transform_measurement1_;
        m1.transform_data_m_t = test_data.transform_data_m_t_1_;
        m2.time_stamp = ii;
        m2.pose = tmp2.pose;
        m2.twist = tmp2.twist;
        m2.transform_state = test_data.transform_state2_;
        m2.transform_data_t_m = test_data.transform_data_t_m_2_;
        m2.transform_meas = test_data.transform_measurement2_;
        m2.transform_data_m_t = test_data.transform_data_m_t_2_;
        m3.time_stamp = ii;
        m3.pose = tmp3.pose;
        m3.transform_state = test_data.transform_state0_;
        m3.transform_data_t_m = test_data.transform_data_t_m_0_;
        m3.transform_meas = test_data.transform_measurement0_;
        m3.transform_data_m_t = test_data.transform_data_m_t_0_;
        sys.current_time_ = ii;

        // std::cout << "state: " << std::endl << track.state_.g_.data_ << std::endl;
        // std::cout << "m0: " << std::endl << m0.pose << std::endl;
        // std::cout << "m1: " << std::endl << m1.pose << std::endl;
        // std::cout << "m3: " << std::endl << m3.pose << std::endl;

        sys.data_tree_.AddMeasurement(sys, m0);
        sys.data_tree_.AddMeasurement(sys, m1);
        sys.data_tree_.AddMeasurement(sys, m2);
        sys.data_tree_.AddMeasurement(sys, m3);
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

using MyTypes = ::testing::Types<Test1,Test2,Test3,Test4,Test5,Test7,Test8>;
// using MyTypes = ::testing::Types<Test7,Test8>;
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
for (auto& sim_track: this->tracks) {

    bool found = false;
    for (auto& created_track: this->sys.models_) {

        if (sim_track.state_.OMinus(created_track.state_).norm() < 2) {
            found = true;
            ASSERT_LT(created_track.err_cov_.determinant(), 1);
            // std::cout << "det: " << created_track.err_cov_.determinant() << std::endl;
            // std::cout << "norm: " << created_track.err_cov_.norm() << std::endl;
            ASSERT_GE(created_track.cs_.Size(), this->sys.params_.RANSAC_score_minimum_requirement_);
            ASSERT_GE(created_track.model_likelihood_, 0.5);
        }
        
    }
    ASSERT_TRUE(found);

}

// // Make sure that the consensus set is not empty and that the model likelihood is greater than 0.5.
// for (auto& track : this->sys.models_) {
//     ASSERT_GE(track.cs_.Size(), this->sys.params_.RANSAC_score_minimum_requirement_);
//     ASSERT_GE(track.model_likelihood_, 0.5);
// }


// make sure that all of the measurement are removed
for (auto track_iter = this->sys.models_.begin(); track_iter != this->sys.models_.end(); ++track_iter) {
for (auto track_meas_outer = track_iter->cs_.consensus_set_.begin(); track_meas_outer != track_iter->cs_.consensus_set_.end(); ++track_meas_outer) {
for (auto track_meas_inner = track_meas_outer->begin(); track_meas_inner != track_meas_outer->end(); ++ track_meas_inner) {
for (auto cluster_iter = this->sys.data_tree_.data_.begin(); cluster_iter != this->sys.data_tree_.data_.end(); ++cluster_iter) {
for (auto cluster_meas_outer = cluster_iter->data_.begin(); cluster_meas_outer != cluster_iter->data_.end(); ++ cluster_meas_outer) {
for (auto cluster_meas_inner = cluster_meas_outer->begin(); cluster_meas_inner != cluster_meas_outer->end(); ++ cluster_meas_inner) {

    
    if (cluster_meas_inner->source_index == track_meas_inner->source_index && cluster_meas_inner->time_stamp == track_meas_inner->time_stamp) {
        
        ASSERT_GT(this->sys.source_container_.OMinus(cluster_meas_inner->source_index,*cluster_meas_inner, *track_meas_inner).norm(), 1e-8 );
        
    }


}}}}}}


}