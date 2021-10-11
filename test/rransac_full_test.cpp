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
#include "rransac/common/transformations/trans_homography.h"
#include "rransac/track_initialization/ransac.h"
#include "rransac/common/data_association/validation_region_policies/validation_region_innov_policy.h"
#include "rransac/common/data_association/track_likelihood_info_policies/tli_ipdaf_policy.h"
#include "rransac/common/data_association/measurement_weight_policies/mw_ipdaf_policy.h"
#include "rransac/rransac.h"
#include "rransac/common/utilities.h"



using namespace lie_groups;
using namespace rransac;

struct Test1 {
    public:

    typedef SourceRN<R2_r2,MeasurementTypes::RN_POS,TransformNULL> SourceR2Pos;
    typedef SourceRN<R2_r2,MeasurementTypes::RN_POS_VEL,TransformNULL> SourceR2PosVEL;
    typedef SourceContainer<SourceR2Pos,SourceR2PosVEL,SourceR2PosVEL> SC;


    typedef ModelRN<SC> Model_;
    typedef typename Model_::Transformation Transformation_;
    typedef typename Model_::Transformation::TransformDataType TransformDataType;
    typedef typename Model_::State State_;
    typedef typename State_::Algebra Algebra_;
    typedef RRANSACTemplateParameters<SC,ModelRN,NULLSeedPolicy,LinearLMLEPolicy,ValidationRegionInnovPolicy, TLI_IPDAFPolicy, MW_IPDAFPolicy> RRANSACParameters;
    typedef typename RRANSACParameters::_Ransac RANSAC_;
    typedef RRANSAC<RRANSACParameters> RRANSAC_;
    std::vector<State_> states_;
    std::string test_name_ = "R2 Test";



    bool transform_tracking_frame_ = false;
    std::vector<bool> transform_state_;
    std::vector<bool> transform_meas_;
    TransformDataType transform_tracking_frame_data_;
    std::vector<TransformDataType> transform_data_t_m_;
    std::vector<TransformDataType> transform_data_m_t_;

    Test1() : transform_state_(SC::num_sources_,false), transform_meas_(SC::num_sources_,false) {
        double pos = 5;
        double vel = 0.5;
        states_.resize(4);
        states_[0].g_.data_ << pos,pos;
        states_[0].u_.data_ << 0, -vel;
        states_[1].g_.data_ << pos, -pos;
        states_[1].u_.data_ << -vel,0;
        states_[2].g_.data_ << -pos, -pos;
        states_[2].u_.data_ << 0, vel;
        states_[3].g_.data_ << -pos, pos;
        states_[3].u_.data_ << vel,0;

        transform_data_t_m_.resize(SC::num_sources_);
        transform_data_m_t_.resize(SC::num_sources_);

        
    }

  
};

//---------------------------------------------------------------------------------------------------------

struct Test2 {
    public:

    typedef SourceRN<R3_r3,MeasurementTypes::RN_POS,TransformNULL> SourceRnPos;
    typedef SourceRN<R3_r3,MeasurementTypes::RN_POS_VEL,TransformNULL> SourceRnPosVEL;
    typedef SourceContainer<SourceRnPos,SourceRnPosVEL,SourceRnPosVEL> SC;

    typedef ModelRN<SC> Model_;
    typedef typename Model_::Transformation Transformation_;
    typedef typename Model_::Transformation::TransformDataType TransformDataType;
    typedef typename Model_::State State_;
    typedef typename State_::Algebra Algebra_;
    typedef RRANSACTemplateParameters<SC,ModelRN,NULLSeedPolicy,LinearLMLEPolicy,ValidationRegionInnovPolicy, TLI_IPDAFPolicy, MW_IPDAFPolicy> RRANSACParameters;
    typedef RRANSAC<RRANSACParameters> RRANSAC_;
    typedef typename RRANSACParameters::_Ransac RANSAC_;
    std::vector<State_> states_;
    std::string test_name_ = "R3 Test";


    bool transform_tracking_frame_ = false;
    std::vector<bool> transform_state_;
    std::vector<bool> transform_meas_;
    TransformDataType transform_tracking_frame_data_;
    std::vector<TransformDataType> transform_data_t_m_;
    std::vector<TransformDataType> transform_data_m_t_;

    Test2() : transform_state_(SC::num_sources_,false), transform_meas_(SC::num_sources_,false) {
        double pos = 20;
        double vel = 0.5;
        states_.resize(4);
        states_[0].g_.data_ << pos,pos, pos;
        states_[0].u_.data_ << 0, -vel, 0;
        states_[1].g_.data_ << pos, -pos, -pos;
        states_[1].u_.data_ << -vel,0,-vel;
        states_[2].g_.data_ << -pos, -pos, -pos;
        states_[2].u_.data_ << 0, vel, vel;
        states_[3].g_.data_ << -pos, pos, pos;
        states_[3].u_.data_ << vel,0,0;

        transform_data_t_m_.resize(SC::num_sources_);
        transform_data_m_t_.resize(SC::num_sources_);
    }

  
};

//---------------------------------------------------------------------------------------------------------

struct Test3 {
    public:

    typedef SourceSENPoseTwist<SE2_se2,MeasurementTypes::SEN_POSE,TransformNULL> SourceSE2Pose;
    typedef SourceSENPoseTwist<SE2_se2,MeasurementTypes::SEN_POSE_TWIST,TransformNULL> SourceSE2PoseTwist;
    typedef SourceContainer<SourceSE2Pose,SourceSE2PoseTwist,SourceSE2PoseTwist> SC;

    typedef ModelSENPoseTwist<SC> Model_;
    typedef typename Model_::Transformation Transformation_;
    typedef typename Model_::Transformation::TransformDataType TransformDataType;
    typedef typename Model_::State State_;
    typedef typename State_::Algebra Algebra_;
    typedef RRANSACTemplateParameters<SC,ModelSENPoseTwist,NULLSeedPolicy,NonLinearLMLEPolicy,ValidationRegionInnovPolicy, TLI_IPDAFPolicy, MW_IPDAFPolicy> RRANSACParameters;
    typedef RRANSAC<RRANSACParameters> RRANSAC_;
    typedef typename RRANSACParameters::_Ransac RANSAC_;
    std::vector<State_> states_;
    std::string test_name_ = "SE2 Pose Test";

    bool transform_tracking_frame_ = false;
    std::vector<bool> transform_state_;
    std::vector<bool> transform_meas_;
    TransformDataType transform_tracking_frame_data_;
    std::vector<TransformDataType> transform_data_t_m_;
    std::vector<TransformDataType> transform_data_m_t_;

    Test3() : transform_state_(SC::num_sources_,false), transform_meas_(SC::num_sources_,false) {
        double pos = 5;
        double rot = 0.1;
        double t_vel = 0.5;
        double a_vel = 0.1;
        Eigen::Matrix<double,3,1> pose;
        states_.resize(4);
        pose << pos, pos, 0;
        states_[0].g_.data_ = State_::Algebra::Exp(pose);
        states_[0].u_.data_ << t_vel, t_vel,0;
        pose << pos, -pos, 0;
        states_[1].g_.data_ = State_::Algebra::Exp(pose);
        states_[1].u_.data_ << -t_vel,-t_vel,0;
        pose << -pos, -pos, rot;
        states_[2].g_.data_ = State_::Algebra::Exp(pose);
        states_[2].u_.data_ << t_vel, 0, a_vel;
        pose << -pos, pos, -rot;
        states_[3].g_.data_ = State_::Algebra::Exp(pose);
        states_[3].u_.data_ << t_vel,0,-a_vel;

        transform_data_t_m_.resize(SC::num_sources_);
        transform_data_m_t_.resize(SC::num_sources_);
    }

  
};

//---------------------------------------------------------------------------------------------------------

struct Test4 {
    public:

    typedef SourceSENPoseTwist<SE3_se3,MeasurementTypes::SEN_POSE,TransformNULL> SourceSE3Pose;
    typedef SourceSENPoseTwist<SE3_se3,MeasurementTypes::SEN_POSE_TWIST,TransformNULL> SourceSE3PoseTwist;
    typedef SourceContainer<SourceSE3Pose,SourceSE3PoseTwist,SourceSE3PoseTwist> SC;

    typedef ModelSENPoseTwist<SC> Model_;
    typedef typename Model_::Transformation Transformation_;
    typedef typename Model_::Transformation::TransformDataType TransformDataType;
    typedef typename Model_::State State_;
    typedef typename State_::Algebra Algebra_;
    typedef RRANSACTemplateParameters<SC,ModelSENPoseTwist,NULLSeedPolicy,NonLinearLMLEPolicy,ValidationRegionInnovPolicy, TLI_IPDAFPolicy, MW_IPDAFPolicy> RRANSACParameters;
    typedef RRANSAC<RRANSACParameters> RRANSAC_;
    typedef typename RRANSACParameters::_Ransac RANSAC_;
    std::vector<State_> states_;

    std::string test_name_ = "SE3 Pose Test";


    bool transform_tracking_frame_ = false;
    std::vector<bool> transform_state_;
    std::vector<bool> transform_meas_;
    TransformDataType transform_tracking_frame_data_;
    std::vector<TransformDataType> transform_data_t_m_;
    std::vector<TransformDataType> transform_data_m_t_;

    Test4() : transform_state_(SC::num_sources_,false), transform_meas_(SC::num_sources_,false) {
        double pos = 10;
        double rot1 = 0.2;
        double rot2 = -0.2;
        double rot3 = 0;
        double t_vel = 0.5;
        double a_vel = 0.1;
        Eigen::Matrix<double,6,1> pose;
        states_.resize(4);
        pose << pos, pos, pos, rot1, rot2, rot3;
        states_[0].g_.data_ = State_::Algebra::Exp(pose);
        states_[0].u_.data_ << t_vel, t_vel, t_vel, 0, 0, 0;
        pose << -pos, pos, -pos, rot1, -rot2, rot3;
        states_[1].g_.data_ = State_::Algebra::Exp(pose);
        states_[1].u_.data_ << -t_vel, t_vel, t_vel, a_vel, -a_vel, a_vel;
        pose << pos, -pos, pos, -rot1, rot2, -rot3;
        states_[2].g_.data_ = State_::Algebra::Exp(pose);
        states_[2].u_.data_ << t_vel, -t_vel, t_vel, -a_vel, a_vel, -a_vel;
        pose << -pos, -pos, -pos, -rot1, -rot2, -rot3;
        states_[3].g_.data_ = State_::Algebra::Exp(pose);
        states_[3].u_.data_ << -t_vel, -t_vel, -t_vel, 0, 0, 0;

        transform_data_t_m_.resize(SC::num_sources_);
        transform_data_m_t_.resize(SC::num_sources_);
    }

  
};

//---------------------------------------------------------------------------------------------------------

struct Test5 {
    public:

    typedef SourceSENPosVel<SE2_se2,MeasurementTypes::SEN_POS,TransformNULL> SourceSE2Pos;
    typedef SourceSENPosVel<SE2_se2,MeasurementTypes::SEN_POS_VEL,TransformNULL> SourceSE2PosVel;
    typedef SourceContainer<SourceSE2Pos,SourceSE2PosVel,SourceSE2PosVel> SC;


    typedef ModelSENPosVel<SC> Model_;
    typedef typename Model_::Transformation Transformation_;
    typedef typename Model_::Transformation::TransformDataType TransformDataType;
    typedef typename Model_::State State_;
    typedef typename State_::Algebra Algebra_;
    typedef RRANSACTemplateParameters<SC,ModelSENPosVel,SE2PosSeedPolicy,NonLinearLMLEPolicy,ValidationRegionInnovPolicy, TLI_IPDAFPolicy, MW_IPDAFPolicy> RRANSACParameters;
    typedef RRANSAC<RRANSACParameters> RRANSAC_;
    typedef typename RRANSACParameters::_Ransac RANSAC_;
    std::vector<State_> states_;
    std::string test_name_ = "SE2 Pos Test";
    
    bool transform_tracking_frame_ = false;
    std::vector<bool> transform_state_;
    std::vector<bool> transform_meas_;
    TransformDataType transform_tracking_frame_data_;
    std::vector<TransformDataType> transform_data_t_m_;
    std::vector<TransformDataType> transform_data_m_t_;

    Test5() : transform_state_(SC::num_sources_,false), transform_meas_(SC::num_sources_,false) {
        double pos = 5;
        double rot = 0.1;
        double t_vel = 0.5;
        double a_vel = 0.1;
        double th = 0.2;
        State_ state;
        Eigen::Matrix<double,3,1> pose;
        for (int ii = 0; ii < 4; ++ii)
            states_.push_back(state);
        pose << pos, pos, 0;
        states_[0].g_.data_ = State_::Algebra::Exp(pose);
        states_[0].u_.data_ << t_vel, 0,0;
        pose << pos, -pos, 0;
        states_[1].g_.data_ = State_::Algebra::Exp(pose);
        states_[1].u_.data_ << t_vel,0,0;
        pose << -pos, -pos, rot;
        // pose << -pos, -pos, 0;
        states_[2].g_.data_ = State_::Algebra::Exp(pose);
        states_[2].u_.data_ << t_vel, 0, a_vel;
        // states[2].u_.data_ << t_vel, 0, 0;
        pose << -pos, pos, -rot;
        // pose << -pos, pos, 0;
        states_[3].g_.data_ = State_::Algebra::Exp(pose);
        states_[3].u_.data_ << t_vel,0,-a_vel;

        transform_data_t_m_.resize(SC::num_sources_);
        transform_data_m_t_.resize(SC::num_sources_);

    }

  
};

//--------------------------------------------------------------------------------------------------------


struct Test6 {
    public:

    typedef lie_groups::State<Rn,double,3,3> StateR3_3;
    typedef SourceRN<StateR3_3,MeasurementTypes::RN_POS,TransformNULL> SourceRnPos;
    typedef SourceRN<StateR3_3,MeasurementTypes::RN_POS,TransformNULL> SourceRnPosVEL;
    typedef SourceContainer<SourceRnPos,SourceRnPosVEL,SourceRnPosVEL> SC;

    typedef ModelRN<SC> Model_;
    typedef typename Model_::Transformation Transformation_;
    typedef typename Model_::Transformation::TransformDataType TransformDataType;
    typedef typename Model_::State State_;
    typedef typename State_::Algebra Algebra_;
    typedef RRANSACTemplateParameters<SC,ModelRN,NULLSeedPolicy,LinearLMLEPolicy,ValidationRegionInnovPolicy, TLI_IPDAFPolicy, MW_IPDAFPolicy> RRANSACParameters;
    typedef RRANSAC<RRANSACParameters> RRANSAC_;
    typedef typename RRANSACParameters::_Ransac RANSAC_;
    std::vector<State_> states_;
    std::string test_name_ = "R3 Test";


    bool transform_tracking_frame_ = false;
    std::vector<bool> transform_state_;
    std::vector<bool> transform_meas_;
    TransformDataType transform_tracking_frame_data_;
    std::vector<TransformDataType> transform_data_t_m_;
    std::vector<TransformDataType> transform_data_m_t_;

    Test6() : transform_state_(SC::num_sources_,false), transform_meas_(SC::num_sources_,false) {
        double pos = 20;
        double vel = 0.5;
        double accel = 0.1*0;
        double jerk = -0.01*0;
        states_.resize(4);
        states_[0].g_.data_ << pos,pos, pos;
        states_[0].u_.data_ << 0, -vel, 0,accel,0,0,jerk,0,0;
        states_[1].g_.data_ << pos, -pos, -pos;
        states_[1].u_.data_ << -vel,0,-vel,0,accel,0,0,jerk,0;
        states_[2].g_.data_ << -pos, -pos, -pos;
        states_[2].u_.data_ << 0, vel, vel,0,0,accel,0,0,jerk;
        states_[3].g_.data_ << -pos, pos, pos;
        states_[3].u_.data_ << vel,0,0,-accel,0,0,-jerk,0,0;

        transform_data_t_m_.resize(SC::num_sources_);
        transform_data_m_t_.resize(SC::num_sources_);
    }

  
};


//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------

template<typename T> 
class RRANSACTest : public testing::Test {

public:

typedef typename T::Model_ Model_;
typedef typename T::State_ State_;
typedef typename T::RANSAC_ RANSAC_;
typedef typename T::RRANSAC_ RRANSAC_;
typedef typename T::Transformation_ Transformation_;
typedef typename T::SC::Source0 Source0_;
typedef typename T::SC::Source1 Source1_;
typedef typename T::SC::Source2 Source2_;
typedef typename T::Model_::Measurement Measurement_;
typedef typename Transformation_::TransformDataType TransformDataType_;

void SetUp() override {

    sys_ = rransac_.GetSystemInformation();
    transformation_.Init();

    // Setup sources
    std::vector<SourceParameters> source_params;
    source_params.resize(T::SC::num_sources_);
    
    source_params[0].type_ = Source0_::measurement_type_;
    source_params[0].source_index_ = 0;
    source_params[0].meas_cov_ = Source0_::MatMeasCov::Identity()*noise_;
    source_params[0].gate_probability_ = 0.9;

    source_params[1].type_ = Source1_::measurement_type_;
    source_params[1].source_index_ = 1;
    source_params[1].meas_cov_ = Source1_::MatMeasCov::Identity()*noise_;
    source_params[1].gate_probability_ = 0.9;

    source_params[2].type_ = Source2_::measurement_type_;
    source_params[2].source_index_ = 2;
    source_params[2].meas_cov_ = Source2_::MatMeasCov::Identity()*noise_;
    source_params[2].gate_probability_ = 0.9;

    // Source_ source1, source2, source3;
    // source1.Init(source_params1);
    // source2.Init(source_params2);
    // source3.Init(source_params3);
    // sources_.push_back(source1);
    // sources_.push_back(source2);
    // sources_.push_back(source3);
    for (int ii = 0; ii < T::SC::num_sources_; ++ii){
        rransac_.AddSource(source_params[ii]);
    }


    // Setup system
    Parameters params;
    params.process_noise_covariance_ = Model_::MatModelCov::Identity()*noise_*0.5;
    params.RANSAC_max_iters_ = 5;
    params.RANSAC_minimum_subset_ = 5; 
    params.RANSAC_score_stopping_criteria_ = 15;
    params.RANSAC_score_minimum_requirement_ = 13;
    params.meas_time_window_ = end_time_ - start_time_;                   // 5 seconds
    params.cluster_time_threshold_ = 5;
    params.cluster_velocity_threshold_ = 1.2;
    params.cluster_position_threshold_ = 1.2;
    params.cluster_min_size_requirement_ = 30;
    params.track_max_num_tracks_ = 5;
    params.track_similar_tracks_threshold_ = 1;
    params.track_good_model_threshold_ = 0.8;
    // params.sequential_else_parallel_fusion_ = false;
    // params.track_max_missed_detection_time_ = 0.5;
    // params.nonlinear_innov_cov_id_ = true;

    rransac_.SetSystemParameters(params);



    // Setup Measurements
    meas_.resize(T::SC::num_sources_);
    for (int ii = 0; ii < T::SC::num_sources_;++ii) {
        meas_[ii].source_index = source_params[ii].source_index_;
        meas_[ii].type = source_params[ii].type_;
        meas_[ii].transform_state = test_data_.transform_state_[ii];
        meas_[ii].transform_meas = test_data_.transform_meas_[ii];
        meas_[ii].transform_data_t_m = test_data_.transform_data_t_m_[ii];
        meas_[ii].transform_data_m_t = test_data_.transform_data_m_t_[ii];
    }


    // Setup tracks
    tracks_.resize(4);
    for (int ii = 0; ii < 4; ++ii) {
        tracks_[ii].Init(sys_->params_);
        tracks_[ii].state_ = test_data_.states_[ii];
    }


}


//---------------------------------------------------------------------------------------------
void Propagate(double start_time, double end_time, std::vector<int>& track_indices) {


    Measurement_ tmp;
    std::list<Measurement_> new_measurements;
    Eigen::Matrix<double,1,1> rand_num;
    bool transform_state = 0;
    TransformDataType_ EmptyMat;

    for (double ii =start_time; ii < end_time; ii += this->dt_) {

        new_measurements.clear();

        // Only produce measurements for the first three targets
        for (auto track_index : track_indices ) {


            // std::cerr << "track index: " << track_index << std::endl;

            auto& track = this->tracks_[track_index];


            // std::cerr << "propagate " << std::endl;
            if (ii !=this->start_time_) {
                Model_::OPlus(track.state_, rransac::utilities::GaussianRandomGenerator(Model_::cov_dim_)*std::sqrt(this->dt_*this->noise_)*0 );
                track.PropagateModel(this->dt_);
            }

            // std::cerr << "transform " << std::endl;

            if (test_data_.transform_tracking_frame_) {
                transformation_.SetData(test_data_.transform_tracking_frame_data_);
                transformation_.TransformTrack(track.state_, track.err_cov_);
            }

            for (int jj = 0; jj < T::SC::num_sources_; ++jj) {

                // Get true measurement
                rand_num.setRandom();
                if ( ii + dt_+1 >= end_time) // Ensure there is a measurement at the last time step
                    rand_num << 0;

                // Generate true measurement
                if (fabs(rand_num(0,0)) < this->sys_->source_container_.GetParams(jj).probability_of_detection_) {

                    tmp = this->sys_->source_container_.GenerateRandomMeasurement(this->meas_[jj].source_index,this->sys_->source_container_.GetParams(jj).meas_cov_.sqrt()*0,track.state_,meas_[jj].transform_state,meas_[jj].transform_data_t_m);

                    meas_[jj].time_stamp= ii;
                    meas_[jj].pose= tmp.pose;
                    meas_[jj].twist = tmp.twist;

                    new_measurements.push_back(meas_[jj]);
                }

                // Generate False measurement
                rand_num.setRandom();
                if (fabs(rand_num(0,0)) < prob_false_meas_) {
                    State_ rand_state = Model_::GetRandomState(this->fov_);
                    tmp = this->sys_->source_container_.GenerateRandomMeasurement(this->meas_[jj].source_index,this->sys_->source_container_.GetParams(jj).meas_cov_.sqrt(),rand_state,meas_[jj].transform_state,meas_[jj].transform_data_t_m);

                    meas_[jj].time_stamp= ii;
                    meas_[jj].pose= tmp.pose;
                    meas_[jj].twist = tmp.twist;

                    new_measurements.push_back(meas_[jj]);
                }

            }
            // Generates measurements according to the probability of detection



        }

        // std::cerr << "here0 " << std::endl;

        if (test_data_.transform_tracking_frame_) {
            this->rransac_.AddMeasurements(new_measurements,this->meas_[0].time_stamp,test_data_.transform_tracking_frame_data_);
        } else {
            this->rransac_.AddMeasurements(new_measurements,this->meas_[0].time_stamp);
        }
        this->rransac_.RunTrackInitialization();

        this->rransac_.RunTrackManagement();


    }
}

//---------------------------------------------------------------------------------------------


std::vector<Measurement_> meas_;
double noise_ = 1e-1;
T test_data_;
std::vector<Model_> tracks_;
RRANSAC_ rransac_;
const System<Model_>* sys_;
// std::vector<Source_> sources_;
Transformation_ transformation_;

// Simulation Parameters
double dt_ = 0.1;
double end_time_ = 5; // seconds;
double start_time_ = 0; // seconds;
double fov_ = 50;  // The surveillance region is a square centered at zero with side length 20
double prob_false_meas_ = 0.1;


};

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------

// using MyTypes = ::testing::Types<Test1,Test2,Test3,Test4,Test5, Test6>;
// using MyTypes = ::testing::Types<Test1,Test2,Test3,Test4>;
using MyTypes = ::testing::Types<Test4>;
TYPED_TEST_SUITE(RRANSACTest, MyTypes);


TYPED_TEST(RRANSACTest, FullTest) {

std::vector<int> track_indices = {0,1,2};
this->Propagate(this->start_time_,this->end_time_,track_indices);
// Get the current score of the model likelihoods of the good models
std::vector<double> model_likelihood(this->sys_->good_models_.size(),0);
for (auto& created_track: this->sys_->good_models_) {
    model_likelihood[created_track->label_] = created_track->model_likelihood_;
}

for (auto& created_track: this->sys_->models_) {
    std::cout << "created_track g: " << std::endl << created_track.state_.g_.data_ << std::endl;
    std::cout << "created_track u: " << std::endl << created_track.state_.u_.data_ << std::endl << std::endl;

}

for (auto& sim_track: this->tracks_) {
    std::cout << "sim_track g: " << std::endl << sim_track.state_.g_.data_ << std::endl;
    std::cout << "sim_track u: " << std::endl << sim_track.state_.u_.data_ << std::endl << std::endl;

}


// make sure that the tracks were created
ASSERT_GE(this->sys_->models_.size(), 3 );

for (auto index : track_indices) {
    bool found = false;
    for (auto& created_track: this->sys_->models_) {

        if (this->tracks_[index].state_.OMinus(created_track.state_).norm() < 2) {
            found = true;
            ASSERT_LT(created_track.err_cov_.determinant(), 1); // error covariance should have gotten smaller
            ASSERT_GT(created_track.model_likelihood_, this->sys_->params_.track_good_model_threshold_);
        }

    }
    ASSERT_TRUE(found);

}



// there should be three good models
ASSERT_GE(this->sys_->good_models_.size(),3);
EXPECT_EQ(this->sys_->good_models_.size(),3) << "There could be more than one good model depending on the noise";


track_indices = {1,2,3};
// std::cerr << "here5 " << std::endl;

this->Propagate(this->end_time_+this->dt_,this->end_time_*2.0,track_indices);

// std::cerr << "here6 " << std::endl;


for (auto& created_track: this->sys_->models_) {
    std::cout << "created_track g: " << std::endl << created_track.state_.g_.data_ << std::endl;
    std::cout << "created_track u: " << std::endl << created_track.state_.u_.data_ << std::endl << std::endl;

}

for (auto& sim_track: this->tracks_) {
    std::cout << "sim_track g: " << std::endl << sim_track.state_.g_.data_ << std::endl;
    std::cout << "sim_track u: " << std::endl << sim_track.state_.u_.data_ << std::endl << std::endl;

}

// make sure that the tracks were created
int num_models = this->sys_->models_.size();
ASSERT_GE(this->sys_->models_.size(), 4 );


for (auto index : track_indices) {
    bool found = false;
    for (auto& created_track: this->sys_->models_) {

        if (this->tracks_[index].state_.OMinus(created_track.state_).norm() < 2) {
            found = true;
            ASSERT_LT(created_track.err_cov_.determinant(), 1); // error covariance should have gotten smaller
            ASSERT_GT(created_track.model_likelihood_, this->sys_->params_.track_good_model_threshold_);
        }

    }
    ASSERT_TRUE(found);

}



// One of the good tracks stopped receiving measurements so it's likelihood should be less while the other two should be 
// larger. Since all of the tracks received measurements, none of them should have stayed the same. 
double gate = 1e-6;
int num_increase = 0;
int num_decrease = 0;
int num_constant = 0;
for (auto& created_track: this->sys_->models_) {
    if (created_track.label_ < model_likelihood.size() && created_track.label_ >=0) {
        if ( created_track.model_likelihood_ < model_likelihood[created_track.label_]-gate) {
            num_decrease++;
        } else if (created_track.model_likelihood_ > model_likelihood[created_track.label_]-gate) {
            num_increase++;
        } else {
            num_constant++;
        }
    }
    
}

ASSERT_GE(num_increase+num_constant,2); // They should've increased or stayed the same. They would stay the same if their likelihood is already really close to 1. 
ASSERT_GE(num_decrease,1);


// std::cerr << "here7 " << std::endl;

// Propagate the tracks some more in order to kill the track that hasn't been receiving measurements
this->Propagate(this->end_time_*2.0+this->dt_,this->end_time_*2.0+this->dt_*5,track_indices);
// std::cerr << "here8 " << std::endl;

ASSERT_LT(this->sys_->models_.size(), num_models );
// std::cerr << "here9 " << std::endl;



}