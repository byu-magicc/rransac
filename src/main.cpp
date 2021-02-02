
#include <stdlib.h>
#include <time.h>
#include <Eigen/Dense>
#include <string.h>
#include <chrono>
#include <vector>

#include "system.h"
#include "state.h"
// #include "common/models/model_RN.h"
#include "common/models/model_SEN_pos_vel.h"
// #include "common/models/model_SEN_pose_twist.h"
#include "common/sources/source_RN.h"
#include "common/sources/source_SEN_pos_vel.h"
#include "common/sources/source_SEN_pose_twist.h"
#include "common/measurement/measurement_base.h"
#include "parameters.h"
#include "track_initialization/seed_policies/null_policy.h"
#include "track_initialization/seed_policies/SE2_pos_policy.h"
#include "track_initialization/lmle_policies/linear_lmle_policy.h"
#include "track_initialization/lmle_policies/nonlinear_lmle_policy.h"
#include "common/transformations/transformation_null.h"
#include "common/transformations/trans_homography.h"
#include "track_initialization/ransac.h"
#include "common/data_association/model_policies/model_pdf_policy.h"
#include "common/data_association/cluster_data_tree_policies/data_tree_cluster_association_policy.h"
#include "rransac.h"
#include "common/utilities.h"



using namespace lie_groups;
using namespace rransac;



//---------------------------------------------------------------------------------------------------------

struct Test5 {
    public:

    typedef ModelSENPosVel<SE2_se2, TransformHomography> Model_;
    typedef typename Model_::Transformation Transformation_;
    typedef typename Model_::Transformation::MatData TransformMatData_;
    typedef typename Model_::State State_;
    typedef typename State_::Algebra Algebra_;
    typedef typename Model_::Source Source_;
    typedef Ransac<Model_, SE2PosSeedPolicy, NonLinearLMLEPolicy, ModelPDFPolicy> RANSAC_;
    typedef RRANSACTemplateParameters<Model_,ModelPDFPolicy,DataTreeClusterAssociationPolicy,SE2PosSeedPolicy,NonLinearLMLEPolicy> RRANSACParameters;
    typedef RRANSAC<RRANSACParameters> RRANSAC_;
    TransformMatData_ transform_data;
    static constexpr bool transform_data_ = true;
    // static constexpr bool transform_data_ = false;



    typedef Eigen::Matrix<double,2,2> MatR_;
    typedef Eigen::Matrix<double,4,4> MatR2_;
    static constexpr MeasurementTypes MeasurementType1= MeasurementTypes::SEN_POS;
    static constexpr MeasurementTypes MeasurementType2= MeasurementTypes::SEN_POS_VEL;
    typedef Eigen::Matrix<double,5,5> ProcessNoiseCov_;
    std::vector<State_> states;
    typedef Eigen::Matrix<double,3,1> VecU_;
    std::string test_name = "SE2 Pos Test";
    Eigen::Matrix<double,Algebra_::dim_,Algebra_::dim_> noise_mat;
    

    Test5() {
        double pos = 5;
        double rot = 0.1;
        double t_vel = 0.5;
        double a_vel = 0.1;
        double th = 0.2;
        Eigen::Matrix<double,3,1> pose;
        State_ state;
        for (int ii = 0; ii < 4; ++ii)
            states.push_back(state);
        pose << pos, pos, 0;
        states[0].g_.data_ = State_::Algebra::Exp(pose);
        states[0].u_.data_ << t_vel, 0,0;
        pose << pos, -pos, 0;
        states[1].g_.data_ = State_::Algebra::Exp(pose);
        states[1].u_.data_ << t_vel,0,0;
        pose << -pos, -pos, rot;
        // pose << -pos, -pos, 0;
        states[2].g_.data_ = State_::Algebra::Exp(pose);
        states[2].u_.data_ << t_vel, 0, a_vel;
        // states[2].u_.data_ << t_vel, 0, 0;
        pose << -pos, pos, -rot;
        // pose << -pos, pos, 0;
        states[3].g_.data_ = State_::Algebra::Exp(pose);
        states[3].u_.data_ << t_vel,0,-a_vel;
        // states[3].u_.data_ << t_vel,0,0;
        transform_data << cos(th), -sin(th), 0, sin(th), cos(th), 0, 0, 0, 1;
        noise_mat.setZero();
        noise_mat(0,0) = 1;
    }

  
};


//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------

template<typename T> 
class RRANSACTest  {

public:

typedef typename T::Model_ Model_;
typedef typename T::State_ State_;
typedef typename T::Source_ Source_;
typedef typename T::RANSAC_ RANSAC_;
typedef typename T::RRANSAC_ RRANSAC_;
typedef typename T::Transformation_ Transformation_;
static constexpr bool transform_data_ = T::transform_data_;

void SetUp() {

    sys_ = rransac_.GetSystemInformation();
    transformation_.Init();

    // Setup sources
    SourceParameters source_params1, source_params2, source_params3;
    source_params1.type_ = T::MeasurementType1;
    source_params1.source_index_ = 0;
    source_params1.meas_cov_ = T::MatR_::Identity()*noise_;
    source_params1.RANSAC_inlier_probability_ = 0.95;
    source_params1.gate_probability_ = 0.9;

    source_params2.type_ = T::MeasurementType2;
    source_params2.source_index_ = 1;
    source_params2.meas_cov_ = T::MatR2_::Identity()*noise_;
    source_params2.RANSAC_inlier_probability_ = 0.95;
    source_params2.gate_probability_ = 0.9;

    source_params3.type_ = T::MeasurementType2;
    source_params3.source_index_ = 2;
    source_params3.meas_cov_ = T::MatR2_::Identity()*noise_;
    source_params3.RANSAC_inlier_probability_ = 0.95;
    source_params3.gate_probability_ = 0.9;

    Source_ source1, source2, source3;
    source1.Init(source_params1);
    source2.Init(source_params2);
    source3.Init(source_params3);
    sources_.push_back(source1);
    sources_.push_back(source2);
    sources_.push_back(source3);

    rransac_.AddSource(source_params1);
    rransac_.AddSource(source_params2);
    rransac_.AddSource(source_params3);

    // Setup system
    Parameters params;
    params.process_noise_covariance_ = T::ProcessNoiseCov_::Identity()*noise_;
    params.RANSAC_max_iters_ = 5;
    params.RANSAC_minimum_subset_ = 3;
    params.RANSAC_score_stopping_criteria_ = 15;
    params.RANSAC_score_minimum_requirement_ = 13;
    params.meas_time_window_ = end_time_ - start_time_;                   // 5 seconds
    params.cluster_time_threshold_ = 0.5;
    params.cluster_velocity_threshold_ = 1.2;
    params.cluster_position_threshold_ = 1.2;
    params.max_num_models_ = 5;
    params.similar_tracks_threshold_ = 1;
    params.good_model_threshold_ = 90;
    // params.NonLinearInnovCovId_ = true;

    rransac_.SetSystemParameters(params);


    // Setup Measurements
    m1_.source_index = 0;
    m1_.type = source_params1.type_ ;


    m2_.source_index = 1;
    m2_.type = source_params2.type_;

    // This measurement is noise
    m3_.source_index = 2;
    m3_.type = source_params3.type_;

    m4_.source_index = 0;
    m4_.type =source_params1.type_ ;

    // Setup tracks
    tracks_.resize(4);
    for (int ii = 0; ii < 4; ++ii) {
        tracks_[ii].Init(sys_->params_);
        tracks_[ii].state_ = test_data_.states[ii];
    }


}


//---------------------------------------------------------------------------------------------
void Propagate(double start_time, double end_time, std::vector<int>& track_indices) {


    Meas<double> tmp1, tmp2, tmp3, tmp4;
    std::list<Meas<double>> new_measurements;
    Eigen::Matrix<double,1,1> rand_num;

    for (double ii =start_time; ii < end_time; ii += this->dt_) {

        new_measurements.clear();

        // Only produce measurements for the first three targets
        for (auto track_index : track_indices ) {


            std::cerr << "track index: " << track_index << std::endl;

            auto& track = this->tracks_[track_index];


            // std::cerr << "propagate " << std::endl;
            if (ii !=this->start_time_) {
                track.state_.u_.data_ += test_data_.noise_mat*sqrt(this->noise_)*utilities::GaussianRandomGenerator(T::Algebra_::dim_)*this->dt_;
                track.PropagateModel(this->dt_);
            }

            std::cerr << "transform " << std::endl;

            if (T::transform_data_) {
            std::cerr << "data " << std::endl << test_data_.transform_data << std::endl;

                transformation_.SetData(test_data_.transform_data);
                transformation_.TransformTrack(track.state_, track.err_cov_);
            }

            std::cerr << "random " << std::endl;

            // Generates measurements according to the probability of detection
            rand_num.setRandom();
            if ( ii + dt_ >= end_time) // Ensure there is a measurement at the last time step
                rand_num << 0;


            std::cerr << "meas " << std::endl;

            if (fabs(rand_num(0,0)) < this->sources_[this->m1_.source_index].params_.probability_of_detection_) {

                tmp1 = this->sources_[this->m1_.source_index].GenerateRandomMeasurement(track.state_,T::MatR_ ::Identity()*sqrt(this->noise_*0.5));
                tmp2 = this->sources_[this->m2_.source_index].GenerateRandomMeasurement(track.state_,T::MatR2_::Identity()*sqrt(this->noise_*0.5));
                tmp4 = this->sources_[this->m4_.source_index].GenerateRandomMeasurement(track.state_,T::MatR_ ::Identity()*sqrt(this->noise_*0.5));

                this->m1_.time_stamp = ii;
                this->m1_.pose = tmp1.pose;
                this->m2_.time_stamp = ii;
                this->m2_.pose = tmp2.pose;
                this->m2_.twist = tmp2.twist;
                this->m4_.time_stamp = ii;
                this->m4_.pose = tmp4.pose;

                new_measurements.push_back(this->m1_);
                new_measurements.push_back(this->m2_);
                new_measurements.push_back(this->m4_);
            }

            std::cerr << "tmp3 " << std::endl;

            State_ rand_state;
            rand_state.g_.data_ = T::Algebra_::Exp(Eigen::Matrix<double,State_::Group::dim_,1>::Random()*this->fov_);
            rand_state.u_.data_ = T::VecU_::Random();
            tmp3 = this->sources_[this->m3_.source_index].GenerateRandomMeasurement(rand_state,T::MatR2_::Identity()*sqrt(this->noise_*0.5));
            this->m3_.time_stamp = ii;
            this->m3_.pose = tmp3.pose;
            this->m3_.twist = tmp3.twist;
            new_measurements.push_back(this->m3_);

        }

        std::cerr << "here0 " << std::endl;

        if (T::transform_data_) {
            this->rransac_.AddMeasurements(new_measurements,test_data_.transform_data);
        } else {
            this->rransac_.AddMeasurements(new_measurements);
        }
        std::cerr << "here1 " << std::endl;
        this->rransac_.RunTrackInitialization();
        std::cerr << "here2 " << std::endl;

        this->rransac_.RunTrackManagement();
        std::cerr << "here3 " << std::endl;


    }
}

//---------------------------------------------------------------------------------------------


Meas<double> m1_, m2_, m3_, m4_;
double noise_ = 1e-2;
T test_data_;
std::vector<Model_> tracks_;
RRANSAC_ rransac_;
const System<Model_>* sys_;
std::vector<Source_> sources_;
Transformation_ transformation_;

// Simulation Parameters
double dt_ = 0.1;
double end_time_ = 5; // seconds;
double start_time_ = 0; // seconds;
double fov_ = 50;  // The surveillance region is a square centered at zero with side length 20



};

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------

int main() {


RRANSACTest<Test5> test;
test.SetUp();
std::vector<int> track_indices = {0,1,2};
std::cerr << "here4 " << std::endl;

test.Propagate(test.start_time_,test.end_time_,track_indices);

track_indices = {1,2,3};
std::cerr << "here5 " << std::endl;

test.Propagate(test.end_time_+test.dt_,test.end_time_*2.0,track_indices);
std::cerr << "here6 " << std::endl;


test.Propagate(test.end_time_*2.0+test.dt_,test.end_time_*2.0+test.dt_*5,track_indices);

std::cerr << "here7 " << std::endl;


for (auto& created_track: test.sys_->models_) {
    std::cout << "created_track g: " << std::endl << created_track.state_.g_.data_ << std::endl;
    std::cout << "created_track u: " << std::endl << created_track.state_.u_.data_ << std::endl << std::endl;

}

for (auto& sim_track: test.tracks_) {
    std::cout << "sim_track g: " << std::endl << sim_track.state_.g_.data_ << std::endl;
    std::cout << "sim_track u: " << std::endl << sim_track.state_.u_.data_ << std::endl << std::endl;

}



    return 0;
}