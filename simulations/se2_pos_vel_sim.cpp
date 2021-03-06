#include <gtest/gtest.h>
#include <stdlib.h>
#include <time.h>
#include <Eigen/Dense>
#include <string.h>
#include <chrono>
#include <vector>

#include "lie_groups/state.h"
#include "rransac/system.h"
#include "rransac/common/models/model_SEN_pos_vel.h"
#include "rransac/common/sources/source_SEN_pos_vel.h"
#include "rransac/common/measurement/measurement_base.h"
#include "rransac/parameters.h"
#include "rransac/track_initialization/seed_policies/SE2_pos_seed_policy.h"
#include "rransac/track_initialization/lmle_policies/nonlinear_lmle_policy.h"
#include "rransac/common/transformations/trans_homography.h"
#include "rransac/track_initialization/ransac.h"
#include "rransac/common/data_association/validation_region_policies/validation_region_innov_policy.h"
#include "rransac/common/data_association/track_likelihood_info_policies/tli_ipdaf_policy.h"
#include "rransac/common/data_association/measurement_weight_policies/mw_ipdaf_policy.h"
#include "rransac/rransac.h"
#include "rransac/common/utilities.h"
#include "rransac/visualization/visualization_host.h"
#include "rransac/visualization/draw_meas_policies/draw_meas_R2_SE2_pos_policy.h"
#include "rransac/visualization/draw_track_policies/draw_track_policy_SE2.h"


using namespace lie_groups;
using namespace rransac;


struct Scenario1 {
    public:

    typedef SourceSENPosVel<SE2_se2,MeasurementTypes::SEN_POS,TransformNULL> SourceSE2PosNull;
    typedef SourceSENPosVel<SE2_se2,MeasurementTypes::SEN_POS_VEL,TransformNULL> SourceSE2PosVelNull;
    typedef SourceContainer<SourceSE2PosNull,SourceSE2PosVelNull> SourceContainerSE2PosVelNull;

    typedef ModelSENPosVel<SourceContainerSE2PosVelNull> Model_;
    typedef typename Model_::Transformation Transformation_;
    typedef typename Model_::Transformation::MatData TransformMatData_;
    typedef typename Model_::State State_;
    typedef typename State_::Algebra Algebra_;
    typedef RRANSACTemplateParameters<SourceContainerSE2PosVelNull,ModelSENPosVel,SE2PosSeedPolicy,NonLinearLMLEPolicy,ValidationRegionInnovPolicy, TLI_IPDAFPolicy, MW_IPDAFPolicy> RRANSACParameters;
    typedef typename RRANSACParameters::TRansac RANSAC_;
    typedef RRANSAC<RRANSACParameters> RRANSAC_;
    TransformMatData_ transform_data;
    // static constexpr bool transform_data_ = true;
    static constexpr bool transform_data_ = false;



    typedef Eigen::Matrix<double,2,2> MatR_;
    typedef Eigen::Matrix<double,4,4> MatR2_;
    static constexpr MeasurementTypes MeasurementType1= MeasurementTypes::SEN_POS;
    static constexpr MeasurementTypes MeasurementType2= MeasurementTypes::SEN_POS_VEL;
    typedef Eigen::Matrix<double,5,5> ProcessNoiseCov_;
    std::vector<State_> states;
    typedef Eigen::Matrix<double,3,1> VecU_;
    std::string test_name = "SE2 Pos Test";
    Eigen::Matrix<double,Algebra_::dim_,Algebra_::dim_> noise_mat;
    

    Scenario1() {
        double pos = 15;
        double rot = 0.1;
        double t_vel = 2;
        // double a_vel = 0.3;
        double a_vel = 0.6;
        double th = 0.01;
        State_ state;
        Eigen::Matrix<double,3,1> pose;
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
        if (transform_data_){
            transform_data << cos(th), -sin(th), 0, sin(th), cos(th), 0, 0, 0, 1;
        }
        noise_mat.setZero();
        noise_mat(0,0) = 1;
    }

  
};


//-------------------------------------------------------------------------------------------------------------

template<typename T> 
class RRANSACSimulation {

public:

typedef typename T::Model_ Model_;
typedef typename T::State_ State_;
typedef typename T::RANSAC_ RANSAC_;
typedef typename T::RRANSAC_ RRANSAC_;
typedef typename T::Transformation_ Transformation_;
static constexpr bool transform_data_ = T::transform_data_;


RRANSACSimulation(const std::vector<int>& img_dimensions, const DrawInfo& draw_info ) : viz_(img_dimensions, draw_info) {}
RRANSACSimulation(const std::vector<int>& img_dimensions, const DrawInfo& draw_info, const std::string file_path ) : viz_(img_dimensions, draw_info, file_path) {}
RRANSACSimulation(const std::vector<int>& img_dimensions, const DrawInfo& draw_info, const std::string file_name, const double fps ) : viz_(img_dimensions, draw_info, file_name, fps) {}


void SetUp() {

    sys_ = rransac_.GetSystemInformation();
    transformation_.Init();

    // Setup sources
    SourceParameters source_params1, source_params2, source_params3;
    source_params1.type_ = T::MeasurementType1;
    source_params1.source_index_ = 0;
    source_params1.meas_cov_ = T::MatR_::Identity()*noise_;
    source_params1.gate_probability_ = 0.95;
    source_params1.spacial_density_of_false_meas_ = 0.001;

    source_params2.type_ = T::MeasurementType2;
    source_params2.source_index_ = 1;
    source_params2.meas_cov_ = T::MatR2_::Identity()*noise_;
    source_params2.gate_probability_ = 0.95;
    source_params2.spacial_density_of_false_meas_ = 0.001;




    rransac_.AddSource(source_params1);
    rransac_.AddSource(source_params2);

    // Setup system
    Parameters params;
    params.process_noise_covariance_ = T::ProcessNoiseCov_::Identity()*noise_;
    params.RANSAC_max_iters_ = 5;
    params.RANSAC_minimum_subset_ = 3;
    params.RANSAC_score_stopping_criteria_ = 19;
    params.RANSAC_score_minimum_requirement_ = 15;
    params.meas_time_window_ = end_time_ - start_time_;                   // 5 seconds
    // params.meas_time_window_ = 2;                   // 5 seconds
    params.cluster_time_threshold_ = 0.5;
    params.cluster_velocity_threshold_ = 2.5;
    params.cluster_position_threshold_ = 1.2;
    params.cluster_min_size_requirement_ = 5;
    params.track_max_num_tracks_ = 5;
    params.track_similar_tracks_threshold_ = 1;
    params.track_good_model_threshold_ = 0.8;
    params.track_max_missed_detection_time_ = 2;
    // params.nonlinear_innov_cov_id_ = true;

    rransac_.SetSystemParameters(params);


    // Setup Measurements
    m1_.source_index = 0;
    m1_.type = source_params1.type_ ;


    m2_.source_index = 1;
    m2_.type = source_params2.type_;



    // Setup tracks
    tracks_.resize(4);
    for (int ii = 0; ii < 4; ++ii) {
        tracks_[ii].Init(sys_->params_);
        tracks_[ii].state_ = test_data_.states[ii];
    }


}


//---------------------------------------------------------------------------------------------
void Propagate(double start_time, double end_time, std::vector<int>& track_indices) {


    Measurement tmp1, tmp2;
    std::list<Measurement> new_measurements;
    Eigen::Matrix<double,1,1> rand_num;
    bool transform_state = false;
    Eigen::MatrixXd EmptyMat;

    for (double ii =start_time; ii < end_time; ii += this->dt_) {

        new_measurements.clear();

        // Only produce measurements for the first three targets
        for (auto track_index : track_indices ) {


            // std::cerr << "track index: " << track_index << std::endl;

            auto& track = this->tracks_[track_index];


            // std::cerr << "propagate " << std::endl;
            if (ii !=this->start_time_) {
                track.state_.u_.data_ += test_data_.noise_mat*sqrt(this->noise_)*rransac::utilities::GaussianRandomGenerator(T::Algebra_::dim_)*this->dt_;
                track.PropagateModel(this->dt_);
            }

            // std::cerr << "transform " << std::endl;

            // if (T::transform_data_) {
            //     transformation_.SetData(test_data_.transform_data);
            //     transformation_.TransformTrack(track.state_, track.err_cov_);
            // }

            // Generates measurements according to the probability of detection
            rand_num.setRandom();
            if ( ii + dt_ >= end_time) // Ensure there is a measurement at the last time step
                rand_num << 0;


            // std::cerr << "meas " << std::endl;

            if (fabs(rand_num(0,0)) < this->sys_->source_container_.GetParams(m1_.source_index).probability_of_detection_) {

                tmp1 = this->sys_->source_container_.GenerateRandomMeasurement(this->m1_.source_index, T::MatR_ ::Identity()*sqrt(this->noise_),track.state_,transform_state,EmptyMat);
                tmp2 = this->sys_->source_container_.GenerateRandomMeasurement(this->m2_.source_index, T::MatR2_::Identity()*sqrt(this->noise_),track.state_,transform_state,EmptyMat);

                this->m1_.time_stamp = ii;
                this->m1_.pose = tmp1.pose;
                this->m2_.time_stamp = ii;
                this->m2_.pose = tmp2.pose;
                this->m2_.twist = tmp2.twist;

                new_measurements.push_back(this->m1_);
                new_measurements.push_back(this->m2_);
            }


        }

        // State_ rand_state;
        // for (int jj =0; jj < this->num_false_meas_; ++jj) {

        //     rand_state.g_.R_ = so2<double>::Exp(Eigen::Matrix<double,1,1>::Random()*3);
        //     rand_state.g_.t_ = Eigen::Matrix<double,2,1>::Random()*this->fov_;
        //     rand_state.u_.data_ = T::VecU_::Random();
        //     tmp1 = this->sources_[this->m1_.source_index].GenerateRandomMeasurement(rand_state,T::MatR_ ::Identity()*sqrt(this->noise_));


        //     rand_state.g_.R_ = so2<double>::Exp(Eigen::Matrix<double,1,1>::Random()*3);
        //     rand_state.g_.t_ = Eigen::Matrix<double,2,1>::Random()*this->fov_;
        //     rand_state.u_.data_ = T::VecU_::Random();
        //     tmp2 = this->sources_[this->m2_.source_index].GenerateRandomMeasurement(rand_state,T::MatR2_::Identity()*sqrt(this->noise_));


        //     this->m1_.time_stamp = ii;
        //     this->m1_.pose = tmp1.pose;
        //     this->m2_.time_stamp = ii;
        //     this->m2_.pose = tmp2.pose;
        //     this->m2_.twist = tmp2.twist;
        //     new_measurements.push_back(this->m1_);
        //     new_measurements.push_back(this->m2_);
        // }


        

        if (T::transform_data_) {
            this->rransac_.AddMeasurements(new_measurements,new_measurements.front().time_stamp,test_data_.transform_data);
        } else {
            this->rransac_.AddMeasurements(new_measurements,new_measurements.front().time_stamp);
        }
        
        this->rransac_.RunTrackInitialization();
        // viz_.DrawSystem(sys_);
        this->rransac_.RunTrackManagement();

        // viz_.DrawSystem(sys_);
        std::vector<Model_> tracks_to_draw;
        for (auto track_index : track_indices ) {
            tracks_to_draw.push_back(this->tracks_[track_index]);
        }

        viz_.DrawClusters(sys_,true);
        viz_.DrawTrueTracks(tracks_to_draw,sys_,false);
        viz_.DrawUnAssociatedMeasurements(sys_,false);
        viz_.DrawEstimatedTracks(sys_,false);
        viz_.DrawNewMeasurements(new_measurements,sys_,false);
        viz_.RecordImage();

    }
}
//---------------------------------------------------------------------------------------------


Measurement m1_, m2_, m3_, m4_;
double noise_ = 0.5;
T test_data_;
std::vector<Model_> tracks_;
RRANSAC_ rransac_;
const System<Model_>* sys_;
Transformation_ transformation_;

unsigned int num_false_meas_ = 200;

// Simulation Parameters
double dt_ = 0.1;
double end_time_ = 5; // seconds;
double start_time_ = 0; // seconds;
double fov_ = 50;  // The surveillance region is a square centered at zero with side length 20
VisualizationHost<Model_, DrawMeasR2SE2PosPolicy, DrawTrackPolicySE2> viz_;



};


//-------------------------------------------------------------------------------------------------------------


int main(int argc, char *argv[]) {

    int dim1 = 1000;
    int dim2 = 1000;
    double scale = 12;
    std::vector<int> img_dimensions = {1920,1080};

    if (argc > 1) {
        dim1 = std::atof(argv[1]);
        dim2 = std::atof(argv[2]);
        scale = std::atof(argv[3]);
        img_dimensions[0] = dim1;
        img_dimensions[1] = dim2;
    }

    DrawInfo draw_info;
    draw_info.scale_drawing = scale;
    draw_info.draw_validation_region = true;
    draw_info.draw_measurment_velocity_position_threshold = false;
    draw_info.draw_cluster_velocity_position_threshold = true;
    draw_info.draw_poor_tracks = true;
    draw_info.flip_image_x_axis = true;



    RRANSACSimulation<Scenario1> sim(img_dimensions, draw_info, "/home/mark/Videos/sim1.mp4", 10);

    sim.SetUp();

    // std::vector<int> track_indices = {0,1,2,3};
    // sim.Propagate(sim.start_time_,10,track_indices);


    std::vector<int> track_indices = {0,1,2};

    sim.Propagate(sim.start_time_,sim.end_time_,track_indices);

    track_indices = {1,2,3};
// std::cerr << "here5 " << std::endl;

    sim.Propagate(sim.end_time_+sim.dt_,sim.end_time_*2.0,track_indices);

    sim.Propagate(sim.end_time_*2.0+sim.dt_,sim.end_time_*3.0+sim.dt_*5,track_indices);




    
    

    return 0;

}
