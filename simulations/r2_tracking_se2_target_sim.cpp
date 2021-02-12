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
#include "rransac/common/models/model_RN.h"
#include "rransac/common/sources/source_SEN_pos_vel.h"
#include "rransac/common/sources/source_RN.h"
#include "rransac/common/measurement/measurement_base.h"
#include "rransac/parameters.h"
#include "rransac/track_initialization/lmle_policies/linear_lmle_policy.h"
#include "rransac/track_initialization/seed_policies/null_seed_policy.h"
#include "rransac/common/transformations/trans_homography.h"
#include "rransac/common/transformations/transformation_null.h"
#include "rransac/track_initialization/ransac.h"
#include "rransac/common/data_association/model_policies/model_pdf_policy.h"
#include "rransac/common/data_association/cluster_data_tree_policies/data_tree_cluster_association_policy.h"
#include "rransac/rransac.h"
#include "rransac/common/utilities.h"
#include "rransac/visualization/visualization_host.h"
#include "rransac/visualization/draw_meas_policies/draw_meas_R2_SE2_pos_policy.h"
#include "rransac/visualization/draw_track_policies/draw_track_policy_SE2.h"
#include "rransac/visualization/draw_track_policies/draw_track_policy_R2.h"


using namespace lie_groups;
using namespace rransac;


struct Scenario1 {
    public:

    typedef ModelSENPosVel<SE2_se2, TransformNULL> TargetModel_;
    typedef typename TargetModel_::State TargetState_;
    typedef typename TargetModel_::Source TargetSource_;
    typedef typename TargetState_::Algebra TargetAlgebra_;

    typedef ModelRN<R2_r2,TransformNULL,SourceRN> TrackingModel_;
    typedef typename TrackingModel_::Transformation Transformation_;
    typedef typename TrackingModel_::Transformation::MatData TransformMatData_;
    typedef typename TrackingModel_::State TrackingState_;
    typedef typename TrackingState_::Algebra TrackingAlgebra_;
    typedef typename TrackingModel_::Source TrackingSource_;
    typedef RRANSACTemplateParameters<R2_r2,SourceRN,TransformNULL,ModelRN,NULLSeedPolicy,LinearLMLEPolicy,ModelPDFPolicy,DataTreeClusterAssociationPolicy> RRANSACParameters;
    typedef RRANSAC<RRANSACParameters> RRANSAC_;
    typedef typename RRANSAC_::tRansac RANSAC_;
    TransformMatData_ transform_data;
    // static constexpr bool transform_data_ = true;
    static constexpr bool transform_data_ = false;



    typedef Eigen::Matrix<double,2,2> MatR_;
    typedef Eigen::Matrix<double,4,4> MatR2_;
    static constexpr MeasurementTypes MeasurementType1= MeasurementTypes::RN_POS;
    static constexpr MeasurementTypes MeasurementType2= MeasurementTypes::RN_POS_VEL;
    typedef Eigen::Matrix<double,4,4> TrackingProcessNoiseCov_;
    typedef Eigen::Matrix<double,5,5> TargetProcessNoiseCov_;
    std::vector<TargetState_> target_states;
    typedef Eigen::Matrix<double,2,1> VecU_;
    std::string test_name = "R2 tracking SE2 target Test";
    Eigen::Matrix<double,TargetAlgebra_::dim_,TargetAlgebra_::dim_> target_noise_mat;
    Eigen::Matrix<double,TrackingAlgebra_::dim_,TrackingAlgebra_::dim_> track_noise_mat;
    

    Scenario1() {
        double pos = 15;
        double rot = 0.1;
        // double t_vel = 2;
        // double a_vel = 0.3;        
        double t_vel = 2;
        double a_vel = 0.6;
        double th = 0.01;
        TargetState_ state;
        Eigen::Matrix<double,3,1> pose;
        for (int ii = 0; ii < 4; ++ii)
            target_states.push_back(state);
        pose << pos, pos, 0;
        target_states[0].g_.data_ = TargetState_::Algebra::Exp(pose);
        target_states[0].u_.data_ << t_vel, 0,0;
        pose << pos, -pos, 0;
        target_states[1].g_.data_ = TargetState_::Algebra::Exp(pose);
        target_states[1].u_.data_ << t_vel,0,0;
        pose << -pos, -pos, rot;
        // pose << -pos, -pos, 0;
        target_states[2].g_.data_ = TargetState_::Algebra::Exp(pose);
        target_states[2].u_.data_ << t_vel, 0, a_vel;
        // states[2].u_.data_ << t_vel, 0, 0;
        pose << -pos, pos, -rot;
        // pose << -pos, pos, 0;
        target_states[3].g_.data_ = TargetState_::Algebra::Exp(pose);
        target_states[3].u_.data_ << t_vel,0,-a_vel;
        // states[3].u_.data_ << t_vel,0,0;
        transform_data << cos(th), -sin(th), 0, sin(th), cos(th), 0, 0, 0, 1;
        target_noise_mat.setZero();
        target_noise_mat(0,0) = 1;
        track_noise_mat.setIdentity();
    }

  
};


//-------------------------------------------------------------------------------------------------------------

template<typename T> 
class RRANSACSimulation {

public:

typedef typename T::TrackingModel_ TrackingModel_;
typedef typename T::TrackingState_  TrackingState_;
typedef typename T::TrackingSource_ TrackingSource_;
typedef typename T::TargetModel_ TargetModel_;
typedef typename T::TargetState_  TargetState_;
typedef typename T::TargetSource_ TargetSource_;
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
    source_params1.RANSAC_inlier_probability_ = 0.95;
    source_params1.gate_probability_ = 0.9;
    source_params1.spacial_density_of_false_meas_ = 0.03125;

    source_params2.type_ = T::MeasurementType2;
    source_params2.source_index_ = 1;
    source_params2.meas_cov_ = T::MatR2_::Identity()*noise_;
    source_params2.RANSAC_inlier_probability_ = 0.95;
    source_params2.gate_probability_ = 0.9;
    source_params2.spacial_density_of_false_meas_ = 0.03125;


    TrackingSource_ source1, source2;
    source1.Init(source_params1);
    source2.Init(source_params2);
    sources_.push_back(source1);
    sources_.push_back(source2);

    rransac_.AddSource(source_params1);
    rransac_.AddSource(source_params2);

    // Setup system
    Parameters params;
    params.process_noise_covariance_ = T::TrackingProcessNoiseCov_::Identity()*noise_;
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
    params.track_good_model_threshold_ = 100;
    params.track_max_missed_detection_time_ = 2;
    // params.nonlinear_innov_cov_id_ = true;

    rransac_.SetSystemParameters(params);


    // Setup Measurements
    m1_.source_index = 0;
    m1_.type = source_params1.type_ ;


    m2_.source_index = 1;
    m2_.type = source_params2.type_;


    Parameters target_params;
    target_params.process_noise_covariance_ = T::TargetProcessNoiseCov_::Identity()*noise_;

    // Setup tracks
    tracks_.resize(4);
    for (int ii = 0; ii < 4; ++ii) {
        tracks_[ii].Init(target_params);
        tracks_[ii].state_ = test_data_.target_states[ii];
    }


}


//---------------------------------------------------------------------------------------------
void Propagate(double start_time, double end_time, std::vector<int>& track_indices) {


    Meas<double> tmp1, tmp2;
    std::list<Meas<double>> new_measurements;
    Eigen::Matrix<double,1,1> rand_num;

    for (double ii =start_time; ii < end_time; ii += this->dt_) {

        new_measurements.clear();

        // Only produce measurements for the first three targets
        for (auto track_index : track_indices ) {


            // std::cerr << "track index: " << track_index << std::endl;

            auto& track = this->tracks_[track_index];


            // std::cerr << "propagate " << std::endl;
            if (ii !=this->start_time_) {
                track.state_.u_.data_ += test_data_.target_noise_mat*sqrt(this->noise_)*rransac::utilities::GaussianRandomGenerator(T::TargetAlgebra_::dim_)*this->dt_;
                track.PropagateModel(this->dt_);
            }

            // std::cerr << "transform " << std::endl;

            // if (T::transform_data_) {
            //     transformation_.SetData(test_data_.transform_data);
            //     transformation_.TransformTrack(track.state_, track.err_cov_);
            // }

            tmp_track_.state_.g_.data_ = track.state_.g_.t_;
            tmp_track_.state_.u_.data_ = track.state_.g_.R_*track.state_.u_.p_;

            // Generates measurements according to the probability of detection
            rand_num.setRandom();
            if ( ii + dt_ >= end_time) // Ensure there is a measurement at the last time step
                rand_num << 0;


            // std::cerr << "meas " << std::endl;

            if (fabs(rand_num(0,0)) < this->sources_[this->m1_.source_index].params_.probability_of_detection_) {

                tmp1 = this->sources_[this->m1_.source_index].GenerateRandomMeasurement(tmp_track_.state_,T::MatR_ ::Identity()*sqrt(this->noise_));
                tmp2 = this->sources_[this->m2_.source_index].GenerateRandomMeasurement(tmp_track_.state_,T::MatR2_::Identity()*sqrt(this->noise_));

                this->m1_.time_stamp = ii;
                this->m1_.pose = tmp1.pose;
                this->m2_.time_stamp = ii;
                this->m2_.pose = tmp2.pose;
                this->m2_.twist = tmp2.twist;

                new_measurements.push_back(this->m1_);
                new_measurements.push_back(this->m2_);
            }


        }

        // TrackingState_ rand_state;
        // for (int jj =0; jj < this->num_false_meas_; ++jj) {

        //     rand_state.g_.data_ = Eigen::Matrix<double,2,1>::Random()*this->fov_;
        //     rand_state.u_.data_ = T::VecU_::Random();
        //     tmp1 = this->sources_[this->m1_.source_index].GenerateRandomMeasurement(rand_state,T::MatR_ ::Identity()*sqrt(this->noise_));


        //     rand_state.g_.data_ = Eigen::Matrix<double,2,1>::Random()*this->fov_;
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
            this->rransac_.AddMeasurements(new_measurements,test_data_.transform_data);
        } else {
            this->rransac_.AddMeasurements(new_measurements);
        }
        
        this->rransac_.RunTrackInitialization();
        // viz_.DrawSystem(sys_);
        this->rransac_.RunTrackManagement();

        // viz_.DrawSystem(sys_);
        std::vector<TrackingModel_> tracks_to_draw;
        for (auto track_index : track_indices ) {
            tmp_track_.state_.g_.data_ = this->tracks_[track_index].state_.g_.t_;
            tmp_track_.state_.u_.data_ = this->tracks_[track_index].state_.g_.R_*this->tracks_[track_index].state_.u_.p_;
            tracks_to_draw.push_back(tmp_track_);
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


Meas<double> m1_, m2_, m3_, m4_;
double noise_ = 0.5;
T test_data_;
std::vector<TargetModel_> tracks_;
RRANSAC_ rransac_;
const System<TrackingModel_>* sys_;
std::vector<TrackingSource_> sources_;
Transformation_ transformation_;
TrackingModel_ tmp_track_;

unsigned int num_false_meas_ = 200;

// Simulation Parameters
double dt_ = 0.1;
double end_time_ = 5; // seconds;
double start_time_ = 0; // seconds;
double fov_ = 50;  // The surveillance region is a square centered at zero with side length 20
VisualizationHost<TrackingModel_, DrawMeasR2SE2PosPolicy, DrawTrackPolicyR2> viz_;



};


//-------------------------------------------------------------------------------------------------------------


int main(int argc, char *argv[]) {

    int dim1 = 1000;
    int dim2 = 1000;
    double scale = 12;
    std::vector<int> img_dimensions = {1080,1920};

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

    RRANSACSimulation<Scenario1> sim(img_dimensions, draw_info, "/home/mark/Videos/sim1.mp4", 10);

    sim.SetUp();



    std::vector<int> track_indices = {0,1,2};

    sim.Propagate(sim.start_time_,sim.end_time_,track_indices);

    track_indices = {1,2,3};

    sim.Propagate(sim.end_time_+sim.dt_,sim.end_time_*2.0,track_indices);

    sim.Propagate(sim.end_time_*2.0+sim.dt_,sim.end_time_*3.0+sim.dt_*5,track_indices);




    
    

    return 0;

}
