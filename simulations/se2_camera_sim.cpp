
#include "se2_camera_sim.h"



using namespace lie_groups;
using namespace rransac;

CameraSimSE2::CameraSimSE2(const CameraData& camera_data, double dt, double end_time, int num_tracks) : camera_data_(camera_data) {

    DrawInfo draw_info;
    double scale = camera_data_.K(0,0);
    int fps = 10;
    draw_info.scale_drawing = scale;
    draw_info.draw_validation_region = true;
    draw_info.draw_measurment_velocity_position_threshold = true;
    draw_info.draw_cluster_velocity_position_threshold = true;
    draw_info.draw_poor_tracks = true;
    draw_info.flip_image_x_axis = true;
    std::string file_name = "/home/mark/Videos/sim_cam001.mp4";
    std::vector<int> img_dimensions = {camera_data_.width,camera_data_.height};
    viz_.Setup(img_dimensions,draw_info,file_name,fps);

    start_time_ =0;
    dt_ = dt;
    end_time_ = end_time;

    sys_ = rransac_.GetSystemInformation();
    noise_ = 0.05;
    transformation_.Init();

    // Setup sources
    SourceParameters source_params1, source_params2, source_params3;
    source_params1.type_ = MeasurementType;
    source_params1.source_index_ = 0;
    source_params1.meas_cov_ =camera_data_.R;
    source_params1.RANSAC_inlier_probability_ = 0.80;
    source_params1.gate_probability_ = 0.80;
    source_params1.spacial_density_of_false_meas_ = camera_data_.lambda;

    Source_ source1;
    source1.Init(source_params1, std::bind(&CameraSimSE2::InsurveillanceRegion, this, std::placeholders::_1));
    sources_.push_back(source1);
    rransac_.AddSource(source_params1, std::bind(&CameraSimSE2::InsurveillanceRegion, this, std::placeholders::_1));

    Parameters params;
    params.process_noise_covariance_ = ProcessNoiseCov_::Identity()*noise_;
    params.RANSAC_max_iters_ = 5;
    params.RANSAC_minimum_subset_ = 3;
    params.RANSAC_score_stopping_criteria_ = 10;
    params.RANSAC_score_minimum_requirement_ = 5;
    params.meas_time_window_ = 2;                   
    params.cluster_time_threshold_ = 1;
    params.cluster_velocity_threshold_ = 0.4;
    params.cluster_position_threshold_ = 0.1;
    params.cluster_min_size_requirement_ = 5;
    params.track_max_num_tracks_ = num_tracks+5;
    params.track_similar_tracks_threshold_ = 0.1;
    params.track_good_model_threshold_ = 100;
    params.track_max_missed_detection_time_ = 1;

    rransac_.SetSystemParameters(params);

    // Setup Measurements
    m_.source_index = 0;
    m_.type = source_params1.type_ ;


    // Setup tracks
    tracks_.resize(num_tracks);
    for (int ii = 0; ii < num_tracks; ++ii) {
        State_ state = GenerateRandomState();
        tracks_[ii].Init(sys_->params_);
        tracks_[ii].state_ = state;
    }

    double th = 0.01;
    t_data_ << cos(th), -sin(th), 0, sin(th), cos(th), 0, 0, 0, 1;
    noise_mat_.setZero();
    noise_mat_(0,0) = 1;

}

//----------------------------------------------------------------------

typename CameraSimSE2::State_ CameraSimSE2::GenerateRandomState() {

    State_ x;
    std::default_random_engine gen;
    std::uniform_real_distribution<double> dist_x(camera_data_.minx,camera_data_.maxx);
    std::uniform_real_distribution<double> dist_y(camera_data_.miny,camera_data_.maxy);
    std::uniform_real_distribution<double> dist_th(-M_PI, M_PI);
    std::uniform_real_distribution<double> dist_w(-0.2,0.2);
    std::uniform_real_distribution<double> dist_v(-1/camera_data_.altitude,1/camera_data_.altitude);

    Eigen::Matrix<double,3,1> pose;
    pose << 0,0,dist_th(gen);
    x.g_.data_ = State_::Algebra::Exp(pose);
    x.g_.t_ << dist_x(gen), dist_y(gen);
    x.u_.data_ << dist_v(gen),0,dist_w(gen);
    return x;

}

//---------------------------------------------------------------------------

bool CameraSimSE2::InsurveillanceRegion(const State_& state) {

    if (state.g_.t_(0) <=  camera_data_.maxx && state.g_.t_(0) >=  camera_data_.minx && state.g_.t_(1) <=  camera_data_.maxy && state.g_.t_(1) >=  camera_data_.miny) {
        return true;
    } else {
        return false;
    }
}

//-----------------------------------------------------------------------------

void CameraSimSE2::Propagate(double start_time, double end_time) {


    Meas<double> tmp1, tmp2;
    std::list<Meas<double>> new_measurements;
    Eigen::Matrix<double,1,1> rand_num;

    for (double ii =start_time; ii < end_time; ii += this->dt_) {

        new_measurements.clear();

        // Only produce measurements for the first three targets
        for (auto& track : tracks_ ) {


            // std::cerr << "propagate " << std::endl;
            if (ii !=this->start_time_) {
                // track.state_.u_.data_ += noise_mat_*sqrt(this->noise_)*rransac::utilities::GaussianRandomGenerator(Algebra_::dim_)*this->dt_;
                track.PropagateModel(this->dt_);
            }

            // std::cerr << "transform " << std::endl;

            if (transform_data_) {
                transformation_.SetData(t_data_);
                transformation_.TransformTrack(track.state_, track.err_cov_);
            }

            // Generates measurements according to the probability of detection
            rand_num.setRandom();
            if ( ii + dt_ >= end_time) // Ensure there is a measurement at the last time step
                rand_num << 0;


            // std::cerr << "meas " << std::endl;

            if (fabs(rand_num(0,0)) < this->sources_[this->m_.source_index].params_.probability_of_detection_) {

                tmp1 = this->sources_[this->m_.source_index].GenerateRandomMeasurement(track.state_,MatR_ ::Identity()*sqrt(this->noise_));

                this->m_.time_stamp = ii;
                this->m_.pose = tmp1.pose;

                new_measurements.push_back(this->m_);
            }


        }

        State_ rand_state;
        for (int jj =0; jj < this->camera_data_.num_false_meas; ++jj) {

            rand_state = GenerateRandomState();
            tmp1 = this->sources_[this->m_.source_index].GenerateRandomMeasurement(rand_state,MatR_ ::Identity()*sqrt(this->noise_));

            this->m_.time_stamp = ii;
            this->m_.pose = tmp1.pose;
            new_measurements.push_back(this->m_);
        }


        

        if (transform_data_) {
            this->rransac_.AddMeasurements(new_measurements,t_data_);
        } else {
            this->rransac_.AddMeasurements(new_measurements);
        }
        
        this->rransac_.RunTrackInitialization();
        // viz_.DrawSystem(sys_);
        this->rransac_.RunTrackManagement();

        // viz_.DrawSystem(sys_);
        // std::vector<Model_> tracks_to_draw;
        // for (auto track_index : track_indices ) {
        //     tracks_to_draw.push_back(this->tracks_[track_index]);
        // }

        viz_.DrawClusters(sys_,true);
        viz_.DrawTrueTracks(tracks_,sys_,false);
        viz_.DrawUnAssociatedMeasurements(sys_,false);
        viz_.DrawEstimatedTracks(sys_,false);
        viz_.DrawNewMeasurements(new_measurements,sys_,false);
        viz_.RecordImage();

    }
}

//----------------------------------------------------------------------------------------------------------------------------

int main() {

double noise = 0.0001;
Eigen::Matrix2d R;
R << noise, 0, 0, noise;
int width = 1920;
int height = 1080;
double fov = 70;
double altitude = 30;
int num_false_meas = 0;
double dt = 0.1;
double end_time = 10;
int num_tracks = 5;

CameraData cam_data( width,  height,  fov,  altitude, R,  num_false_meas);

CameraSimSE2 sim(cam_data, dt, end_time, num_tracks);

sim.Propagate(0, end_time);

return 0;

}

