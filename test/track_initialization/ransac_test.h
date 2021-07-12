#ifndef RRANSAC_TEST_TRACK_INITIALIZATION_RANSAC_TEST_
#define RRANSAC_TEST_TRACK_INITIALIZATION_RANSAC_TEST_
#pragma once


#include <gtest/gtest.h>
#include <stdlib.h>
#include <time.h>
#include <Eigen/Dense>
#include <string.h>
#include <chrono>
#include <vector>

#include "rransac/parameters.h"
#include "rransac/common/sources/source_base.h"
#include "rransac/system.h"


namespace rransac
{
    


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
typedef typename Model::Base::Measurement Measurmeent;




void SetUp() override {
srand((unsigned int) time(0));
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
params.RANSAC_minimum_subset_ = 7;
params.RANSAC_score_stopping_criteria_ = 15;
params.RANSAC_score_minimum_requirement_ = 6;
params.meas_time_window_ = 5;                   // 5 seconds
params.cluster_time_threshold_ = 0.5;
params.cluster_velocity_threshold_ = 1.5;
params.cluster_position_threshold_ = 0.6;
params.track_max_num_tracks_ = 5;
params.nonlinear_innov_cov_id_ = true;

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
        m0.twist = tmp0.twist;
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
        m3.twist = tmp3.twist;
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
// std::cout << "cluster label: " << sys.cluster_label_ << std::endl;

// for (auto& track: tracks) {
//     std::cout << "track pose: " << std::endl << track.state_.g_.data_ << std::endl;
//     std::cout << "track twist: " << std::endl << track.state_.u_.data_ << std::endl << std::endl;
// }


}

double noise = 1e-2;
System<Model> sys;
T test_data;
std::vector<Model> tracks;
RANSAC ransac;



};


} // namespace rransac

#endif //RRANSAC_TEST_TRACK_INITIALIZATION_RANSAC_TEST_