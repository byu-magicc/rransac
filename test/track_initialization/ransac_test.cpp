#include <gtest/gtest.h>
#include <stdlib.h>
#include <time.h>
#include <Eigen/Dense>

#include "system.h"
#include "state.h"
#include "data_containers/cluster.h"
#include "track_initialization/ransac.h"
#include "common/models/model_RN.h"
#include "common/transformations/transformation_null.h"
#include "common/sources/source_RN.h"
#include "common/data_association/model_policies/model_pdf_policy.h"
#include "track_initialization/lmle_policies/linear_lmle_policy.h"

using namespace rransac;
using namespace lie_groups;

template<typename tModel, template<typename > typename tSeed> 
struct LMLEDummy : tSeed<tModel>{
    static void GenerateHypotheticalStateEstimatePolicy(const std::vector<Cluster::IteratorPair>& meas_subset, const System<tModel>& sys);
};

template<typename tModel>
struct StatisticDummy{
    static void CalculateMeasurmentAndLikelihoodDataPolicy(const System<tModel>& sys, tModel& model);
};

template<typename tModel>
struct SeedDummy {
    typedef tModel Model;
};

TEST(RANSAC_TEST, GenerateMinimumSubsetTest) {

typedef ModelRN<R2_r2, TransformNULL> Model;




Cluster cluster;
Ransac<Model, SeedDummy, LMLEDummy, StatisticDummy> ransac;

srand(time(NULL));

Meas m;
m.pose = Eigen::Matrix<double,2,1>::Random();
m.type = MeasurementTypes::RN_POS;
unsigned int max_times = 20;
unsigned int max_meas_per_time = 100;

unsigned int num_times = 0;



while (num_times <= 1) {
    num_times = abs(rand()) % max_times;
}



for (unsigned int ii = 0; ii < num_times; ++ii) {
    m.time_stamp = ii;
    unsigned int num_meas = 0;
    while (num_meas < 1) {
        num_meas = rand() % max_meas_per_time;
    }

    for (int jj = 0; jj < num_meas; ++jj) {
        m.pose << ii, jj;
        cluster.AddMeasurement(m);
    }
}


unsigned int num_min_subset = 0;


while (num_min_subset < 1) {
    num_min_subset = abs(rand()) % num_times;
}


std::vector<Cluster::IteratorPair> meas_index = ransac.GenerateMinimumSubset(num_min_subset, cluster);

// Make sure there is a measurement from the current time step and that they are all from different time steps
int times[meas_index.size()];
bool meas_from_curr_time_step_found = false;
int index = 0;
for (auto iter_pair_iter = meas_index.begin(); iter_pair_iter != meas_index.end(); ++iter_pair_iter) {

    if(iter_pair_iter->inner_it->time_stamp == num_times -1) {
        meas_from_curr_time_step_found = true;
    }
    times[index] = iter_pair_iter->inner_it->time_stamp;
    // std::cout << "time " << index << ": " << iter_pair_iter->inner_it->time_stamp << std::endl;
    // std::cout << "meas " << index << ": " << iter_pair_iter->inner_it->pose(1,0) << std::endl;
    index++;
    
}

ASSERT_TRUE(meas_from_curr_time_step_found);

// Make sure that they are all from different time steps
for (unsigned int ii = 0; ii < meas_index.size(); ++ii) {
    for (unsigned int jj = ii+1; jj < meas_index.size(); ++jj) {
        ASSERT_FALSE(times[ii]==times[jj]);
    }
}



}


// ----------------------------------------------------------------------------------------------------------------------

TEST(RANSAC_TEST, ScoreHypotheticalStateEstimateTest) {

// typedef ModelRN<R2_r2, TransformNULL> Model;
typedef ModelRN<R2_r2, TransformNULL> Model;

// Setup sources
double noise = 1e-2;
SourceParameters source_params1, source_params2, source_params3;
source_params1.type_ = MeasurementTypes::RN_POS;
source_params1.source_index_ = 0;
source_params1.meas_cov_ = Eigen::Matrix2d::Identity()*noise;
source_params1.RANSAC_inlier_probability_ = 0.9;

source_params2.type_ = MeasurementTypes::RN_POS_VEL;
source_params2.source_index_ = 1;
source_params2.meas_cov_ = Eigen::Matrix4d::Identity()*noise;
source_params2.RANSAC_inlier_probability_ = 0.9;

source_params3.type_ = MeasurementTypes::RN_POS_VEL;
source_params3.source_index_ = 2;
source_params3.meas_cov_ = Eigen::Matrix4d::Identity()*noise;
source_params3.RANSAC_inlier_probability_ = 0.9;


SourceR2 source1,source2, source3;
source1.Init(source_params1);
source2.Init(source_params2);
source3.Init(source_params3);

// Setup system
Parameters params;
params.process_noise_covariance_ = Eigen::Matrix4d::Identity()*noise;
System<Model> sys;
sys.sources_.push_back(source1);
sys.sources_.push_back(source2);
sys.sources_.push_back(source3);
sys.params_ = params;

// Setup cluster 
Cluster cluster;

// Setup the model
Model x;
x.Init(params);
x.state_.g_.data_.setRandom();
x.state_.u_.data_.setRandom();

// Setup Measurements
Meas m1, m2, m3, m4;
m1.source_index = 0;
m1.type = MeasurementTypes::RN_POS;


m2.source_index = 1;
m2.type = MeasurementTypes::RN_POS_VEL;

// This measurement is noise
m3.source_index = 2;
m3.type = MeasurementTypes::RN_POS;

// Another measurement of source 2
m4.source_index = 1;
m4.type = MeasurementTypes::RN_POS_VEL;



// Propagate model and add in two true measurements and one false measurement per time step

int steps = 10;
double dt = 0.1;
double start_time = 0;

// std::cerr << "here 1" << std::endl;

for (double ii = start_time; ii < steps*dt; ii += dt ) {
    x.PropagateModel(dt);
    Meas tmp1 = sys.sources_[m1.source_index].GenerateRandomMeasurement(x.state_,Eigen::Matrix2d::Identity()*sqrt(noise));
    Meas tmp2 = sys.sources_[m2.source_index].GenerateRandomMeasurement(x.state_,Eigen::Matrix<double,4,4>::Identity()*sqrt(noise));
    Meas tmp4 = sys.sources_[m2.source_index].GenerateRandomMeasurement(x.state_,Eigen::Matrix<double,4,4>::Identity()*sqrt(noise));
    m1.time_stamp = ii + dt;
    m1.pose = tmp1.pose;

    m2.time_stamp = ii +dt;
    m2.pose = tmp2.pose;
    m2.twist = tmp2.twist;

    m3.time_stamp = ii+dt;
    m3.pose = x.state_.g_.data_ + Eigen::Matrix<double,2,1>::Random()*20;
    m3.twist = x.state_.u_.data_ + Eigen::Matrix<double,2,1>::Random()*20;

    m4.time_stamp = ii+dt;
    m4.pose = tmp4.pose;
    m4.twist = tmp4.twist;

    sys.current_time_ = ii+dt;

    cluster.AddMeasurement(m1);
    cluster.AddMeasurement(m2);
    cluster.AddMeasurement(m3);
    cluster.AddMeasurement(m4);
}

// std::cerr << "here 2" << std::endl;


// Get score and inliers
std::vector<Cluster::IteratorPair> inliers;
Ransac<Model, SeedDummy, LMLEDummy, ModelPDFPolicy> ransac;
int score = ransac.ScoreHypotheticalStateEstimate(x.state_,cluster,sys,inliers);

// std::cout << "score: " << score << std::endl;

ASSERT_LE(score, 2*(steps+1)+1);
ASSERT_GT(score, steps+1) << "There is a very slight chance that the test can fail here. Just run it again.";

// Make sure the inliers are in chronological order
for (auto iter = inliers.begin(); iter != std::prev(inliers.end()); ++iter) {
    ASSERT_TRUE( iter->inner_it->time_stamp <= std::next(iter)->inner_it->time_stamp);
}

Model track = ransac.GenerateTrack(x.state_,sys,inliers);

// std::cerr << " x.g: " << std::endl << x.state_.g_.data_ << std::endl;
// std::cerr << " x.u: " << std::endl << x.state_.u_.data_ << std::endl;
// std::cerr << " track.g: " << std::endl << track.state_.g_.data_ << std::endl;
// std::cerr << " track.u: " << std::endl << track.state_.u_.data_ << std::endl;
// std::cerr << "err cov: " << std::endl << track.err_cov_ << std::endl;
// std::cerr << "likelihood: " << std::endl << track.model_likelihood_ << std::endl;

// On average it is about 100
ASSERT_GT(track.model_likelihood_ , 0);

ASSERT_LT( (x.state_.g_.data_-track.state_.g_.data_).norm(), 1e-1  );
ASSERT_LT( (x.state_.u_.data_-track.state_.u_.data_).norm(), 1e-1  );

// The error covariance of x is only the initialized error covariance while track's error covariance has been updated
// with measurements and should be smaller than the inital error covariance. 
ASSERT_LT(track.err_cov_.determinant(), x.err_cov_.determinant());


}

// ---------------------------------------------------------------------------------------------

TEST(RANSAC_TEST, RUN_TEST) {

// typedef ModelRN<R2_r2, TransformNULL> Model;
typedef ModelRN<R2_r2, TransformNULL> Model;

// Setup sources
double noise = 1e-2;
SourceParameters source_params1, source_params2, source_params3;
source_params1.type_ = MeasurementTypes::RN_POS;
source_params1.source_index_ = 0;
source_params1.meas_cov_ = Eigen::Matrix2d::Identity()*noise;
source_params1.RANSAC_inlier_probability_ = 0.9;

source_params2.type_ = MeasurementTypes::RN_POS_VEL;
source_params2.source_index_ = 1;
source_params2.meas_cov_ = Eigen::Matrix4d::Identity()*noise;
source_params2.RANSAC_inlier_probability_ = 0.9;

source_params3.type_ = MeasurementTypes::RN_POS_VEL;
source_params3.source_index_ = 2;
source_params3.meas_cov_ = Eigen::Matrix4d::Identity()*noise;
source_params3.RANSAC_inlier_probability_ = 0.9;


SourceR2 source1,source2, source3;
source1.Init(source_params1);
source2.Init(source_params2);
source3.Init(source_params3);

// Setup system
Parameters params;
params.process_noise_covariance_ = Eigen::Matrix4d::Identity()*noise;
params.RANSAC_max_iters_ = 100;
params.RANSAC_minimum_subset_ = 3;
params.RANSAC_score_stopping_criteria_ = 10;
params.RANSAC_score_minimum_requirement_ = 6;
params.meas_time_window_ = 5;                   // 5 seconds
params.cluster_time_threshold_ = 2;
params.cluster_velocity_threshold_ = 1;
params.cluster_position_threshold_ = 0.5;
params.max_num_models_ = 5;

System<Model> sys;
sys.sources_.push_back(source1);
sys.sources_.push_back(source2);
sys.sources_.push_back(source3);
sys.params_ = params;

// Setup Measurements
Meas m1, m2, m3, m4;
m1.source_index = 0;
m1.type = MeasurementTypes::RN_POS;


m2.source_index = 1;
m2.type = MeasurementTypes::RN_POS_VEL;

// This measurement is noise
m3.source_index = 2;
m3.type = MeasurementTypes::RN_POS_VEL;

m4.source_index = 0;
m4.type = MeasurementTypes::RN_POS;

// setup ransac
Ransac<Model, SeedDummy, LinearLMLEPolicy, ModelPDFPolicy> ransac;

// Setup the models that will produce the measurements in the clusters
std::vector<Model> tracks(4);
tracks[0].Init(sys.params_);
tracks[1].Init(sys.params_);
tracks[2].Init(sys.params_);
tracks[3].Init(sys.params_);

double pos = 5;
double vel = 0.1;

tracks[0].state_.g_.data_ << pos,pos;
tracks[0].state_.u_.data_ << 0, -vel;
tracks[1].state_.g_.data_ << pos, -pos;
tracks[1].state_.u_.data_ << -vel,0;
tracks[2].state_.g_.data_ << -pos, -pos;
tracks[2].state_.u_.data_ << 0, vel;
tracks[3].state_.g_.data_ << -pos, pos;
tracks[3].state_.u_.data_ << vel,0;

// Create simulation data
double dt = 0.1;
double end_time = 5; // seconds;
double start_time = 0; // seconds;
double fov = 10;  // The surveillance region is a square centered at zero with side length 20
Meas tmp1, tmp2, tmp3, tmp4;

for (double ii =start_time; ii < end_time; ii += dt) {

    for (auto& track: tracks) {

        if (ii !=start_time) {
            track.PropagateModel(dt);
        }

        tmp1 = sys.sources_[m1.source_index].GenerateRandomMeasurement(track.state_,Eigen::Matrix2d::Identity()*sqrt(noise));
        tmp2 = sys.sources_[m2.source_index].GenerateRandomMeasurement(track.state_,Eigen::Matrix<double,4,4>::Identity()*sqrt(noise));
        tmp3.pose = Eigen::Matrix<double,2,1>::Random()*fov;
        tmp3.twist = Eigen::Matrix<double,2,1>::Random();
        tmp4 = sys.sources_[m4.source_index].GenerateRandomMeasurement(track.state_,Eigen::Matrix2d::Identity()*sqrt(noise));

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

        sys.data_tree_.AddMeasurement(sys, m1);
        sys.data_tree_.AddMeasurement(sys, m2);
        sys.data_tree_.AddMeasurement(sys, m3);
        sys.data_tree_.AddMeasurement(sys, m4);
    }
}

sys.data_tree_.ConstructClusters(sys);

ransac.Run(sys);

// make sure that the tracks were created
for (auto& created_track: sys.models_) {

    bool found = false;
    for (auto& sim_track: tracks) {

        if (sim_track.state_.OMinus(created_track.state_).norm() < 1e-1) {
            found = true;
            ASSERT_LT(created_track.err_cov_.norm(), 1);
        }
        
    }
    ASSERT_TRUE(found);

}

// Make sure that the consensus set is not empty and that the model likelihood is greater than 0.
for (auto& track : sys.models_) {
    ASSERT_GE(track.cs_.Size(), sys.params_.RANSAC_score_minimum_requirement_);
    ASSERT_GE(track.model_likelihood_, 10);
}


// make sure that all of the measurement are removed
for (auto track_iter = sys.models_.begin(); track_iter != sys.models_.end(); ++track_iter) {
for (auto track_meas_outer = track_iter->cs_.consensus_set_.begin(); track_meas_outer != track_iter->cs_.consensus_set_.end(); ++track_meas_outer) {
for (auto track_meas_inner = track_meas_outer->begin(); track_meas_inner != track_meas_outer->end(); ++ track_meas_inner) {
for (auto cluster_iter = sys.data_tree_.data_.begin(); cluster_iter != sys.data_tree_.data_.end(); ++cluster_iter) {
for (auto cluster_meas_outer = cluster_iter->data_.begin(); cluster_meas_outer != cluster_iter->data_.end(); ++ cluster_meas_outer) {
for (auto cluster_meas_inner = cluster_meas_outer->begin(); cluster_meas_inner != cluster_meas_outer->end(); ++ cluster_meas_inner) {

    
    if (cluster_meas_inner->source_index == track_meas_inner->source_index && cluster_meas_inner->time_stamp == track_meas_inner->time_stamp) {
        
        ASSERT_GT(sys.sources_[cluster_meas_inner->source_index].OMinus(*cluster_meas_inner, *track_meas_inner).norm(), 1e-8 );
        
    }


}}}}}}

}