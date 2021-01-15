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

using namespace rransac;
using namespace lie_groups;

template<typename tModel>
struct LMLEDummy {
    static void GenerateStateEstimatePolicy(const std::vector<Cluster::IteratorPair>& meas_subset, const System<tModel>& sys);
};

template<typename tModel>
struct StatisticDummy{
    static void CalculateMeasurmentAndLikelihoodDataPolicy(const System<tModel>& sys, tModel& model);
};

TEST(RANSAC_TEST, GenerateMinimumSubsetTest) {

typedef ModelRN<R2_r2, TransformNULL> Model;




Cluster cluster;
Ransac<Model, LMLEDummy, StatisticDummy> ransac;

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
Ransac<Model, LMLEDummy, ModelPDFPolicy> ransac;
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
m3.type = MeasurementTypes::RN_POS;


// Setup the models that will produce the measurements in the clusters
std::vector<tracks>;
Model track1, track2, track3, track4;
track1.Init(sys.params_);
track2.Init(sys.params_);
track3.Init(sys.params_);
track4.Init(sys.params_);

double pos = 5;
double vel = -0.1;

track1.state_.g_.data_ << pos,pos;
track1.state_.u_.data_ << 0, -vel;
track2.state_.g_.data_ << pos, -pos;
track2.state_.u_.data_ << -vel,0;
track3.state_.g_.data_ << -pos, -pos;
track3.state_.u_.data_ << 0, vel;
track4.state_.g_.data_ << -pos, pos;
track4.state_.u_.data_ << vel,0;

// Create simulation data
double dt = 0.1;
double time = 10; // seconds;
Meas temp1;

for (double ii =0; ii < time; ii += dt) {

    Meas tmp1 = sys.sources_[m1.source_index].GenerateRandomMeasurement(x.state_,Eigen::Matrix2d::Identity()*sqrt(noise));
    Meas tmp2 = sys.sources_[m2.source_index].GenerateRandomMeasurement(x.state_,Eigen::Matrix<double,4,4>::Identity()*sqrt(noise));
    Meas tmp4 = sys.sources_[m2.source_index].GenerateRandomMeasurement(x.state_,Eigen::Matrix<double,4,4>::Identity()*sqrt(noise));



}


}