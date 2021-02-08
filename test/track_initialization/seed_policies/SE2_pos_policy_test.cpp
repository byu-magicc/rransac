#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <math.h> 
#include <time.h>
#include <stdlib.h>

#include "lie_groups/state.h"
#include "rransac/system.h"
#include "rransac/common/models/model_SEN_pos_vel.h"
#include "rransac/common/sources/source_SEN_pos_vel.h"
#include "rransac/common/transformations/transformation_null.h"
#include "rransac/track_initialization/seed_policies/SE2_pos_seed_policy.h"
#include "rransac/data_containers/cluster.h"

using namespace rransac;
using namespace lie_groups;

TEST(SE2_pos_policy_test, SE3PoseTest) {

srand((unsigned int) time(0));

typedef SE2_se2 State;
typedef SourceSENPosVel<State> Source;
typedef ModelSENPosVel<State,TransformNULL> Model;

SE2PosSeedPolicy<Model> seed;
System<Model> sys;


double noise = 1e-4;
SourceParameters source_params1, source_params2;
source_params1.source_index_ = 0;
source_params1.type_ = MeasurementTypes::SEN_POS;
source_params1.meas_cov_ = Eigen::Matrix2d::Identity()*noise;

source_params2.source_index_ = 1;
source_params2.type_ = MeasurementTypes::SEN_POS_VEL;
source_params2.meas_cov_ = Eigen::Matrix4d::Identity()*noise;

Source source1, source2;
source1.Init(source_params1);
source2.Init(source_params2);

Parameters params;
params.process_noise_covariance_ = Eigen::Matrix<double,5,5>::Identity();
sys.params_ = params;
sys.sources_.push_back(source1);
sys.sources_.push_back(source2);

Meas<double> m1,m2;
m1.source_index = 0;
m1.type = MeasurementTypes::SEN_POS;
m2.source_index = 1;
m2.type = MeasurementTypes::SEN_POS_VEL;

Model track;
track.state_ = State::Random();
track.state_.g_.t_*=10;
track.state_.u_.data_(0) = fabs(track.state_.u_.data_(0));
while (fabs(track.state_.u_.data_(0) ) < 0.5) {
    track.state_.u_.data_(0)*=2;
}

track.state_.u_.data_(1) = 0;


double start_time = 0;
double end_time = 1;
double dt = 0.1;
Meas<double> tmp1, tmp2;
std::list<std::list<Meas<double>>> measurements1, measurements2;
std::list<Meas<double>> meas_time1,meas_time2;
std::vector<Cluster<double>::IteratorPair> meas_subset1;
std::vector<Cluster<double>::IteratorPair> meas_subset2;
Cluster<double>::IteratorPair iter_pair;

for (double ii = start_time; ii < end_time; ii += dt) {

    if (ii != start_time) {
        track.PropagateModel(dt);
    }
    meas_time1.clear();
    meas_time2.clear();

    tmp1 = sys.sources_[m1.source_index].GenerateRandomMeasurement(track.state_,Eigen::Matrix<double,2,2>::Identity()*sqrt(noise));
    // tmp2 = sys.sources_[m2.source_index].GenerateRandomMeasurement(track.state_,Eigen::Matrix<double,6,6>::Identity()*sqrt(noise));
    tmp2 = sys.sources_[m2.source_index].GenerateRandomMeasurement(track.state_,Eigen::Matrix<double,4,4>::Identity()*sqrt(noise));
    m1.pose = tmp1.pose;
    m1.time_stamp = ii;
    m2.pose = tmp2.pose;
    m2.twist = tmp2.twist;
    m2.time_stamp = ii;
    meas_time1.push_back(m1);
    meas_time2.push_back(m2);
    measurements1.push_back(meas_time1);
    measurements2.push_back(meas_time2);
    sys.current_time_ = ii;

}

for (auto outer_iter = measurements1.begin(); outer_iter != measurements1.end(); ++outer_iter) {
    for (auto inner_iter = outer_iter->begin(); inner_iter != outer_iter->end(); ++ inner_iter) {
        iter_pair.inner_it = inner_iter;
        iter_pair.outer_it = outer_iter;
        meas_subset1.push_back(iter_pair);
    }
}

for (auto outer_iter = measurements2.begin(); outer_iter != measurements2.end(); ++outer_iter) {
    for (auto inner_iter = outer_iter->begin(); inner_iter != outer_iter->end(); ++ inner_iter) {
        iter_pair.inner_it = inner_iter;
        iter_pair.outer_it = outer_iter;
        meas_subset2.push_back(iter_pair);
    }
}

std::random_shuffle(meas_subset1.begin(), meas_subset1.end());
std::random_shuffle(meas_subset2.begin(), meas_subset2.end());


double x[] = {0,0,0,0,0};
Eigen::Matrix<double,6,1> x_vec;
State est_track1, est_track2;

// std::cout << " track g: " << std::endl << track.state_.g_.data_ << std::endl;
// std::cout << " track u: " << std::endl << track.state_.u_.data_ << std::endl;

seed.GenerateSeedPolicy(meas_subset1, sys, x, 6);


for(int ii = 0; ii < 4; ++ii) {
    x_vec(ii) = x[ii];
}
x_vec(4) = 0;
x_vec(5) = x[4];

est_track1.g_.data_ = State::Algebra::Exp(x_vec.block(0,0,3,1));
est_track1.u_.data_ = x_vec.block(3,0,3,1);

// std::cout << " est g: " << std::endl << est_track1.g_.data_ << std::endl;
// std::cout << " est u: " << std::endl << est_track1.u_.data_ << std::endl;



seed.GenerateSeedPolicy(meas_subset2, sys, x, 6);

for(int ii = 0; ii < 4; ++ii) {
    x_vec(ii) = x[ii];
}
x_vec(4) = 0;
x_vec(5) = x[4];

est_track2.g_.data_ = State::Algebra::Exp(x_vec.block(0,0,3,1));
est_track2.u_.data_ = x_vec.block(3,0,3,1);

// std::cout << " est g: " << std::endl << est_track2.g_.data_ << std::endl;
// std::cout << " est u: " << std::endl << est_track2.u_.data_ << std::endl;



EXPECT_LT( (track.state_.g_.data_-est_track1.g_.data_).norm(), 1    ) << "This test can fail on rare occasions depending on the random samples. I ran it 100+ times in a row without it failing so it is very rare. So run it again";
EXPECT_LT( (track.state_.u_.data_-est_track1.u_.data_).norm(), 2    ) << "This test can fail on rare occasions depending on the random samples. I ran it 100+ times in a row without it failing so it is very rare. So run it again";
EXPECT_LT( (track.state_.g_.data_-est_track2.g_.data_).norm(), 0.5  ) << "This test can fail on rare occasions depending on the random samples. I ran it 100+ times in a row without it failing so it is very rare. So run it again";
EXPECT_LT( (track.state_.u_.data_-est_track2.u_.data_).norm(), 0.5  ) << "This test can fail on rare occasions depending on the random samples. I ran it 100+ times in a row without it failing so it is very rare. So run it again";





}



        // std::cerr << "newest meas: " << std::endl;
        // std::cerr << "time stamp: " << meas_subset_ordered[newest_index].inner_it->time_stamp << std::endl;
        // std::cerr << "pose: " << std::endl;
        // std::cerr << meas_subset_ordered[newest_index].inner_it->pose << std::endl;
        // std::cerr << "twist: " << std::endl;
        // std::cerr << meas_subset_ordered[newest_index].inner_it->twist << std::endl;

        // std::cerr << "middle meas: " << std::endl;
        // std::cerr << "time stamp: " << meas_subset_ordered[middle_index].inner_it->time_stamp << std::endl;
        // std::cerr << "pose: " << std::endl;
        // std::cerr << meas_subset_ordered[middle_index].inner_it->pose << std::endl;
        // std::cerr << "twist: " << std::endl;
        // std::cerr << meas_subset_ordered[middle_index].inner_it->twist << std::endl;

        // std::cerr << "oldest meas: " << std::endl;
        // std::cerr << "time stamp: " << meas_subset_ordered[oldest_index].inner_it->time_stamp << std::endl;
        // std::cerr << "pose: " << std::endl;
        // std::cerr << meas_subset_ordered[oldest_index].inner_it->pose << std::endl;
        // std::cerr << "twist: " << std::endl;
        // std::cerr << meas_subset_ordered[oldest_index].inner_it->twist << std::endl;