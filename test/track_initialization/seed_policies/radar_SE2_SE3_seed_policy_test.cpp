#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <math.h> 
#include <time.h>
#include <stdlib.h>

#include "lie_groups/state.h"
#include "rransac/system.h"
#include "rransac/common/measurement/measurement_base.h"
#include "rransac/common/models/model_SEN_pos_vel.h"
#include "rransac/common/sources/source_SE2_SE3_radar.h"
#include "rransac/common/sources/source_container.h"
#include "rransac/common/transformations/trans_radar_SE2_SE3_with_SE2_SE3.h"
#include "rransac/track_initialization/seed_policies/radar_SE2_SE3_seed_policy.h"


using namespace rransac;
using namespace lie_groups;


template<typename _State>
class RadarSE2SE3SeedPolicyTest : public testing::Test {

typedef _State State;
typedef SourceRadarSE2SE3<State,MeasurementTypes::SE2_SE3_RADAR,TransRadarSE2SE3WithSE2SE3> Source;
typedef SourceContainer<Source> SC;
typedef ModelSENPosVel<SC> Model;
typedef RadarSE2SE3SeedPolicy<Model> SeedPolicy;
typedef System<Model> Sys;
typedef typename Model::Measurement Measurement;
typedef typename System<Model>::ClusterT ClusterT;
typedef typename SeedPolicy::MatRot MatRot;
typedef typename Model::Base::TransformDataType TransformDataType;
typedef typename Model::Base::Transformation Transformation;

static constexpr bool is_SE3_ = Model::cov_dim_ ==10 ? true : false;


public:

Sys sys_;
double noise_ = 1e-4;
double start_time_ = 0;
double end_time_ = 1;
double tolerance_ = 1*State::Group::dim_pos_;
double dt_ = 0.1;
double x_[Model::cov_dim_];
typename Model::VecCov x_vec_;
typename Model::State est_state_;
std::list<std::list<Measurement>> measurements_;
std::list<Measurement> meas_same_time_step_;
std::vector<typename ClusterT::IteratorPair> meas_subset_;
SourceParameters source_params_;
Parameters params_;
int num_meas_per_timestep_ = 1;
Measurement meas_, tmp_;
typename ClusterT::IteratorPair iter_pair_;
Model current_track_, final_track_;
SeedPolicy seed_;


MatRot ConstructRotation(double roll, double pitch, double yaw) {

    MatRot R;
    double c_th = cos(pitch);
    double c_phi = cos(roll);
    double c_psi = cos(yaw);
    double s_th = sin(pitch);
    double s_phi = sin(roll);
    double s_psi = sin(yaw);

    R << c_th*c_psi, c_th*s_psi, -s_th, s_phi*s_th*c_psi-c_phi*s_psi, s_phi*s_th*s_psi+c_phi*c_psi, s_phi*c_th, c_phi*s_th*c_psi+s_phi*s_psi, c_phi*s_th*s_psi-s_phi*c_psi, c_phi*c_th;
    return R.transpose();
}

protected:

void SetUp() override {

    srand((unsigned int) time(0));
    est_state_ = State::Identity();

    // setup sources
    source_params_.source_index_ = 0;
    source_params_.type_ = SC::Source0::measurement_type_;
    source_params_.meas_cov_ = SC::Source0::MatMeasCov::Identity()*noise_;

    // setup system
    params_.process_noise_covariance_ = Model::MatModelCov::Identity();
    sys_.params_ = params_;
    sys_.source_container_.AddSource(source_params_);

    // setup measurements
    meas_.source_index = source_params_.source_index_;
    meas_.type = source_params_.type_;

    double should_transform_meas = rand();
    bool transform_meas = should_transform_meas > 0.5;
    if (transform_meas) {
        meas_.transform_state = true;
        meas_.transform_meas  = true;
        meas_.transform_data_m_t = Transformation::GetRandomTransform(10);
        meas_.transform_data_t_m = meas_.transform_data_m_t.inverse();
    }


    // double should_add_meas_noise = rand();
    // bool add_meas_noise = should_add_meas_noise > 0.5;
    // double meas_noise = 0;
    // if(add_meas_noise) {
    //     meas_noise = 1;
    //     tolerance_ = 1e-1;
    // }

    current_track_.Init(params_);
    current_track_.state_ = Model::GetRandomState(10);
    if(is_SE3_) {
        Eigen::Matrix<double,2,1> rand_numbers = Eigen::Matrix<double,2,1>::Random();
        current_track_.state_.g_.R_ = ConstructRotation(0,rand_numbers(0)*M_PI/2,rand_numbers(1));
    }
    current_track_.state_.u_.th_*=0.1;

    for (double ii = start_time_; ii < end_time_; ii += dt_) {

        if (ii != start_time_) {
            current_track_.PropagateModel(dt_);
        }
        // std::cout << "track pose: " << std::endl <<  current_track_.state_.g_.data_ << std::endl;
        meas_same_time_step_.clear();

        for (int jj = 0; jj < num_meas_per_timestep_; ++jj){
            tmp_ = sys_.source_container_.GenerateRandomMeasurement(meas_.source_index,source_params_.meas_cov_*0,current_track_.state_,meas_.transform_state,meas_.transform_data_t_m);

            meas_.pose = tmp_.pose;
            meas_.twist = tmp_.twist;
            meas_.time_stamp = ii;
            meas_same_time_step_.push_back(meas_);
            // std::cout << "meas pose: " << std::endl <<  meas_.pose<< std::endl;

        }

        measurements_.push_back(meas_same_time_step_);
        sys_.current_time_ = ii;

    }

    for (auto outer_iter = measurements_.begin(); outer_iter != measurements_.end(); ++outer_iter) {
        for (auto inner_iter = outer_iter->begin(); inner_iter != outer_iter->end(); ++ inner_iter) {
            iter_pair_.inner_it = inner_iter;
            iter_pair_.outer_it = outer_iter;
            meas_subset_.push_back(iter_pair_);
        }
    }

    seed_.GenerateSeedPolicy(meas_subset_, sys_, x_, Model::cov_dim_);

    for(int ii = 0; ii <  Model::cov_dim_; ++ii) {
        x_vec_(ii) = x_[ii];
    }

    Model::OPlus(est_state_,x_vec_);


}
}; // RadarSE2SE3SeedPolicyTest

//-----------------------------------------------------------------------------
using MyTypes = ::testing::Types<SE2_se2,SE3_se3>;
TYPED_TEST_SUITE(RadarSE2SE3SeedPolicyTest, MyTypes);



TYPED_TEST(RadarSE2SE3SeedPolicyTest, TestGenerateSeedPolicy) {


double error = TypeParam::State::OMinus(this->current_track_.state_,this->est_state_).norm();




// std::cout << "est track pose: " << std::endl << this->est_state_.g_.data_ << std::endl;
// std::cout << "est track twist: " << std::endl << this->est_state_.u_.data_ << std::endl;

// std::cout << "true track pose: " << std::endl << this->current_track_.state_.g_.data_ << std::endl;
// std::cout << "true track twist: " << std::endl << this->current_track_.state_.u_.data_ << std::endl;
ASSERT_LT(error,this->tolerance_);
// std::cout << " est u: " << std::endl << est_track2.u_.data_ << std::endl;



// EXPECT_LT( (current_track.state_.g_.data_-est_state.g_.data_).norm(), 1e-10    ); 
// EXPECT_LT( (current_track.state_.u_.data_-est_state.u_.data_).norm(), 1e-10    ); 
// EXPECT_LT( (track.state_.g_.data_-est_track2.g_.data_).norm(), 0.5  ) << "This test can fail on rare occasions depending on the random samples. I ran it 100+ times in a row without it failing so it is very rare. So run it again";
// EXPECT_LT( (track.state_.u_.data_-est_track2.u_.data_).norm(), 0.5  ) << "This test can fail on rare occasions depending on the random samples. I ran it 100+ times in a row without it failing so it is very rare. So run it again";





}




