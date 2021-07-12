#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <math.h> 
#include <time.h>
#include <stdlib.h>

#include "lie_groups/state.h"
#include "rransac/system.h"
#include "rransac/common/models/model_SEN_pos_vel.h"
#include "rransac/common/sources/source_SE3_cam_depth.h"
#include "rransac/common/transformations/trans_SE3_cam_depth.h"
#include "rransac/track_initialization/seed_policies/SE3_cam_depth_seed_policy.h"
#include "rransac/data_containers/cluster.h"

using namespace rransac;
using namespace lie_groups;

typedef SE3_se3 StateT;
typedef SourceSE3CamDepth<SE3_se3,MeasurementTypes::SE3_CAM_DEPTH,TransformSE3CamDepth> Source;
typedef SourceContainer<Source> SC;
typedef ModelSENPosVel<SC> Model;
typedef typename Source::Measurement Measurement;
typedef typename System<Model>::ClusterT ClusterT;
typedef SE3CamDepthSeedPolicy<Model> SeedPolicy;
typedef typename SeedPolicy::MatRot MatRot;
typedef System<Model> Sys;
typedef typename Model::Base::TransformDataType TransformDataType;
typedef typename Model::Base::Transformation Transformation;


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

//------------------------------------------------------------------------------------------------------------------

struct Test1 {

    // No angular velocity, no roll, no measurement noise, no transformation

    Model track_init_;
    double tolerance_ = 1e-9;
    bool transform_state_ = false;
    TransformDataType transform_data_;
    typename Source::MatMeasCov meas_noise_;

    Test1(){
        srand((unsigned int) time(0));

        MatRot R;

        Eigen::Matrix<double,2,1> rand_numbers = Eigen::Matrix<double,2,1>::Random();   
        R = ConstructRotation(0,rand_numbers(0)*M_PI/2,rand_numbers(1));

        // Start with zero roll and no angular velocity.
        track_init_.state_ = Model::GetRandomState(10);
        track_init_.state_.g_.R_ = R;
        track_init_.state_.u_.th_.setZero();
        transform_data_.setIdentity();
        meas_noise_.setZero();

    }

    double GetError(const Model& track, const StateT& est_state) {

        return StateT::OMinus(track.state_,est_state).norm();

    }
};

//-----------------------------------------------------------------------------------------------------------------

struct Test2 {

    // No angular velocity, no roll, no measurement noise, but transformation

    Model track_init_;
    double tolerance_ = 1e-9;
    bool transform_state_ = true;
    TransformDataType transform_data_;
    typename Source::MatMeasCov meas_noise_;


    Test2(){
        srand((unsigned int) time(0));

        MatRot R;

        Eigen::Matrix<double,2,1> rand_numbers = Eigen::Matrix<double,2,1>::Random();   
        R = ConstructRotation(0,rand_numbers(0)*M_PI/2,rand_numbers(1));

        // Start with zero roll and no angular velocity.
        track_init_.state_ = Model::GetRandomState(10);
        track_init_.state_.g_.R_ = R;
        track_init_.state_.u_.th_.setZero();
        transform_data_ = Transformation::GetRandomTransform(10);
        meas_noise_.setZero();


    }

    double GetError(const Model& track, const StateT& est_state) {

        return StateT::OMinus(track.state_,est_state).norm();

    }
};

//-----------------------------------------------------------------------------------------------------------------

struct Test3 {

    // No angular velocity, no roll, no transformation, but measurement noise.

    Model track_init_;
    double tolerance_ = 1;
    bool transform_state_ = true;
    TransformDataType transform_data_;
    typename Source::MatMeasCov meas_noise_;


    Test3(){
        srand((unsigned int) time(0));

        MatRot R;

        Eigen::Matrix<double,2,1> rand_numbers = Eigen::Matrix<double,2,1>::Random();   
        R = ConstructRotation(0,rand_numbers(0)*M_PI/2,rand_numbers(1));

        // Start with zero roll and no angular velocity.
        track_init_.state_ = Model::GetRandomState(10);
        track_init_.state_.g_.R_ = R;
        track_init_.state_.u_.th_.setZero();
        transform_data_ = Transformation::GetRandomTransform(10);
        meas_noise_.setZero();
        meas_noise_.diagonal() << 1e-1, 1e-3,1e-3,1e-3,1e-2,1e-2,1e-2;


    }

    double GetError(const Model& track, const StateT& est_state) {

        return StateT::OMinus(track.state_,est_state).norm();

    }
};

//-----------------------------------------------------------------------------------------------------------------


struct Test4 {

    // no initial roll, small angular velocities, no transformation, no measurement noise.

    Model track_init_;
    double tolerance_ = 5;
    bool transform_state_ = true;
    TransformDataType transform_data_;
    typename Source::MatMeasCov meas_noise_;


    Test4(){
        srand((unsigned int) time(0));

        MatRot R;

        Eigen::Matrix<double,2,1> rand_numbers = Eigen::Matrix<double,2,1>::Random();   
        R = ConstructRotation(0,rand_numbers(0)*M_PI/2,rand_numbers(1));

        // Start with zero roll and no angular velocity.
        track_init_.state_ = Model::GetRandomState(10);
        track_init_.state_.g_.R_ = R;
        track_init_.state_.u_.th_.setRandom();
        transform_data_ = Transformation::GetRandomTransform(10);
        meas_noise_.setZero();
        // meas_noise_.diagonal() << 1e-1, 1e-3,1e-3,1e-3,1e-2,1e-2,1e-2;


    }

    double GetError(const Model& track, const StateT& est_state) {

        return StateT::OMinus(track.state_,est_state).norm();

    }
};

//-----------------------------------------------------------------------------------------------------------------



template<typename _Test>
class SE3CamDepthSeedPolicyTest :  public testing::Test {

public:

_Test test_;
double noise_ = 1e-4;
double start_time_ = 0;
double end_time_ = 1;
double dt_ = 0.1;
double x_[10];
typename Model::VecCov x_vec_;
typename Model::State est_state_;
std::list<std::list<Measurement>> measurements_;
std::list<Measurement> meas_same_time_step_;
std::vector<typename ClusterT::IteratorPair> meas_subset_;
SourceParameters source_params0_;
Sys sys_;
Parameters params_;
int num_meas_per_timestep_ = 1;
Measurement meas_, tmp_;
typename ClusterT::IteratorPair iter_pair_;
Model current_track_, final_track_;
SeedPolicy seed_;

protected:

void SetUp() override {
    srand((unsigned int) time(0));

    // setup sources
    source_params0_.source_index_ = 0;
    source_params0_.type_ = SC::Source0::measurement_type_;
    source_params0_.meas_cov_ = SC::Source0::MatMeasCov::Identity()*noise_;

    // setup system
    params_.process_noise_covariance_ = Model::MatModelCov::Identity();
    sys_.params_ = params_;
    sys_.source_container_.AddSource(source_params0_);

    // setup measurements
    meas_.source_index = source_params0_.source_index_;
    meas_.type = source_params0_.type_;
    meas_.transform_state = test_.transform_state_;
    meas_.transform_meas  = test_.transform_state_;
    meas_.transform_data_m_t = test_.transform_data_.inverse();
    meas_.transform_data_t_m = test_.transform_data_;

    current_track_ = test_.track_init_;

    for (double ii = start_time_; ii < end_time_; ii += dt_) {

        if (ii != start_time_) {
            current_track_.PropagateModel(dt_);
        }
        // std::cout << "track pose: " << std::endl <<  current_track_.state_.g_.data_ << std::endl;
        meas_same_time_step_.clear();

        for (int jj = 0; jj < num_meas_per_timestep_; ++jj){
            tmp_ = sys_.source_container_.GenerateRandomMeasurement(meas_.source_index,test_.meas_noise_,current_track_.state_,meas_.transform_state,meas_.transform_data_t_m);

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



    std::random_shuffle(meas_subset_.begin(), meas_subset_.end());

    seed_.GenerateSeedPolicy(meas_subset_, sys_, x_, 10);


    for(int ii = 0; ii < 10; ++ii) {
        x_vec_(ii) = x_[ii];
    }

    Model::OPlus(est_state_,x_vec_);



}



};

using MyTypes = ::testing::Types<Test1,Test2,Test3,Test4>;
// using MyTypes = ::testing::Types<Test4>;
TYPED_TEST_SUITE(SE3CamDepthSeedPolicyTest, MyTypes);



TYPED_TEST(SE3CamDepthSeedPolicyTest, TestNoAngularVelocityConstrainedPitch) {

TypeParam test_type;

double error = test_type.GetError(this->current_track_,this->est_state_);




// std::cout << "est track pose: " << std::endl << this->est_state_.g_.data_ << std::endl;
// std::cout << "est track twist: " << std::endl << this->est_state_.u_.data_ << std::endl;

// std::cout << "true track pose: " << std::endl << this->current_track_.state_.g_.data_ << std::endl;
// std::cout << "true track twist: " << std::endl << this->current_track_.state_.u_.data_ << std::endl;
ASSERT_LT(error,test_type.tolerance_);
// std::cout << " est u: " << std::endl << est_track2.u_.data_ << std::endl;



// EXPECT_LT( (current_track.state_.g_.data_-est_state.g_.data_).norm(), 1e-10    ); 
// EXPECT_LT( (current_track.state_.u_.data_-est_state.u_.data_).norm(), 1e-10    ); 
// EXPECT_LT( (track.state_.g_.data_-est_track2.g_.data_).norm(), 0.5  ) << "This test can fail on rare occasions depending on the random samples. I ran it 100+ times in a row without it failing so it is very rare. So run it again";
// EXPECT_LT( (track.state_.u_.data_-est_track2.u_.data_).norm(), 0.5  ) << "This test can fail on rare occasions depending on the random samples. I ran it 100+ times in a row without it failing so it is very rare. So run it again";





}


//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




