#include <gtest/gtest.h>
#include <stdlib.h>
#include <time.h>
#include <Eigen/Dense>
#include <string.h>
#include <chrono>
#include <vector>

#include "lie_groups/state.h"
#include "track_initialization/ransac_test.h"
#include "rransac/common/models/model_SEN_pos_vel.h"
#include "rransac/common/sources/source_SE3_cam_depth.h"
#include "rransac/common/transformations/trans_SE3_cam_depth.h"
#include "rransac/track_initialization/lmle_policies/nonlinear_lmle_policy.h"
#include "rransac/track_initialization/seed_policies/SE3_cam_depth_seed_policy.h"
#include "rransac/track_initialization/ransac.h"
#include "rransac/common/data_association/measurement_weight_policies/mw_ipdaf_policy.h"
#include "rransac/common/data_association/track_likelihood_info_policies/tli_ipdaf_policy.h"
#include "rransac/common/data_association/validation_region_policies/validation_region_innov_policy.h"

using namespace lie_groups;
using namespace rransac;

struct Test1 {
    public:
    typedef SourceSE3CamDepth<SE3_se3,MeasurementTypes::SE3_CAM_DEPTH,TransformSE3CamDepth> SourceSE3Depth;
    typedef SourceContainer<SourceSE3Depth,SourceSE3Depth,SourceSE3Depth> SC;

    typedef ModelSENPosVel<SC> Model;
    typedef typename Model::State State;
    typedef typename State::Algebra Algebra;
    typedef Ransac<Model, SE3CamDepthSeedPolicy, NonLinearLMLEPolicy, ValidationRegionInnovPolicy, TLI_IPDAFPolicy, MW_IPDAFPolicy> RANSAC;
    typedef typename Model::Base::Measurement Measurement;
    typedef typename Model::Base::Transformation Transformation;
    typedef typename Model::Base::TransformDataType TransformDataType;
    typedef Eigen::Matrix<double,3,3> MatRot;

    bool transform_state0_ = false;
    bool transform_state1_ = false;
    bool transform_state2_ = false;
    bool transform_measurement0_ = false;
    bool transform_measurement1_ = false;
    bool transform_measurement2_ = false;
    TransformDataType transform_data_t_m_0_;
    TransformDataType transform_data_t_m_1_;
    TransformDataType transform_data_t_m_2_;
    TransformDataType transform_data_m_t_0_;
    TransformDataType transform_data_m_t_1_;
    TransformDataType transform_data_m_t_2_;
    std::vector<State> states;
    std::string test_name = "SE3 Cam Test";


    Test1() {
        srand((unsigned int) time(0));
        double pos = 20;
        double rot1 = 0.2;
        double rot2 = -0.2;
        double rot3 = 0;
        double t_vel = 0.5;
        double a_vel = 0.1;
        Eigen::Matrix<double,6,1> pose;
        states.resize(4);
        pose << pos, pos, pos, rot1, rot2, rot3;
        // states[0].g_.data_ = State::Algebra::Exp(pose);
        // states[0].u_.data_ << t_vel, t_vel, t_vel, 0, 0, 0;
        // pose << -pos, pos, -pos, rot1, -rot2, rot3;
        // states[1].g_.data_ = State::Algebra::Exp(pose);
        // states[1].u_.data_ << -t_vel, t_vel, t_vel, a_vel, -a_vel, a_vel;
        // pose << pos, -pos, pos, -rot1, rot2, -rot3;
        // states[2].g_.data_ = State::Algebra::Exp(pose);
        // states[2].u_.data_ << t_vel, -t_vel, t_vel, -a_vel, a_vel, -a_vel;
        // pose << -pos, -pos, -pos, -rot1, -rot2, -rot3;
        // states[3].g_.data_ = State::Algebra::Exp(pose);
        // states[3].u_.data_ << -t_vel, -t_vel, -t_vel, 0, 0, 0;
        for (auto& state: states) {
            Eigen::Matrix<double,2,1> rand_num;
            rand_num.setRandom();
            state = Model::GetRandomState(10);
            state.g_.R_ = ConstructRotation(0,rand_num(0)*M_PI/2,rand_num(1)*M_PI);
            state.u_.th_.setRandom();
            state.u_.p_ << t_vel,0,0;
        }

        states[0].g_.t_ << pos, pos, 0;
        states[1].g_.t_ << -pos, 0, 0;
        states[2].g_.t_ << pos/2, pos, 0;
        states[3].g_.t_ << -pos*2, pos, 0;
    }

    bool SimilarStates(const State& state1, const State& state2){
        double error1, error2;
        State tmp1 = state1;
        State tmp2 = state2;
        error1 = (tmp1.g_.t_ - tmp2.g_.t_).norm();
        Model::PropagateState(tmp1,2);
        Model::PropagateState(tmp2,2);
        error2 = (tmp1.g_.t_ - tmp2.g_.t_).norm();

        if (error1 < 7 && error2 < 7) {
            return true;
        } else {
            return false;
        }

    }

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

  
};

//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------




//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------

using MyTypes = ::testing::Types<Test1>;
TYPED_TEST_SUITE(RANSACTest, MyTypes);

TYPED_TEST(RANSACTest, FullTest) {

    TypeParam test;

    auto start = std::chrono::system_clock::now();
    this->ransac.Run(this->sys);
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()/1e9;
    std::cout << "RANSAC test benchmark for " << this->test_data.test_name << " : " <<  elapsed << '\n';


// for (auto& created_track: this->sys.models_) {
//     std::cout << "created_track g: " << std::endl << created_track.state_.g_.data_ << std::endl;
//     std::cout << "created_track u: " << std::endl << created_track.state_.u_.data_ << std::endl << std::endl;

// }

// for (auto& sim_track: this->tracks) {
//     std::cout << "sim_track g: " << std::endl << sim_track.state_.g_.data_ << std::endl;
//     std::cout << "sim_track u: " << std::endl << sim_track.state_.u_.data_ << std::endl << std::endl;

// }

// make sure that the tracks were created
ASSERT_GE(this->sys.models_.size(), 4 );
for (auto& sim_track: this->tracks) {

    bool found = false;
    for (auto& created_track: this->sys.models_) {

        if (test.SimilarStates(sim_track.state_,created_track.state_)) {
            found = true;
            ASSERT_LT(created_track.err_cov_.determinant(), 1);
            // std::cout << "det: " << created_track.err_cov_.determinant() << std::endl;
            // std::cout << "norm: " << created_track.err_cov_.norm() << std::endl;
            ASSERT_GE(created_track.cs_.Size(), this->sys.params_.RANSAC_score_minimum_requirement_);
            ASSERT_GE(created_track.model_likelihood_, 0.5);
        }
        
    }
    ASSERT_TRUE(found);

}

// // Make sure that the consensus set is not empty and that the model likelihood is greater than 0.5.
// for (auto& track : this->sys.models_) {
//     ASSERT_GE(track.cs_.Size(), this->sys.params_.RANSAC_score_minimum_requirement_);
//     ASSERT_GE(track.model_likelihood_, 0.5);
// }


// make sure that all of the measurement are removed
for (auto track_iter = this->sys.models_.begin(); track_iter != this->sys.models_.end(); ++track_iter) {
for (auto track_meas_outer = track_iter->cs_.consensus_set_.begin(); track_meas_outer != track_iter->cs_.consensus_set_.end(); ++track_meas_outer) {
for (auto track_meas_inner = track_meas_outer->begin(); track_meas_inner != track_meas_outer->end(); ++ track_meas_inner) {
for (auto cluster_iter = this->sys.data_tree_.data_.begin(); cluster_iter != this->sys.data_tree_.data_.end(); ++cluster_iter) {
for (auto cluster_meas_outer = cluster_iter->data_.begin(); cluster_meas_outer != cluster_iter->data_.end(); ++ cluster_meas_outer) {
for (auto cluster_meas_inner = cluster_meas_outer->begin(); cluster_meas_inner != cluster_meas_outer->end(); ++ cluster_meas_inner) {

    
    if (cluster_meas_inner->source_index == track_meas_inner->source_index && cluster_meas_inner->time_stamp == track_meas_inner->time_stamp) {
        
        ASSERT_GT(this->sys.source_container_.OMinus(cluster_meas_inner->source_index,*cluster_meas_inner, *track_meas_inner).norm(), 1e-8 );
        
    }


}}}}}}


}