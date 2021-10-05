#include <gtest/gtest.h>
#include <Eigen/Core>

#include "rransac/system.h"
#include "rransac/common/models/model_manager.h"
#include "rransac/common/transformations/transformation_null.h"
#include "rransac/common/transformations/trans_homography.h"
#include "rransac/common/models/model_RN.h"
#include "rransac/common/sources/source_RN.h"
#include "rransac/common/sources/source_container.h"
#include "rransac/common/sources/source_SEN_pose_twist.h"
#include "rransac/common/models/model_SEN_pose_twist.h"

namespace rransac
{
    
using namespace lie_groups;

/**
 * Test the AddModel function. When a model is added and there are more models than the maximum number of models, the model with the 
 * lowest model likelihood is removed. We will verify that it can remove a model from the begining, middle, and end
//  */ 
TEST(ModelManagerTest, AddModel ) {

typedef lie_groups::R3_r3 State;
typedef SourceRN<State,MeasurementTypes::RN_POS,TransformNULL> Source;
typedef SourceContainer<Source> SourceContainer;
typedef ModelRN<SourceContainer> Model;

System<Model> sys;
ModelManager<Model> model_manager;
Model model;

model_manager.AddModel(sys,model);

ASSERT_EQ(sys.models_.size(), 1);



}

/**
 *  Ensures that PropagateModel is called on all of the models
 */ 
TEST(ModelManagerTest, PropagateModel ) {

typedef lie_groups::R3_r3 State;
typedef SourceRN<State,MeasurementTypes::RN_POS,TransformNULL> Source;
typedef SourceContainer<Source> SourceContainer;
typedef ModelRN<SourceContainer> Model;

System<Model> sys;
ModelManager<Model> model_manager;  

Model model;
model.state_.g_.data_ << 1,1,1;
model.state_.u_.data_ << 2,3,4;

sys.params_.track_max_num_tracks_ = 2;
double dt = 0.1;


model_manager.AddModel(sys, model);
model_manager.AddModel(sys, model);

for (int ii = 0; ii < 10; ++ii) {
    model_manager.PropagateModels(sys,dt);
}


ASSERT_LE( (sys.models_.front().state_.u_.data_- model.state_.u_.data_).norm(), 1e-12 );
ASSERT_LE( (sys.models_.front().state_.g_.data_- model.state_.g_.data_ -model.state_.u_.data_).norm(), 1e-12 );
ASSERT_LE( (sys.models_.back().state_.u_.data_ - model.state_.u_.data_).norm(), 1e-12 );
ASSERT_LE( (sys.models_.back().state_.g_.data_ - model.state_.g_.data_ -model.state_.u_.data_).norm(), 1e-12 );


}

TEST(ModelManagerTest, ManageModels) {

typedef lie_groups::R3_r3 State;
typedef SourceRN<State,MeasurementTypes::RN_POS,TransformNULL> Source;
typedef SourceContainer<Source> SourceContainer;
typedef ModelRN<SourceContainer> Model;
typedef typename Model::Base::Measurement Measurement;

System<Model> sys;
ModelManager<Model> model_manager;
Model model, merge1, merge2; 

Measurement m1, m2;
m1.time_stamp = 0.5; 
m2.time_stamp = 1;

std::vector<Measurement> meas{m1,m2};

model.cs_.AddMeasurementsToConsensusSet(meas);
model.err_cov_.setIdentity();
model.label_ = -1;


sys.params_.track_max_num_tracks_ = 10;
sys.params_.track_good_model_threshold_ = 5;
sys.params_.track_similar_tracks_threshold_ = 2;
sys.params_.track_max_missed_detection_time_ = 0.9;
sys.current_time_ = 1;
 
// setup the models

for(int ii = 0; ii < sys.params_.track_max_num_tracks_+4; ++ii) {
    
    model.state_.g_.data_ << ii, 3*(ii+1), ii+2;
    model.state_.u_.data_ << 0.5*ii, ii, 2*ii;

    if (ii == 5)
        model.model_likelihood_ = 0.11;
    else
    {
        model.model_likelihood_ = ii;
    }
    
    if (ii == 7) {
        model.newest_measurement_time_stamp = 0;
    } else {
        model.newest_measurement_time_stamp = 1;
    }
    
    model_manager.AddModel(sys, model);

    
}

sys.models_.back().model_likelihood_ = 0.1;

// Make a model similar to model 8 so that it will be merged. Give this model the better label, newest measurment time stamp and likelihood
model.model_likelihood_ = 20;
model.newest_measurement_time_stamp = 2;
model.state_.g_.data_ << 8.5, 27.5, 10.5;
model.state_.u_.data_ << 4.5,8.5,16.5;
model.err_cov_ = Eigen::Matrix<double,6,6>::Identity()*2;
model.label_ = 0;
sys.model_label_++;

model_manager.AddModel(sys, model);

// std::cout << "sys model size: " << sys.models_.size() << std::endl;


// Get the models that will be merged together
merge2 = model;
auto iter = sys.models_.begin();
iter = std::next(iter,8);
merge1 = *iter;

// std::cout << "merge2: " << std::endl << merge2.state_.g_.data_ << std::endl << merge2.state_.u_.data_ << std::endl << merge2.err_cov_ << std::endl;
// std::cout << "merge1: " << std::endl << merge1.state_.g_.data_ << std::endl << merge1.state_.u_.data_ << std::endl << merge1.err_cov_ << std::endl;

model_manager.ManageModels(sys,0.7);

// Get the merged model
iter = sys.models_.begin();
iter = std::next(iter,5);

// std::cout << "iter: " << std::endl << iter->state_.g_.data_ << std::endl << iter->state_.u_.data_ << std::endl << iter->err_cov_ << std::endl;


// std::cout << "sys model size: " << sys.models_.size() << std::endl;


// Check size of consensus set. One measurement from each model should have been removed.
// This is to ensure that the function gets called and not to check the correctness of the function
ASSERT_EQ(sys.models_.front().cs_.consensus_set_.size(), 1);
ASSERT_EQ(sys.models_.back().cs_.consensus_set_.size(), 1);


//////////////////////////
// Verify that models were fused correctly
/////////////////////////

// The fused state should be between the two original states
ASSERT_LT(iter->state_.g_.data_(0), merge2.state_.g_.data_(0));
ASSERT_LT(iter->state_.g_.data_(1), merge2.state_.g_.data_(1));
ASSERT_LT(iter->state_.g_.data_(2), merge2.state_.g_.data_(2));
ASSERT_LT(iter->state_.u_.data_(0), merge2.state_.u_.data_(0));
ASSERT_LT(iter->state_.u_.data_(1), merge2.state_.u_.data_(1));
ASSERT_LT(iter->state_.u_.data_(2), merge2.state_.u_.data_(2));

ASSERT_GT(iter->state_.g_.data_(0), merge1.state_.g_.data_(0));
ASSERT_GT(iter->state_.g_.data_(1), merge1.state_.g_.data_(1));
ASSERT_GT(iter->state_.g_.data_(2), merge1.state_.g_.data_(2));
ASSERT_GT(iter->state_.u_.data_(0), merge1.state_.u_.data_(0));
ASSERT_GT(iter->state_.u_.data_(1), merge1.state_.u_.data_(1));
ASSERT_GT(iter->state_.u_.data_(2), merge1.state_.u_.data_(2));

// The fused error covariance should be smaller than the other two covariances
ASSERT_LT(iter->err_cov_(0,0), merge2.err_cov_(0,0));
ASSERT_LT(iter->err_cov_(1,1), merge2.err_cov_(1,1));
ASSERT_LT(iter->err_cov_(2,2), merge2.err_cov_(2,2));
ASSERT_LT(iter->err_cov_(3,3), merge2.err_cov_(3,3));
ASSERT_LT(iter->err_cov_(4,4), merge2.err_cov_(4,4));
ASSERT_LT(iter->err_cov_(5,5), merge2.err_cov_(5,5));

ASSERT_LT(iter->err_cov_(0,0), merge1.err_cov_(0,0));
ASSERT_LT(iter->err_cov_(1,1), merge1.err_cov_(1,1));
ASSERT_LT(iter->err_cov_(2,2), merge1.err_cov_(2,2));
ASSERT_LT(iter->err_cov_(3,3), merge1.err_cov_(3,3));
ASSERT_LT(iter->err_cov_(4,4), merge1.err_cov_(4,4));
ASSERT_LT(iter->err_cov_(5,5), merge1.err_cov_(5,5));

ASSERT_GT(iter->err_cov_(0,0), 0);
ASSERT_GT(iter->err_cov_(1,1), 0);
ASSERT_GT(iter->err_cov_(2,2), 0);
ASSERT_GT(iter->err_cov_(3,3), 0);
ASSERT_GT(iter->err_cov_(4,4), 0);
ASSERT_GT(iter->err_cov_(5,5), 0);

ASSERT_EQ(iter->newest_measurement_time_stamp,merge2.newest_measurement_time_stamp);
ASSERT_EQ(iter->label_,merge2.label_);
ASSERT_EQ(iter->model_likelihood_, merge2.model_likelihood_);


// Verify that models were pruned correctly and labels were assigned properly
ASSERT_EQ(sys.models_.size(), sys.params_.track_max_num_tracks_);
iter = sys.models_.begin();

ASSERT_EQ(iter->model_likelihood_, 1);
++iter;
ASSERT_EQ(iter->model_likelihood_, 2);
++iter;
ASSERT_EQ(iter->model_likelihood_, 3);
++iter;
ASSERT_EQ(iter->model_likelihood_, 4);
++iter;
ASSERT_EQ(iter->model_likelihood_, 6);
ASSERT_EQ(iter->label_,1);
++iter;
ASSERT_EQ(iter->model_likelihood_, 20);
ASSERT_EQ(iter->label_,0);
++iter;
ASSERT_EQ(iter->model_likelihood_, 9);
ASSERT_EQ(iter->label_,2);
++iter;
ASSERT_EQ(iter->model_likelihood_, 10);
ASSERT_EQ(iter->label_,3);
++iter;
ASSERT_EQ(iter->model_likelihood_, 11);
ASSERT_EQ(iter->label_,4);
++iter;
ASSERT_EQ(iter->model_likelihood_, 12);
ASSERT_EQ(iter->label_,5);

// Verify the good models vector
ASSERT_EQ(sys.good_models_[0]->model_likelihood_, 6);
ASSERT_EQ(sys.good_models_[1]->model_likelihood_, 20);
ASSERT_EQ(sys.good_models_[2]->model_likelihood_, 9);
ASSERT_EQ(sys.good_models_[3]->model_likelihood_, 10);
ASSERT_EQ(sys.good_models_[4]->model_likelihood_, 11);
ASSERT_EQ(sys.good_models_[5]->model_likelihood_, 12);


}






/**
 *  Ensures that UpdateModel is called on all of the models
 */ 
TEST(ModelManagerTest, UpdateModel ) {

typedef lie_groups::R3_r3 State;
typedef SourceRN<State,MeasurementTypes::RN_POS,TransformNULL> Source;
typedef SourceContainer<Source> SourceContainer;
typedef ModelRN<SourceContainer> Model;
typedef typename Model::Base::Measurement Measurement;

System<Model> sys;
ModelManager<Model> model_manager;  

Source source;
SourceParameters source_params;
source_params.meas_cov_ = Eigen::Matrix3d::Identity();
source_params.spacial_density_of_false_meas_ = 0.1;
source_params.type_ = MeasurementTypes::RN_POS;
source_params.probability_of_detection_ = 0.9;
source_params.gate_probability_ = 0.8;
source_params.source_index_ = 0;


sys.source_container_.AddSource(source_params);

sys.params_.track_max_num_tracks_ = 2;
sys.params_.process_noise_covariance_ = Eigen::Matrix<double,6,6>::Identity();

Measurement m;
m.source_index = 0;
m.weight = 1;
m.probability = 1;
m.time_stamp = 0.1;
m.type = MeasurementTypes::RN_POS;
m.pose = Eigen::Matrix<double,3,1>::Random();


Model model;
model.Init(sys.params_);
model.state_.g_.data_ << 1,1,1;
model.state_.u_.data_ << 2,3,4;

model.new_assoc_meas_ = std::vector<std::vector<Measurement>>{std::vector<Measurement>{m}};

double dt = 0.1;

model_manager.AddModel(sys, model);
model_manager.AddModel(sys, model);

for (int ii = 0; ii < 10; ++ii) {
    model_manager.PropagateModels(sys,dt);
}
ASSERT_DOUBLE_EQ(sys.models_.begin()->newest_measurement_time_stamp,0);

model_manager.UpdateModels(sys);


ASSERT_EQ(sys.models_.front().newest_measurement_time_stamp, 0.1);
ASSERT_EQ(sys.models_.back().newest_measurement_time_stamp, 0.1);


}

// Test the transform model

TEST(ModelManagerTest, TransformModelTest ) {

typedef lie_groups::R2_r2 State;
typedef SourceRN<State,MeasurementTypes::RN_POS,TransformHomography> Source;
typedef SourceContainer<Source> SourceContainer;
typedef ModelRN<SourceContainer> Model;
typedef typename Model::Base::Measurement Measurement;

// TransformHomography<State> trans;

// The homography data is set so that the state and the measurements will all be transformed to zero
Eigen::Matrix3d homography;
homography.setZero();
homography.block(2,2,1,1)<< 1;
// trans.Init();
// trans.SetData(homography);

System<Model> sys;
ModelManager<Model> model_manager;  

Source source;
SourceParameters source_params;
source_params.meas_cov_ = Eigen::Matrix2d::Identity();
source_params.spacial_density_of_false_meas_ = 0.1;
source_params.type_ = MeasurementTypes::RN_POS;
source_params.probability_of_detection_ = 0.9;
source_params.gate_probability_ = 0.8;
source_params.source_index_ = 0;


sys.source_container_.AddSource(source_params);
sys.params_.track_max_num_tracks_ = 2;
sys.params_.process_noise_covariance_ = Eigen::Matrix<double,4,4>::Identity();
sys.params_.transform_consensus_set_ = true;
sys.transformaion_.Init();
sys.transformaion_.SetData(homography);

Measurement m;
m.source_index = 0;
m.weight = 1;
m.probability = 1;
m.time_stamp = 0;
m.type = MeasurementTypes::RN_POS;
m.pose = Eigen::Matrix<double,2,1>::Random();


Model model;
model.Init(sys.params_);
model.state_.g_.data_ << 1,1;
model.state_.u_.data_ << 2,3;
model.err_cov_ == Eigen::Matrix4d::Identity();

model.new_assoc_meas_ = std::vector<std::vector<Measurement>>{std::vector<Measurement>{m}};

for (int ii = 0; ii < 1000; ++ii) {
    m.time_stamp = ii;
    model.cs_.AddMeasToConsensusSet(m);
}

double dt = 0.1;

model_manager.AddModel(sys, model);
model_manager.AddModel(sys, model);

model_manager.TransformModels(sys);

for (auto model_iter = sys.models_.begin(); model_iter != sys.models_.end(); ++model_iter) {

    // After the transformation it should be zero
    ASSERT_DOUBLE_EQ(model_iter->state_.g_.data_.norm(),0);

    for(auto outer_iter = model_iter->cs_.consensus_set_.begin(); outer_iter != model_iter->cs_.consensus_set_.end(); ++ outer_iter) {
        for(auto inner_iter = outer_iter->begin(); inner_iter != outer_iter->end(); ++inner_iter) {
            ASSERT_DOUBLE_EQ(inner_iter->pose.norm(), 0);
        }
    }

}

}

//-------------------------------------------------------------------------------------

TEST(ModelManagerTest, Track2TrackAssociationFusion) {

typedef lie_groups::SE3_se3 State;
typedef SourceSENPoseTwist<State, MeasurementTypes::SEN_POSE_TWIST, TransformNULL > Source;
typedef SourceContainer<Source> SC;
typedef ModelSENPoseTwist<SC> Model;

System<Model> sys;
ModelManager<Model> model_manager;

sys.params_.track_max_num_tracks_ = 10;
sys.params_.track_similar_tracks_threshold_ = 200;
sys.current_time_ = 1;

State true_state = State::Random(10);
double dt = 0.1;

typename Model::MatModelCov P1 = Model::MatModelCov::Identity()*0.05;
typename Model::MatModelCov P2 = Model::MatModelCov::Identity()*0.01;
typename Model::MatModelCov F = Model::GetLinTransFuncMatState(true_state,dt);
P1 = F*(F*(F*P1*F.transpose() + P1)*F.transpose() + P1)*F.transpose();
P2 = F*(F*(F*P2*F.transpose() + P2)*F.transpose() + P2)*F.transpose();

State est_state1 = true_state;
est_state1.OPlusEQ( P1.sqrt()* rransac::utilities::GaussianRandomGenerator(Model::cov_dim_).cwiseAbs() );
State est_state2 = true_state;
est_state2.OPlusEQ( -P2.sqrt()* rransac::utilities::GaussianRandomGenerator(Model::cov_dim_).cwiseAbs() );

// std::cout << "true state: " << std::endl << true_state.g_.data_ << std::endl << true_state.u_.data_ << std::endl;
// std::cout << "est_state1: " << std::endl << est_state1.g_.data_ << std::endl << est_state1.u_.data_ << std::endl << P1 << std::endl;
// std::cout << "est_state2: " << std::endl << est_state2.g_.data_ << std::endl << est_state2.u_.data_ << std::endl << P2 << std::endl;

Model model1, model2;
model1.state_ = est_state1;
model1.err_cov_ = P1;
model1.newest_measurement_time_stamp = 1;
model1.model_likelihood_ = 0.7;
model2.state_ = est_state2;
model2.err_cov_ = P2;
model2.newest_measurement_time_stamp = 0.5;
model2.model_likelihood_ = 0.9;
model_manager.AddModel(sys,model1);
model_manager.AddModel(sys,model2);
model_manager.ManageModels(sys,10);

ASSERT_EQ(sys.models_.size(),1);
auto iter = sys.models_.begin();
// std::cout << "merged state: " << std::endl << iter->state_.g_.data_ << std::endl << iter->state_.u_.data_ << std::endl << iter->err_cov_ << std::endl;

double d1 = State::OMinus(true_state,est_state1).norm();
double d2 = State::OMinus(true_state,est_state2).norm();
double d3 = State::OMinus(true_state,iter->state_).norm();
ASSERT_LT(d3,d1);
EXPECT_LT(d3,d2) << "This might not always be true. Run again.";
ASSERT_LT(iter->err_cov_.determinant(), P1.determinant());
ASSERT_LT(iter->err_cov_.determinant(), P2.determinant());

// std::cout << d1 << std::endl << d2 << std::endl << d3;

}


} // namespace rransac



