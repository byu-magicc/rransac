#include <gtest/gtest.h>
#include <Eigen/Core>

#include "system.h"
#include "common/models/model_manager.h"
#include "common/transformations/transformation_null.h"
#include "common/transformations/trans_homography.h"
#include "common/models/model_RN.h"
#include "common/sources/source_RN.h"

namespace rransac
{
    
using namespace lie_groups;

/**
 * Test the AddModel function. When a model is added and there are more models than the maximum number of models, the model with the 
 * lowest model likelihood is removed. We will verify that it can remove a model from the begining, middle, and end
 */ 
TEST(ModelManagerTest, AddModel ) {

typedef lie_groups::R3_r3 State;
typedef ModelRN<State, TransformNULL,SourceRN> Model;

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
typedef ModelRN<State, TransformNULL,SourceRN> Model;

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

ASSERT_DOUBLE_EQ(sys.models_.begin()->missed_detection_time_, 1.0);
ASSERT_LE( (sys.models_.front().state_.u_.data_- model.state_.u_.data_).norm(), 1e-12 );
ASSERT_LE( (sys.models_.front().state_.g_.data_- model.state_.g_.data_ -model.state_.u_.data_).norm(), 1e-12 );
ASSERT_LE( (sys.models_.back().state_.u_.data_ - model.state_.u_.data_).norm(), 1e-12 );
ASSERT_LE( (sys.models_.back().state_.g_.data_ - model.state_.g_.data_ -model.state_.u_.data_).norm(), 1e-12 );


}

TEST(ModelManagerTest, ManageModels) {

typedef lie_groups::R3_r3 State;
typedef ModelRN<State, TransformNULL,SourceRN> Model;
typedef Meas<double> Measurement;

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
sys.params_.track_max_missed_detection_time_ = 30;
 
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
        model.missed_detection_time_ = 40;
    } else {
        model.missed_detection_time_ = 1;
    }
    
    model_manager.AddModel(sys, model);

    
}

sys.models_.back().model_likelihood_ = 0.1;

// Make a model similar to model 8 so that it will be merged. Give this model the better label, missed detection time and likelihood
model.model_likelihood_ = 20;
model.missed_detection_time_ = 0;
model.state_.g_.data_ << 8.5, 27.5, 10.5;
model.state_.u_.data_ << 4.5,8.5,16.5;
model.err_cov_ = Eigen::Matrix<double,6,6>::Identity()*2;
model.label_ = 0;
sys.model_label_++;

model_manager.AddModel(sys, model);

// Get the models that will be merged together
merge2 = model;
auto iter = sys.models_.begin();
iter = std::next(iter,8);
merge1 = *iter;

model_manager.ManageModels(sys,0.7);

// Get the merged model
iter = sys.models_.begin();
iter = std::next(iter,5);

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

ASSERT_EQ(iter->missed_detection_time_,merge2.missed_detection_time_);
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
typedef ModelRN<State, TransformNULL,SourceRN> Model;
typedef Meas<double> Measurement;

System<Model> sys;
ModelManager<Model> model_manager;  

SourceR3 source;
SourceParameters source_params;
source_params.meas_cov_ = Eigen::Matrix3d::Identity();
source_params.spacial_density_of_false_meas_ = 0.1;
source_params.type_ = MeasurementTypes::RN_POS;
source_params.probability_of_detection_ = 0.9;
source_params.gate_probability_ = 0.8;
source_params.source_index_ = 0;
source.Init(source_params);

sys.sources_.push_back(source);

sys.params_.track_max_num_tracks_ = 2;
sys.params_.process_noise_covariance_ = Eigen::Matrix<double,6,6>::Identity();

Measurement m;
m.source_index = 0;
m.weight = 1;
m.likelihood = 1;
m.time_stamp = 0;
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
ASSERT_DOUBLE_EQ(sys.models_.begin()->missed_detection_time_, 1.0);

model_manager.UpdateModels(sys);


ASSERT_EQ(sys.models_.front().missed_detection_time_, 0);
ASSERT_EQ(sys.models_.back().missed_detection_time_, 0);


}

// Test the transform model

TEST(ModelManagerTest, TransformModelTest ) {

typedef lie_groups::R2_r2 State;
typedef ModelRN<State, TransformHomography,SourceRN> Model;
typedef Meas<double> Measurement;

TransformHomography<State> trans;

// The homography data is set so that the state and the measurements will all be transformed to zero
Eigen::Matrix3d homography;
homography.setZero();
homography.block(2,2,1,1)<< 1;
trans.Init();
trans.SetData(homography);

System<Model> sys;
ModelManager<Model> model_manager;  

SourceR2 source;
SourceParameters source_params;
source_params.meas_cov_ = Eigen::Matrix2d::Identity();
source_params.spacial_density_of_false_meas_ = 0.1;
source_params.type_ = MeasurementTypes::RN_POS;
source_params.probability_of_detection_ = 0.9;
source_params.gate_probability_ = 0.8;
source_params.source_index_ = 0;
source.Init(source_params);

sys.sources_.push_back(source);

sys.params_.track_max_num_tracks_ = 2;
sys.params_.process_noise_covariance_ = Eigen::Matrix<double,4,4>::Identity();
sys.params_.transform_consensus_set_ = true;
sys.transformaion_ = trans;

Measurement m;
m.source_index = 0;
m.weight = 1;
m.likelihood = 1;
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


} // namespace rransac



