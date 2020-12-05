#include <gtest/gtest.h>
#include <Eigen/Core>

#include "system.h"
#include "common/models/model_manager.h"
#include "common/transformations/transformation_null.h"
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
typedef ModelRN<State, TransformNULL<State>> Model;

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
typedef ModelRN<State, TransformNULL<State>> Model;

System<Model> sys;
ModelManager<Model> model_manager;  

Model model;
model.state_.g_.data_ << 1,1,1;
model.state_.u_.data_ << 2,3,4;

sys.params_.max_num_models_ = 2;
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
typedef ModelRN<State, TransformNULL<State>> Model;

System<Model> sys;
ModelManager<Model> model_manager;
Model model; 

model.err_cov_.setIdentity();

sys.params_.max_num_models_ = 10;
sys.params_.good_model_threshold_ = 5;
sys.params_.similar_tracks_threshold_ = 2;
sys.params_.max_missed_detection_time_ = 40;
 
// setup the models

for(int ii = 0; ii < sys.params_.max_num_models_+3; ++ii) {
    
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
    }
    
    model_manager.AddModel(sys, model);

    
}

model.model_likelihood_ = 0.1;


model_manager.AddModel(sys, model);


ASSERT_EQ(sys.models_.size(), max_num_models);

auto iter = sys.models_.begin();

ASSERT_EQ((*iter).label_, 1 );
++iter;
ASSERT_EQ((*iter).label_, 2 );
++iter;
ASSERT_EQ((*iter).label_, 3 );
++iter;
ASSERT_EQ((*iter).label_, 4 );
++iter;
ASSERT_EQ((*iter).label_, 6 );
++iter;
ASSERT_EQ((*iter).label_, 7 );
++iter;
ASSERT_EQ((*iter).label_, 8 );
++iter;
ASSERT_EQ((*iter).label_, 9 );
++iter;
ASSERT_EQ((*iter).label_, 10 );
++iter;
ASSERT_EQ((*iter).label_, 11 );
}

// ----------------------------------------------------------------

TEST(ModelManagerTest, PruneConsensusSet ) {

typedef lie_groups::R3_r3 State;
typedef ModelRN<State, TransformNULL<State>> Model;

System<Model> sys;
ModelManager<Model> model_manager;  

Model model;

sys.params_.max_num_models_ = 2;

Meas m1, m2;
m1.time_stamp = 0.5; 
m2.time_stamp = 1;

std::vector<Meas> meas{m1,m2};

model.cs_.AddMeasurementsToConsensusSet(meas);
model_manager.AddModel(sys, model);
model_manager.AddModel(sys, model);
model_manager.PruneConsensusSets(sys, 0.7);

ASSERT_EQ(sys.models_.begin()->cs_.consensus_set_.size(), 1);
ASSERT_EQ(sys.models_.back().cs_.consensus_set_.size(), 1);

}




/**
 *  Ensures that UpdateModel is called on all of the models
 */ 
TEST(ModelManagerTest, UpdateModel ) {

typedef lie_groups::R3_r3 State;
typedef ModelRN<State, TransformNULL<State>> Model;

System<Model> sys;
ModelManager<Model> model_manager;  

SourceR3 source;
SourceParameters source_params;
source_params.meas_cov_fixed_ = true;
source_params.meas_cov_ = Eigen::Matrix3d::Identity();
source_params.expected_num_false_meas_ = 0.1;
source_params.type_ = MeasurementTypes::RN_POS;
source_params.probability_of_detection_ = 0.9;
source_params.gate_probability_ = 0.8;
source_params.source_index_ = 0;
source.Init(source_params);

sys.sources_.push_back(source);

sys.params_.max_num_models_ = 2;
sys.params_.process_noise_covariance_ = Eigen::Matrix<double,6,6>::Identity();

Meas m;
m.source_index = 0;
m.weight = 1;
m.likelihood = 1;
m.time_stamp = 0;
m.type = MeasurementTypes::RN_POS;
m.pose = Eigen::Matrix<double,3,1>::Random();


Model model;
model.Init(sys.sources_,sys.params_);
model.state_.g_.data_ << 1,1,1;
model.state_.u_.data_ << 2,3,4;

model.new_assoc_meas_ = std::vector<std::vector<Meas>>{std::vector<Meas>{m}};

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

TEST(ModelManagerTest, MergeModels ) {

typedef lie_groups::R3_r3 State;
typedef ModelRN<State, TransformNULL<State>> Model;

System<Model> sys;
ModelManager<Model> model_manager; 
Model model1, model2;

Meas m1,m2;
m1.time_stamp = 0;
m2.time_stamp = 1;

// Make model2 with better stats than model1

model1.label_ = 2;
model1.missed_detection_time_ = 10;
model1.model_likelihood_ = 10;
model1.err_cov_.setZero();
model1.err_cov_.diagonal() << 2.5, 2.5, 2.5, 1.5, 1.5, 1.5;
model1.cs_.AddMeasToConsensusSet(m1);

model2.label_ = 3;
model2.missed_detection_time_ = 1;
model2.model_likelihood_ = 1000;
model2.err_cov_.setZero();
model2.err_cov_.diagonal() << 0.5, 0.5, 0.5, 3,3,3;
model2.cs_.AddMeasToConsensusSet(m2);

// Set states so that they are not similar
model1.state_.g_.data_ << 1,1,1;
model2.state_.g_.data_ << 5,5,5;
model1.state_.u_.data_ << 1,1,1;
model2.state_.u_.data_ << 1,1,1;

sys.params_.max_num_models_ = 3;
sys.params_.similar_tracks_threshold_ = 1;

model_manager.AddModel(sys, model1);
model_manager.AddModel(sys, model2);

// The models are not similar to merge so no models should be merged.
model_manager.MergeModels(sys);

ASSERT_EQ(sys.models_.size(),2);

// change the data of model2 so that it will be merged
model2.state_.g_.data_ << 1.5, 1.5, 1.5;
model2.state_.u_.data_ << 1.1, 1.1, 1.1;
model2.label_ = 1;
model_manager.AddModel(sys, model2);


model_manager.MergeModels(sys);

ASSERT_EQ(sys.models_.back().label_,3);
ASSERT_EQ(sys.models_.size(),2);
ASSERT_EQ(sys.models_.front().label_, model2.label_);
ASSERT_EQ(sys.models_.front().missed_detection_time_, model2.missed_detection_time_);
ASSERT_EQ(sys.models_.front().model_likelihood_, model2.model_likelihood_);

// The fused state should be between the two original states
ASSERT_LT(sys.models_.front().state_.g_.data_(0), model2.state_.g_.data_(0));
ASSERT_LT(sys.models_.front().state_.g_.data_(1), model2.state_.g_.data_(1));
ASSERT_LT(sys.models_.front().state_.g_.data_(2), model2.state_.g_.data_(2));
ASSERT_LT(sys.models_.front().state_.u_.data_(0), model2.state_.u_.data_(0));
ASSERT_LT(sys.models_.front().state_.u_.data_(1), model2.state_.u_.data_(1));
ASSERT_LT(sys.models_.front().state_.u_.data_(2), model2.state_.u_.data_(2));

ASSERT_GT(sys.models_.front().state_.g_.data_(0), model1.state_.g_.data_(0));
ASSERT_GT(sys.models_.front().state_.g_.data_(1), model1.state_.g_.data_(1));
ASSERT_GT(sys.models_.front().state_.g_.data_(2), model1.state_.g_.data_(2));
ASSERT_GT(sys.models_.front().state_.u_.data_(0), model1.state_.u_.data_(0));
ASSERT_GT(sys.models_.front().state_.u_.data_(1), model1.state_.u_.data_(1));
ASSERT_GT(sys.models_.front().state_.u_.data_(2), model1.state_.u_.data_(2));

// The fused error covariance should be smaller than the other two covariances
ASSERT_LT(sys.models_.front().err_cov_(0,0), model2.err_cov_(0,0));
ASSERT_LT(sys.models_.front().err_cov_(1,1), model2.err_cov_(1,1));
ASSERT_LT(sys.models_.front().err_cov_(2,2), model2.err_cov_(2,2));
ASSERT_LT(sys.models_.front().err_cov_(3,3), model2.err_cov_(3,3));
ASSERT_LT(sys.models_.front().err_cov_(4,4), model2.err_cov_(4,4));
ASSERT_LT(sys.models_.front().err_cov_(5,5), model2.err_cov_(5,5));

ASSERT_LT(sys.models_.front().err_cov_(0,0), model1.err_cov_(0,0));
ASSERT_LT(sys.models_.front().err_cov_(1,1), model1.err_cov_(1,1));
ASSERT_LT(sys.models_.front().err_cov_(2,2), model1.err_cov_(2,2));
ASSERT_LT(sys.models_.front().err_cov_(3,3), model1.err_cov_(3,3));
ASSERT_LT(sys.models_.front().err_cov_(4,4), model1.err_cov_(4,4));
ASSERT_LT(sys.models_.front().err_cov_(5,5), model1.err_cov_(5,5));

ASSERT_GT(sys.models_.front().err_cov_(0,0), 0);
ASSERT_GT(sys.models_.front().err_cov_(1,1), 0);
ASSERT_GT(sys.models_.front().err_cov_(2,2), 0);
ASSERT_GT(sys.models_.front().err_cov_(3,3), 0);
ASSERT_GT(sys.models_.front().err_cov_(4,4), 0);
ASSERT_GT(sys.models_.front().err_cov_(5,5), 0);



// std::cout << "state g" << sys.models_.front().state_.g_.data_ << std::endl;
// std::cout << "state u" << sys.models_.front().state_.u_.data_ << std::endl;
// std::cout << "cov" << sys.models_.front().err_cov_ << std::endl;

}

} // namespace rransac



