#include <gtest/gtest.h>
#include <Eigen/Core>
#include <math.h>
#include <iostream>


#include "rransac/common/measurement/measurement_base.h"

#include "rransac/common/models/model_RN.h"
#include "rransac/common/sources/source_RN.h"
#include "rransac/parameters.h"
#include "rransac/common/transformations/transformation_null.h"
#include "rransac/system.h"
#include "lie_groups/state.h"
#include "rransac/common/data_association/data_association_host.h"
#include "rransac/common/sources/source_container.h"

using namespace lie_groups;
using namespace rransac;

typedef SourceRN<R2_r2,MeasurementTypes::RN_POS,TransformNULL> SourcePos;
typedef SourceRN<R2_r2,MeasurementTypes::RN_POS_VEL,TransformNULL> SourcePosVel;
typedef SourceContainer<SourcePos,SourcePosVel> SC;
typedef ModelRN<SC> Model;

const int knum_assoc =100;
const double kdelta = 0.6;
const double kmeas_prob = 0.5;
const double kmeas_weight = 0.1;

// Create dummy policies
template<typename tModel>
class ValidationRegionDummyPolicy {
public:
static bool PolicyInValidationRegion(const System<tModel>& sys, const Meas<typename tModel::DataType>& meas, tModel& track)  {return true;}
};

//-------------------------------------------------------------------------------------------------------------------

template<typename tModel>
class TLI_DummyPolicy {
public:
static void PolicyUpdateTrackLikelihood(System<tModel>& sys,DataAssociationInfo& info )  {
    
    for (auto& track : sys.models_) {

          for (int source_index =0; source_index < track.model_likelihood_update_info_.size(); ++source_index) {
              track.model_likelihood_update_info_[source_index].in_lsr_and_produced_meas = true;
              track.model_likelihood_update_info_[source_index].num_assoc_meas = knum_assoc;
              track.model_likelihood_update_info_[source_index].delta = kdelta;

              for (auto& meas: track.new_assoc_meas_[source_index]) {
                  meas.probability = kmeas_prob;
              }

          }

    }
}

};

//---------------------------------------------------------------------

template<typename tModel>
class MW_DummyPolicy {

public:
    static void PolicyCalculateMeasurementWeight(System<tModel>& sys, DataAssociationInfo& info) {

        for (auto& track: sys.models_) {
            for (int source_index =0; source_index < sys.source_container_.num_sources_; ++source_index) {
                for(auto& meas : track.new_assoc_meas_[source_index]) {
                    meas.weight = kmeas_weight;
                }
            }
        }

    }

    template<template<class> typename tValidationRegionPolicy, template<class> typename tUpdateTrackLikelihoodPolicy>
    struct CompatiblityCheck {
        static constexpr bool value = true;
        
    }; 



};



TEST(TargetLikelihoodUpdateTest, IPDAF) {



// Setup sources
double meas_noise0 = 0.013;
double meas_noise1 = 0.028;
SourceParameters source_params0, source_params1;
source_params0.source_index_ = 0;
source_params0.type_ = MeasurementTypes::RN_POS;
source_params0.meas_cov_ = Eigen::Matrix2d::Identity()*meas_noise0;
source_params1.source_index_ = 1;
source_params1.type_ = MeasurementTypes::RN_POS_VEL;
source_params1.meas_cov_ = Eigen::Matrix<double,4,4>::Identity()*meas_noise1;

// Setup parameters
Parameters params;
double process_noise = 0.31;
params.process_noise_covariance_ = Eigen::Matrix<double,4,4>::Identity()*process_noise;

// Setup tracks
Model track0, track1;
track0.Init(params);
track1.Init(params);


// Setup system
System<Model> sys;
sys.source_container_.AddSource(source_params0);
sys.source_container_.AddSource(source_params1);
sys.models_.push_back(track0);
sys.models_.push_back(track1);
sys.current_time_ = 0.1;
sys.dt_ = 0.1;

auto model_iter0 = sys.models_.begin();
auto model_iter1 = std::next(model_iter0);


// Setup measurements
Meas<double> m0, m1;
m0.type = source_params0.type_;
m0.time_stamp = sys.current_time_;
m0.source_index = 0;
m0.probability = 0;
m0.weight = 0;
m1.type = source_params1.type_;
m1.time_stamp = sys.current_time_;
m1.source_index = 1;
m1.probability = 0;
m1.weight = 0;



for (long int ii =0; ii < knum_assoc; ++ii) {
    sys.new_meas_.push_back(m0);
}

for (long int ii =0; ii < knum_assoc; ++ii) {
    sys.new_meas_.push_back(m1);
}

DataAssociationHost<Model,ValidationRegionDummyPolicy,TLI_DummyPolicy,MW_DummyPolicy> data_association_host;

data_association_host.AssociateNewMeasurements(sys);

// All of the measurements shouldve been added.
ASSERT_EQ(model_iter0->new_assoc_meas_[0].size(), knum_assoc);
ASSERT_EQ(model_iter0->new_assoc_meas_[1].size(), knum_assoc);
ASSERT_EQ(model_iter1->new_assoc_meas_[0].size(), knum_assoc);
ASSERT_EQ(model_iter1->new_assoc_meas_[1].size(), knum_assoc);

for (int source_index =0; source_index < 2; ++source_index) {

    ASSERT_DOUBLE_EQ(model_iter0->model_likelihood_update_info_[source_index].delta,kdelta);
    ASSERT_EQ(model_iter0->model_likelihood_update_info_[source_index].num_assoc_meas,knum_assoc);
    ASSERT_TRUE(model_iter0->model_likelihood_update_info_[source_index].in_lsr_and_produced_meas);

    for (auto& meas : model_iter0->new_assoc_meas_[source_index]) {
        meas.weight = kmeas_weight;
        meas.probability = kmeas_prob;
    }

    for (auto& meas : model_iter1->new_assoc_meas_[source_index]) {
        meas.weight = kmeas_weight;
        meas.probability = kmeas_prob;
    }


}




}