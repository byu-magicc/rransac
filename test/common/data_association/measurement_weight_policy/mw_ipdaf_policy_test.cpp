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

#include "rransac/common/data_association/measurement_weight_policies/mw_ipdaf_policy.h"
#include "rransac/common/data_association/validation_region_policies/validation_region_fixed_policy.h"
#include "rransac/common/data_association/track_likelihood_info_policies/tli_ipdaf_policy.h"
#include "rransac/common/data_association/data_association_host.h"

using namespace lie_groups;
using namespace rransac;

typedef ModelRN<R2_r2, TransformNULL, SourceRN> Model;

TEST(TargetLikelihoodUpdateTest, IPDAF) {



// Setup sources
double meas_noise0 = 0.013;
double meas_noise1 = 0.028;
double PG0 = 0.90;
double PD0 = 0.70;
double lambda0 = 0.14;
double PG1 = 0.80;
double PD1 = 0.85;
double lambda1 = 0.11;
SourceParameters source_params0, source_params1;
source_params0.spacial_density_of_false_meas_ = lambda0;
source_params0.gate_probability_ = PG0;
source_params0.probability_of_detection_ = PD0;
source_params0.source_index_ = 0;
source_params0.type_ = MeasurementTypes::RN_POS;
source_params0.meas_cov_ = Eigen::Matrix2d::Identity()*meas_noise0;
source_params1.spacial_density_of_false_meas_ = lambda1;
source_params1.gate_probability_ = PG1;
source_params1.probability_of_detection_ = PD1;
source_params1.source_index_ = 1;
source_params1.type_ = MeasurementTypes::RN_POS_VEL;
source_params1.meas_cov_ = Eigen::Matrix<double,4,4>::Identity()*meas_noise1;

Model::Source source0, source1;
source0.Init(source_params0);
source1.Init(source_params1);

// Setup parameters
Parameters params;
double process_noise = 0.31;
params.process_noise_covariance_ = Eigen::Matrix<double,4,4>::Identity()*process_noise;

// Setup tracks
Model track0, track1;
double track0_init_model_likelihood = 0.5;
double track1_init_model_likelihood = 0.7;
track0.Init(params,2);
track0.state_.g_.data_ << 0,0;
track0.state_.u_.data_ << 0,0;
track0.err_cov_ = Eigen::Matrix<double,4,4>::Identity()*2;
track0.model_likelihood_ = track0_init_model_likelihood;
track0.newest_measurement_time_stamp = 0;
track1.Init(params,2);
track1.state_.g_.data_ << 1,-1;
track1.state_.u_.data_ <<-1, 1;
track1.err_cov_ = Eigen::Matrix<double,4,4>::Identity()*2;
track1.model_likelihood_ = track1_init_model_likelihood;
track1.newest_measurement_time_stamp = 0;


// Setup system
System<Model> sys;
sys.sources_.push_back(source0);
sys.sources_.push_back(source1);
sys.models_.push_back(track0);
sys.models_.push_back(track1);
sys.current_time_ = 0.1;
sys.dt_ = sys.current_time_ - track1.newest_measurement_time_stamp;

auto model_iter0 = sys.models_.begin();
auto model_iter1 = std::next(model_iter0);

// Setup measurements
Meas<double> m0, m1;
m0.type = source_params0.type_;
m0.time_stamp = sys.current_time_;
m0.source_index = 0;
m1.type = source_params1.type_;
m1.time_stamp = sys.current_time_;
m1.source_index = 1;
std::vector<Meas<double>> meas_base {m0,m1};
std::vector<std::list<Model>::iterator> model_iter = {model_iter0, model_iter1};

// Setup data association info
DataAssociationInfo info;
info.source_produced_measurements_.push_back(true);
info.source_produced_measurements_.push_back(true);

model_iter0->model_likelihood_update_info_[0].in_lsr_and_produced_meas = true;
model_iter0->model_likelihood_update_info_[0].num_assoc_meas = 4;
model_iter0->model_likelihood_update_info_[0].delta = 0.5;

model_iter0->model_likelihood_update_info_[1].in_lsr_and_produced_meas = true;
model_iter0->model_likelihood_update_info_[1].num_assoc_meas = 4;
model_iter0->model_likelihood_update_info_[1].delta = 0.4;

model_iter1->model_likelihood_update_info_[0].in_lsr_and_produced_meas = true;
model_iter1->model_likelihood_update_info_[0].num_assoc_meas = 4;
model_iter1->model_likelihood_update_info_[0].delta = 0.7;

model_iter1->model_likelihood_update_info_[1].in_lsr_and_produced_meas = true;
model_iter1->model_likelihood_update_info_[1].num_assoc_meas = 4;
model_iter1->model_likelihood_update_info_[1].delta = 0.1;

// Add the measurements and calculate the weights.

std::vector<std::vector<std::vector<Meas<double>>>> meas(2,std::vector<std::vector<Meas<double>>>(2, std::vector<Meas<double>>(4)));

for (int model_index =0; model_index < 2; ++model_index) {
    for (int source_index =0; source_index <2; ++ source_index) {
        for (int meas_index =0; meas_index < 4; ++ meas_index) {
            meas[model_index][source_index][meas_index] = meas_base[source_index];
            meas[model_index][source_index][meas_index].probability = (meas_index+1)/10;
            model_iter[model_index]->AddNewMeasurement(meas[model_index][source_index][meas_index] );
            meas[model_index][source_index][meas_index].weight = sys.sources_[source_index].params_.probability_of_detection_*meas[model_index][source_index][meas_index].probability/sys.sources_[source_index].params_.spacial_density_of_false_meas_/(1.0-model_iter[model_index]->model_likelihood_update_info_[source_index].delta);
        }
    }
}



// Setup measurement weight policy and check compatibility
MW_IPDAFPolicy<Model> mw_policy;

static_assert(MW_IPDAFPolicy<Model>::CompatiblityCheck<ValidationRegionFixedPolicy,TLI_IPDAFPolicy>::value,"MW_IPDAFPolicy::CompatiblityCheck The policy MW_IPDAFPolicy is only compatible with the update track likelihood policy TLI_IPDAFPolicy." );

mw_policy.PolicyCalculateMeasurementWeight(sys,info);

// Check measurement weights;
for (int model_index =0; model_index < 2; ++model_index) {
    for (int source_index =0; source_index <2; ++ source_index) {
        for (int meas_index =0; meas_index < 4; ++ meas_index) {
            ASSERT_DOUBLE_EQ(meas[model_index][source_index][meas_index].weight , model_iter[model_index]->new_assoc_meas_[source_index][meas_index].weight) << "model index: " << model_index << " source index: " << source_index << " meas index: " << meas_index;
        }
    }
}

}
