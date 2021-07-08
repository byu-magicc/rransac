#include <gtest/gtest.h>
#include <Eigen/Core>
#include <math.h>
#include <iostream>


#include "rransac/common/measurement/measurement_base.h"

#include "rransac/common/models/model_RN.h"
#include "rransac/common/sources/source_RN.h"
#include "rransac/common/sources/source_container.h"
#include "rransac/parameters.h"
#include "rransac/common/transformations/transformation_null.h"
#include "rransac/system.h"
#include "lie_groups/state.h"

#include "rransac/common/data_association/track_likelihood_info_policies/tli_ipdaf_policy.h"
#include "rransac/common/data_association/data_association_host.h"

using namespace lie_groups;
using namespace rransac;


typedef SourceRN<R2_r2,MeasurementTypes::RN_POS,TransformNULL> SourceR2Pos;
typedef SourceRN<R2_r2,MeasurementTypes::RN_POS_VEL,TransformNULL> SourceR2PosVel;
typedef SourceContainer<SourceR2Pos,SourceR2PosVel> SourceContainerR2;
typedef ModelRN<SourceContainerR2> Model;
typedef typename SourceR2Pos::Measurement Measurement;
typedef typename SourceR2Pos::TransformDataType TransformDataType;

static bool StateInsideSurveillanceRegionFalseCallback(const R2_r2& state) {
    return false;
}

double CalculateMeasurementLikelihood(const System<Model>& sys, Model& track, const Measurement& meas) {

    Eigen::MatrixXd err = sys.source_container_.OMinus(meas.source_index, meas, sys.source_container_.GetEstMeas(meas.source_index,track.state_,meas.transform_state,meas.transform_data_t_m));
    Eigen::MatrixXd S = track.GetInnovationCovariance(sys.source_container_,meas.source_index, meas.transform_state, meas.transform_data_t_m);
    double det_inn_cov_sqrt = sqrt(S.determinant());

    return exp( - (err.transpose()*S.inverse()*err)(0,0)/2.0)   /(pow(2.0*M_PI,S.cols()/2.0)*det_inn_cov_sqrt) ;

}

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



// Setup parameters
Parameters params;
double process_noise = 0.31;
params.process_noise_covariance_ = Eigen::Matrix<double,4,4>::Identity()*process_noise;

// Setup tracks
Model track0, track1;
double track0_init_model_likelihood = 0.5;
double track1_init_model_likelihood = 0.7;
track0.Init(params);
track0.state_.g_.data_ << 0,0;
track0.state_.u_.data_ << 0,0;
track0.err_cov_ = Eigen::Matrix<double,4,4>::Identity()*2;
track0.model_likelihood_ = track0_init_model_likelihood;
track0.newest_measurement_time_stamp = 0;
track1.Init(params);
track1.state_.g_.data_ << 1,-1;
track1.state_.u_.data_ <<-1, 1;
track1.err_cov_ = Eigen::Matrix<double,4,4>::Identity()*2;
track1.model_likelihood_ = track1_init_model_likelihood;
track1.newest_measurement_time_stamp = 0;


// Setup system
System<Model> sys;
sys.source_container_.AddSource(source_params0);
sys.source_container_.AddSource(source_params1,StateInsideSurveillanceRegionFalseCallback);
sys.models_.push_back(track0);
sys.models_.push_back(track1);
sys.current_time_ = 0.1;
sys.dt_ = sys.current_time_ - track1.newest_measurement_time_stamp;

auto model_iter0 = sys.models_.begin();
auto model_iter1 = std::next(model_iter0);

// Setup measurements
Measurement m0, m1;
m0.type = source_params0.type_;
m0.time_stamp = sys.current_time_;
m0.source_index = 0;
m1.type = source_params1.type_;
m1.time_stamp = sys.current_time_;
m1.source_index = 1;

// Setup data association info
DataAssociationInfo<TransformDataType> info;
TransformDataType EmptyMat;
int num_sources = SourceContainerR2::num_sources_;
for (int ii = 0; ii < num_sources; ++ii) {
    info.source_produced_measurements_.push_back(false);
    info.transform_state_.push_back(false);
    info.transform_data_t_m_.push_back(EmptyMat);
}



// Verify that the model parameters are set up
ASSERT_EQ(sys.models_.front().innov_cov_set_.size(), num_sources) << "The size of innov_cov_set_ should always be the size of sys.sources_";
ASSERT_EQ(sys.models_.front().innovation_covariances_.size(), num_sources)<< "The size of innovation_covariances_ should always be the size of sys.sources_";
ASSERT_EQ(sys.models_.front().new_assoc_meas_.size(), num_sources)<< "The size of new_assoc_meas_ should always be the size of sys.sources_";
ASSERT_EQ(sys.models_.front().model_likelihood_update_info_.size(),num_sources)<< "The size of model_likelihood_update_info_ should always be the size of sys.sources_";

TLI_IPDAFPolicy<Model> tli_policy;

//Since there are no new measurements test that the likelihood decays with time
tli_policy.PolicyUpdateTrackLikelihood(sys,info);
track0_init_model_likelihood *= (1.0-tli_policy.decay_rate_*sys.dt_);
track1_init_model_likelihood *= (1.0-tli_policy.decay_rate_*sys.dt_);
ASSERT_DOUBLE_EQ(model_iter0->model_likelihood_, track0_init_model_likelihood);
ASSERT_DOUBLE_EQ(model_iter1->model_likelihood_, track1_init_model_likelihood);

ASSERT_FALSE(model_iter0->model_likelihood_update_info_[0].in_lsr_and_produced_meas);
ASSERT_DOUBLE_EQ(model_iter0->model_likelihood_update_info_[0].num_assoc_meas,0);
ASSERT_DOUBLE_EQ(model_iter0->model_likelihood_update_info_[0].delta,0);
ASSERT_FALSE(model_iter1->model_likelihood_update_info_[0].in_lsr_and_produced_meas);
ASSERT_DOUBLE_EQ(model_iter1->model_likelihood_update_info_[0].num_assoc_meas,0);
ASSERT_DOUBLE_EQ(model_iter1->model_likelihood_update_info_[0].delta,0);
ASSERT_FALSE(model_iter0->model_likelihood_update_info_[1].in_lsr_and_produced_meas);
ASSERT_DOUBLE_EQ(model_iter0->model_likelihood_update_info_[1].num_assoc_meas,0);
ASSERT_DOUBLE_EQ(model_iter0->model_likelihood_update_info_[1].delta,0);
ASSERT_FALSE(model_iter1->model_likelihood_update_info_[1].in_lsr_and_produced_meas);
ASSERT_DOUBLE_EQ(model_iter1->model_likelihood_update_info_[1].num_assoc_meas,0);
ASSERT_DOUBLE_EQ(model_iter1->model_likelihood_update_info_[1].delta,0);

// Verify the update when there are no measurements, the target is inside the surveillance region and the source produced measurements
// The tracks are only in the surveillance region of the first source.
std::fill(info.source_produced_measurements_.begin(), info.source_produced_measurements_.end(), true);
tli_policy.PolicyUpdateTrackLikelihood(sys,info);
track0_init_model_likelihood *= (1.0-tli_policy.decay_rate_*sys.dt_);
track1_init_model_likelihood *= (1.0-tli_policy.decay_rate_*sys.dt_);
track0_init_model_likelihood = (1.0 - PD0*PG0)*track0_init_model_likelihood/(1.0-PD0*PG0*track0_init_model_likelihood);
track1_init_model_likelihood = (1.0 - PD0*PG0)*track1_init_model_likelihood/(1.0-PD0*PG0*track1_init_model_likelihood);
ASSERT_DOUBLE_EQ(model_iter0->model_likelihood_, track0_init_model_likelihood);
ASSERT_DOUBLE_EQ(model_iter1->model_likelihood_, track1_init_model_likelihood);

ASSERT_TRUE(model_iter0->model_likelihood_update_info_[0].in_lsr_and_produced_meas);
ASSERT_DOUBLE_EQ(model_iter0->model_likelihood_update_info_[0].num_assoc_meas,0);
ASSERT_DOUBLE_EQ(model_iter0->model_likelihood_update_info_[0].delta,PG0*PD0);
ASSERT_TRUE(model_iter1->model_likelihood_update_info_[0].in_lsr_and_produced_meas);
ASSERT_DOUBLE_EQ(model_iter1->model_likelihood_update_info_[0].num_assoc_meas,0);
ASSERT_DOUBLE_EQ(model_iter1->model_likelihood_update_info_[0].delta,PG0*PD0);
ASSERT_FALSE(model_iter0->model_likelihood_update_info_[1].in_lsr_and_produced_meas);
ASSERT_DOUBLE_EQ(model_iter0->model_likelihood_update_info_[1].num_assoc_meas,0);
ASSERT_DOUBLE_EQ(model_iter0->model_likelihood_update_info_[1].delta,0);
ASSERT_FALSE(model_iter1->model_likelihood_update_info_[1].in_lsr_and_produced_meas);
ASSERT_DOUBLE_EQ(model_iter1->model_likelihood_update_info_[1].num_assoc_meas,0);
ASSERT_DOUBLE_EQ(model_iter1->model_likelihood_update_info_[1].delta,0);

// Change the source1 so that both models are in its surveillance region. However, specify that only the first produced measurements
sys.source_container_ = SourceContainerR2();
sys.source_container_.AddSource(source_params0);
sys.source_container_.AddSource(source_params1);
info.source_produced_measurements_[0] = true;
info.source_produced_measurements_[1] = false;
tli_policy.PolicyUpdateTrackLikelihood(sys,info);
track0_init_model_likelihood *= (1.0-tli_policy.decay_rate_*sys.dt_);
track1_init_model_likelihood *= (1.0-tli_policy.decay_rate_*sys.dt_);
track0_init_model_likelihood = (1.0 - PD0*PG0)*track0_init_model_likelihood/(1.0-PD0*PG0*track0_init_model_likelihood);
track1_init_model_likelihood = (1.0 - PD0*PG0)*track1_init_model_likelihood/(1.0-PD0*PG0*track1_init_model_likelihood);
ASSERT_DOUBLE_EQ(model_iter0->model_likelihood_, track0_init_model_likelihood);
ASSERT_DOUBLE_EQ(model_iter1->model_likelihood_, track1_init_model_likelihood);

ASSERT_TRUE(model_iter0->model_likelihood_update_info_[0].in_lsr_and_produced_meas);
ASSERT_DOUBLE_EQ(model_iter0->model_likelihood_update_info_[0].num_assoc_meas,0);
ASSERT_DOUBLE_EQ(model_iter0->model_likelihood_update_info_[0].delta,PG0*PD0);
ASSERT_TRUE(model_iter1->model_likelihood_update_info_[0].in_lsr_and_produced_meas);
ASSERT_DOUBLE_EQ(model_iter1->model_likelihood_update_info_[0].num_assoc_meas,0);
ASSERT_DOUBLE_EQ(model_iter1->model_likelihood_update_info_[0].delta,PG0*PD0);
ASSERT_FALSE(model_iter0->model_likelihood_update_info_[1].in_lsr_and_produced_meas);
ASSERT_DOUBLE_EQ(model_iter0->model_likelihood_update_info_[1].num_assoc_meas,0);
ASSERT_DOUBLE_EQ(model_iter0->model_likelihood_update_info_[1].delta,0);
ASSERT_FALSE(model_iter1->model_likelihood_update_info_[1].in_lsr_and_produced_meas);
ASSERT_DOUBLE_EQ(model_iter1->model_likelihood_update_info_[1].num_assoc_meas,0);
ASSERT_DOUBLE_EQ(model_iter1->model_likelihood_update_info_[1].delta,0);

// Specify that both produced measurements
sys.source_container_ = SourceContainerR2();
sys.source_container_.AddSource(source_params0);
sys.source_container_.AddSource(source_params1);
info.source_produced_measurements_[0] = true;
info.source_produced_measurements_[1] = true;
tli_policy.PolicyUpdateTrackLikelihood(sys,info);
track0_init_model_likelihood *= (1.0-tli_policy.decay_rate_*sys.dt_);
track1_init_model_likelihood *= (1.0-tli_policy.decay_rate_*sys.dt_);
track0_init_model_likelihood = (1.0 - PD0*PG0)*track0_init_model_likelihood/(1.0-PD0*PG0*track0_init_model_likelihood);
track1_init_model_likelihood = (1.0 - PD0*PG0)*track1_init_model_likelihood/(1.0-PD0*PG0*track1_init_model_likelihood);
track0_init_model_likelihood = (1.0 - PD1*PG1)*track0_init_model_likelihood/(1.0-PD1*PG1*track0_init_model_likelihood);
track1_init_model_likelihood = (1.0 - PD1*PG1)*track1_init_model_likelihood/(1.0-PD1*PG1*track1_init_model_likelihood);
ASSERT_DOUBLE_EQ(model_iter0->model_likelihood_, track0_init_model_likelihood);
ASSERT_DOUBLE_EQ(model_iter1->model_likelihood_, track1_init_model_likelihood);

ASSERT_TRUE(model_iter0->model_likelihood_update_info_[0].in_lsr_and_produced_meas);
ASSERT_DOUBLE_EQ(model_iter0->model_likelihood_update_info_[0].num_assoc_meas,0);
ASSERT_DOUBLE_EQ(model_iter0->model_likelihood_update_info_[0].delta,PG0*PD0);
ASSERT_TRUE(model_iter1->model_likelihood_update_info_[0].in_lsr_and_produced_meas);
ASSERT_DOUBLE_EQ(model_iter1->model_likelihood_update_info_[0].num_assoc_meas,0);
ASSERT_DOUBLE_EQ(model_iter1->model_likelihood_update_info_[0].delta,PG0*PD0);
ASSERT_TRUE(model_iter0->model_likelihood_update_info_[1].in_lsr_and_produced_meas);
ASSERT_DOUBLE_EQ(model_iter0->model_likelihood_update_info_[1].num_assoc_meas,0);
ASSERT_DOUBLE_EQ(model_iter0->model_likelihood_update_info_[1].delta,PG1*PD1);
ASSERT_TRUE(model_iter1->model_likelihood_update_info_[1].in_lsr_and_produced_meas);
ASSERT_DOUBLE_EQ(model_iter1->model_likelihood_update_info_[1].num_assoc_meas,0);
ASSERT_DOUBLE_EQ(model_iter1->model_likelihood_update_info_[1].delta,PG1*PD1);

// Reset the track likelihoods
track0_init_model_likelihood = 0.51;
model_iter0->model_likelihood_ = track0_init_model_likelihood;
track1_init_model_likelihood = 0.63;
model_iter1->model_likelihood_ = track1_init_model_likelihood;

// Let the first source produce measurements and not the second
info.source_produced_measurements_[0] = true;
info.source_produced_measurements_[1] = false;
m0.pose = model_iter0->state_.g_.data_;
model_iter0->AddNewMeasurement(m0);

m0.pose = model_iter1->state_.g_.data_;
model_iter1->AddNewMeasurement(m0);

Eigen::Matrix2d S00, S10;
Eigen::Matrix<double,4,4> S01, S11;
double p00, p10,p01,p11;
double delta00, delta10, delta01, delta11;

S00 = model_iter0->err_cov_.block(0,0,2,2) + sys.source_container_.GetParams(0).meas_cov_;
S10 = model_iter1->err_cov_.block(0,0,2,2) + sys.source_container_.GetParams(0).meas_cov_;
p00 = pow(2*M_PI, -1)*pow(S00.determinant(),-0.5);
p10 = pow(2*M_PI,-1)*pow(S10.determinant(),-0.5);

delta00 = PD0*PG0 - PD0/lambda0*p00;
delta10 = PD0*PG0 - PD0/lambda0*p10;

tli_policy.PolicyUpdateTrackLikelihood(sys,info);
track0_init_model_likelihood *= (1.0-tli_policy.decay_rate_*sys.dt_);
track1_init_model_likelihood *= (1.0-tli_policy.decay_rate_*sys.dt_);
track0_init_model_likelihood = (1.0-delta00)*track0_init_model_likelihood/(1-delta00*track0_init_model_likelihood);
track1_init_model_likelihood = (1.0-delta10)*track1_init_model_likelihood/(1-delta10*track1_init_model_likelihood);

ASSERT_DOUBLE_EQ(model_iter0->model_likelihood_, track0_init_model_likelihood);
ASSERT_DOUBLE_EQ(model_iter1->model_likelihood_, track1_init_model_likelihood);

ASSERT_TRUE(model_iter0->model_likelihood_update_info_[0].in_lsr_and_produced_meas);
ASSERT_DOUBLE_EQ(model_iter0->model_likelihood_update_info_[0].num_assoc_meas,1);
ASSERT_DOUBLE_EQ(model_iter0->model_likelihood_update_info_[0].delta,delta00);
ASSERT_TRUE(model_iter1->model_likelihood_update_info_[0].in_lsr_and_produced_meas);
ASSERT_DOUBLE_EQ(model_iter1->model_likelihood_update_info_[0].num_assoc_meas,1);
ASSERT_DOUBLE_EQ(model_iter1->model_likelihood_update_info_[0].delta,delta10);
ASSERT_FALSE(model_iter0->model_likelihood_update_info_[1].in_lsr_and_produced_meas);
ASSERT_DOUBLE_EQ(model_iter0->model_likelihood_update_info_[1].num_assoc_meas,0);
ASSERT_DOUBLE_EQ(model_iter0->model_likelihood_update_info_[1].delta,0);
ASSERT_FALSE(model_iter1->model_likelihood_update_info_[1].in_lsr_and_produced_meas);
ASSERT_DOUBLE_EQ(model_iter1->model_likelihood_update_info_[1].num_assoc_meas,0);
ASSERT_DOUBLE_EQ(model_iter1->model_likelihood_update_info_[1].delta,0);

// Let both produce measurements but only the first have associated measurements. Also have dt=0
info.source_produced_measurements_[0] = true;
info.source_produced_measurements_[1] = true;
sys.dt_ = 0;
sys.current_time_ = 0.1;

tli_policy.PolicyUpdateTrackLikelihood(sys,info);
track0_init_model_likelihood = (1.0-delta00)*track0_init_model_likelihood/(1.0-delta00*track0_init_model_likelihood);
track1_init_model_likelihood = (1.0-delta10)*track1_init_model_likelihood/(1.0-delta10*track1_init_model_likelihood);
track0_init_model_likelihood = (1.0 - PD1*PG1)*track0_init_model_likelihood/(1.0-PD1*PG1*track0_init_model_likelihood);
track1_init_model_likelihood = (1.0 - PD1*PG1)*track1_init_model_likelihood/(1.0-PD1*PG1*track1_init_model_likelihood);

ASSERT_DOUBLE_EQ(model_iter0->model_likelihood_, track0_init_model_likelihood);
ASSERT_DOUBLE_EQ(model_iter1->model_likelihood_, track1_init_model_likelihood);

ASSERT_TRUE(model_iter0->model_likelihood_update_info_[0].in_lsr_and_produced_meas);
ASSERT_DOUBLE_EQ(model_iter0->model_likelihood_update_info_[0].num_assoc_meas,1);
ASSERT_DOUBLE_EQ(model_iter0->model_likelihood_update_info_[0].delta,delta00);
ASSERT_TRUE(model_iter1->model_likelihood_update_info_[0].in_lsr_and_produced_meas);
ASSERT_DOUBLE_EQ(model_iter1->model_likelihood_update_info_[0].num_assoc_meas,1);
ASSERT_DOUBLE_EQ(model_iter1->model_likelihood_update_info_[0].delta,delta10);
ASSERT_TRUE(model_iter0->model_likelihood_update_info_[1].in_lsr_and_produced_meas);
ASSERT_DOUBLE_EQ(model_iter0->model_likelihood_update_info_[1].num_assoc_meas,0);
ASSERT_DOUBLE_EQ(model_iter0->model_likelihood_update_info_[1].delta,PD1*PG1);
ASSERT_TRUE(model_iter1->model_likelihood_update_info_[1].in_lsr_and_produced_meas);
ASSERT_DOUBLE_EQ(model_iter1->model_likelihood_update_info_[1].num_assoc_meas,0);
ASSERT_DOUBLE_EQ(model_iter1->model_likelihood_update_info_[1].delta,PD1*PG1);

// Let both produce associated measurements with dt=0 and Reset the track likelihoods

track0_init_model_likelihood = 0.51;
model_iter0->model_likelihood_ = track0_init_model_likelihood;
track1_init_model_likelihood = 0.63;
model_iter1->model_likelihood_ = track1_init_model_likelihood;

info.source_produced_measurements_[0] = true;
info.source_produced_measurements_[1] = true;
m0.pose << 1.1, 1.3;
model_iter0->AddNewMeasurement(m0);

m0.pose << -0.5, 0.1;
model_iter1->AddNewMeasurement(m0);

m1.pose = Eigen::Matrix<double,2,1>::Random();
m1.twist = Eigen::Matrix<double,2,1>::Random();
model_iter0->AddNewMeasurement(m1);
m1.pose = Eigen::Matrix<double,2,1>::Random();
m1.twist = Eigen::Matrix<double,2,1>::Random();
model_iter0->AddNewMeasurement(m1);

m1.pose = Eigen::Matrix<double,2,1>::Random();
m1.twist = Eigen::Matrix<double,2,1>::Random();
model_iter1->AddNewMeasurement(m1);

m1.pose = Eigen::Matrix<double,2,1>::Random();
m1.twist = Eigen::Matrix<double,2,1>::Random();
model_iter1->AddNewMeasurement(m1);


p00 = CalculateMeasurementLikelihood( sys, *model_iter0, model_iter0->new_assoc_meas_[0].front());
p01 = CalculateMeasurementLikelihood( sys, *model_iter0, model_iter0->new_assoc_meas_[0].back());

delta00 = PD0*PG0 - PD0/lambda0*(p00+p01);

p00 = CalculateMeasurementLikelihood( sys, *model_iter0, model_iter0->new_assoc_meas_[1].front());
p01 = CalculateMeasurementLikelihood( sys, *model_iter0, model_iter0->new_assoc_meas_[1].back());

delta01 = PD1*PG1 - PD1/lambda1*(p00+p01);

p10 = CalculateMeasurementLikelihood( sys, *model_iter1, model_iter1->new_assoc_meas_[0].front());
p11 = CalculateMeasurementLikelihood( sys, *model_iter1, model_iter1->new_assoc_meas_[0].back());

delta10 = PD0*PG0 - PD0/lambda0*(p10+p11);

p10 = CalculateMeasurementLikelihood( sys, *model_iter1, model_iter1->new_assoc_meas_[1].front());
p11 = CalculateMeasurementLikelihood( sys, *model_iter1, model_iter1->new_assoc_meas_[1].back());

delta11 = PD1*PG1 - PD1/lambda1*(p10+p11);



tli_policy.PolicyUpdateTrackLikelihood(sys,info);
track0_init_model_likelihood = (1.0-delta00)*track0_init_model_likelihood/(1.0-delta00*track0_init_model_likelihood);
track1_init_model_likelihood = (1.0-delta10)*track1_init_model_likelihood/(1.0-delta10*track1_init_model_likelihood);
track0_init_model_likelihood = (1.0-delta01)*track0_init_model_likelihood/(1.0-delta01*track0_init_model_likelihood);
track1_init_model_likelihood = (1.0-delta11)*track1_init_model_likelihood/(1.0-delta11*track1_init_model_likelihood);

ASSERT_TRUE(model_iter0->model_likelihood_update_info_[0].in_lsr_and_produced_meas);
ASSERT_DOUBLE_EQ(model_iter0->model_likelihood_update_info_[0].num_assoc_meas,2);
ASSERT_DOUBLE_EQ(model_iter0->model_likelihood_update_info_[0].delta,delta00);
ASSERT_TRUE(model_iter1->model_likelihood_update_info_[0].in_lsr_and_produced_meas);
ASSERT_DOUBLE_EQ(model_iter1->model_likelihood_update_info_[0].num_assoc_meas,2);
ASSERT_DOUBLE_EQ(model_iter1->model_likelihood_update_info_[0].delta,delta10);
ASSERT_TRUE(model_iter0->model_likelihood_update_info_[1].in_lsr_and_produced_meas);
ASSERT_DOUBLE_EQ(model_iter0->model_likelihood_update_info_[1].num_assoc_meas,2);
ASSERT_DOUBLE_EQ(model_iter0->model_likelihood_update_info_[1].delta,delta01);
ASSERT_TRUE(model_iter1->model_likelihood_update_info_[1].in_lsr_and_produced_meas);
ASSERT_DOUBLE_EQ(model_iter1->model_likelihood_update_info_[1].num_assoc_meas,2);
ASSERT_DOUBLE_EQ(model_iter1->model_likelihood_update_info_[1].delta,delta11);

}