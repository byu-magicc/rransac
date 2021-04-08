#include <gtest/gtest.h>
#include <typeinfo>
#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>


#include "rransac/common/models/model_base.h"
#include "rransac/common/models/model_RN.h"
#include "rransac/common/models/model_SEN_pos_vel.h"
#include "rransac/common/models/model_SEN_pose_twist.h"
#include "rransac/common/sources/source_base.h"
#include "rransac/common/sources/source_RN.h"
#include "rransac/common/sources/source_SEN_pos_vel.h"
#include "rransac/common/sources/source_SEN_pose_twist.h"
#include "rransac/parameters.h"
#include "rransac/common/transformations/transformation_base.h"
#include "rransac/common/transformations/trans_homography.h"
#include "rransac/common/transformations/transformation_null.h"


namespace rransac
{

using namespace lie_groups;

typedef ModelRN<R2_r2, TransformNULL, SourceRN> Model1;
// typedef ModelRN<R3_r3, TransformNULL, SourceRN> Model2;
// typedef ModelSENPosVel<SE2_se2, TransformNULL, SourceSENPosVel> Model3;
// typedef ModelSENPoseTwist<SE2_se2, TransformNULL, SourceSENPoseTwist> Model4;
// typedef ModelSENPosVel<SE3_se3, TransformNULL, SourceSENPosVel> Model5;
// typedef ModelSENPoseTwist<SE3_se3, TransformNULL, SourceSENPoseTwist> Model6;



template<class M, MeasurementTypes MT1, MeasurementTypes MT2>
struct ModelHelper {
    typedef M Model;
    static constexpr MeasurementTypes MeasType1 = MT1;
    static constexpr MeasurementTypes MeasType2 = MT2;
};

typedef ModelHelper<Model1, MeasurementTypes::RN_POS, MeasurementTypes::RN_POS_VEL> ModelHelper1;
// typedef ModelHelper<Model2, MeasurementTypes::RN_POS, MeasurementTypes::RN_POS_VEL> ModelHelper2;
// typedef ModelHelper<Model3, MeasurementTypes::SEN_POS, MeasurementTypes::SEN_POS_VEL> ModelHelper3;
// typedef ModelHelper<Model4, MeasurementTypes::SEN_POSE, MeasurementTypes::SEN_POSE_TWIST> ModelHelper4;
// typedef ModelHelper<Model5, MeasurementTypes::SEN_POS, MeasurementTypes::SEN_POS_VEL> ModelHelper5;
// typedef ModelHelper<Model6, MeasurementTypes::SEN_POSE, MeasurementTypes::SEN_POSE_TWIST> ModelHelper6;


// using MyTypes = ::testing::Types<ModelHelper1, ModelHelper2, ModelHelper3, ModelHelper4, ModelHelper5, ModelHelper6 >;
// using MyTypes = ::testing::Types<ModelHelper3, ModelHelper5 >;
using MyTypes = ::testing::Types< ModelHelper1>;

/////////////////////////////////////////////////////////////////////////////////////////////


template<class Model>
class ModelTest : public ::testing::Test {

protected:

typedef Meas<double> Measurement;

static constexpr unsigned int meas_dim = Model::Model::Source::meas_space_dim_;
static constexpr unsigned int state_dim = Model::Model::State::g_type_::dim_*2;
static constexpr unsigned int cov_dim = Model::Model::cov_dim_;
static constexpr unsigned int a_vel_dim = Model::Model::cov_dim_ - Model::Model::State::g_type_::dim_-1;
static constexpr unsigned int t_vel_dim = state_dim - Model::Model::State::g_type_::dim_ - a_vel_dim;

void SetUp() override {

// This grabs the desired components of the F and G matrices.
if (Model::MeasType1 == MeasurementTypes::SEN_POS|| Model::MeasType1 == MeasurementTypes::SEN_POS_VEL) {
    Filter.setZero();
    Filter.block(0,0, g_dim_+1, g_dim_+1).setIdentity();
    Filter.block(g_dim_+1, g_dim_ + t_vel_dim,  a_vel_dim, a_vel_dim).setIdentity();
// std::cout << Filter << std::endl;
} else {
    Filter.setIdentity();
}



// std::cerr << "setup" << std::endl;
// Setup the parameters
double meas_cov_scale = 0.01;
double system_cov_scale = 0.1;
meas_cov1.setIdentity();
meas_cov1 *= meas_cov_scale;
meas_std1.setIdentity();
meas_std1 *= sqrt(meas_cov_scale);
// source_params1.meas_cov_fixed_ = false;
source_params1.meas_cov_ = meas_cov1;
source_params1.type_ = Model::MeasType1;
source_params1.source_index_ = 0;
source_params1.spacial_density_of_false_meas_ = 0.8;
source_params1.probability_of_detection_ = 0.8;
source_params1.gate_probability_ = 0.8;

meas_cov2.setIdentity();
meas_cov2 *= meas_cov_scale;
// source_params2.meas_cov_fixed_ = true;
// source_params2.meas_cov_ = meas_cov2;
source_params2.type_ = Model::MeasType2;
source_params2.source_index_ = 1;
source_params2.spacial_density_of_false_meas_ = 0.8;
source_params2.probability_of_detection_ = 0.8;
source_params2.gate_probability_ = 0.8;
source_params2.meas_cov_ = meas_cov2;
meas_std2.setIdentity();
meas_std2 *= sqrt(meas_cov_scale);

// std::cerr << "here0" << std::endl;
typename Model::Model::Source source1;
typename Model::Model::Source source2;
source1.Init(source_params1);
source2.Init(source_params2);
// std::cerr << "here00" << std::endl;

sources.push_back(source1);
sources.push_back(source2);
// std::cerr << source1.gsd_ptr_ << std::endl;
// std::cerr << sources.front().gsd_ptr_ << std::endl;
params.process_noise_covariance_ = Eigen::Matrix<double,cov_dim,cov_dim>::Identity()*system_cov_scale;

// std::cerr << "here 1" << std::endl;

// Generate Trajectory and measurements
Measurement m1;
Measurement m2;
std::vector<std::vector<Measurement>> gen_meas(2);
std::vector<Measurement> meas1(num_meas);
std::vector<Measurement> meas2(num_meas);

state = Model::Model::GetRandomState();



state0 = state;                           // Initial state
states.push_back(state);                

// std::cerr << "here 2" << std::endl;


for (unsigned int ii = 0; ii < num_iters; ++ii) {

    // Propagate state
    state.g_.OPlusEq(state.u_.data_*dt);
    states.push_back(state); 

    // Get measurements
    for (unsigned int jj=0; jj < num_meas; ++jj) {

        m1 = source1.GenerateRandomMeasurement(state, meas_std1);
        m2 = source2.GenerateRandomMeasurement(state, meas_std2);
        // m2.meas_cov = meas_cov2;
        m1.weight = 1.0/(num_meas+1);
        m2.weight = 1.0/(num_meas+1);
        m1.time_stamp = ii*dt;
        m2.time_stamp = ii*dt;
        meas1[jj] = m1;
        meas2[jj] = m2;
    }
    gen_meas[0] = meas1;
    gen_meas[1] = meas2;
    new_meas.push_back(gen_meas);

}

// std::cerr << "here 3" << std::endl;
// Construct Jacobians with state0 and dt
Eigen::Matrix<double, g_dim_*2, g_dim_*2> F_tmp, G_tmp;
F_tmp.block(0,0,g_dim_, g_dim_) = typename Model::Model::State::g_type_(Model::Model::State::u_type_::Exp(-state0.u_.data_*dt)).Adjoint();
F_tmp.block(0,g_dim_,g_dim_,g_dim_) = (state0.u_*dt).Jr()*dt;
F_tmp.block(g_dim_,0,g_dim_,g_dim_).setZero();
F_tmp.block(g_dim_,g_dim_,g_dim_,g_dim_).setIdentity(); 

// std::cout << "f_tmp: " << std::endl << F_tmp << std::endl << std::endl;


G_tmp.block(0,0,g_dim_, g_dim_) = (state0.u_*dt).Jr();
G_tmp.block(0,g_dim_,g_dim_,g_dim_) = (state0.u_*dt).Jr()/2;
G_tmp.block(g_dim_,0,g_dim_,g_dim_).setZero(); 
G_tmp.block(g_dim_,g_dim_,g_dim_,g_dim_)= Eigen::Matrix<double,Model::Model::g_dim_,Model::Model::g_dim_>::Identity();

F = Filter*F_tmp*Filter.transpose();
G = Filter*G_tmp*Filter.transpose();

}

//-----------------------------------------------------------------------------------------------
typedef typename Model::Model model;
SourceParameters source_params1;
SourceParameters source_params2;
Eigen::Matrix<double,meas_dim, meas_dim> meas_cov1;
Eigen::Matrix<double,meas_dim*2, meas_dim*2> meas_cov2;
Eigen::Matrix<double,meas_dim, meas_dim> meas_std1;
Eigen::Matrix<double,meas_dim*2, meas_dim*2> meas_std2;
Parameters params;

typename Model::Model::State state;
typename Model::Model::State state0;
std::vector<typename Model::Model::State> states;

std::vector<typename Model::Model::Source> sources;

std::vector<std::vector<std::vector<Measurement>>> new_meas; // time, source, measurements

int num_meas = 2;
int num_iters = 10;
double dt = 0.1;

typename Model::Model track;
typename Model::Model::Mat F;
typename Model::Model::Mat G;
static constexpr unsigned int g_dim_ = Model::Model::g_dim_;
static constexpr unsigned int cov_dim_ = Model::Model::cov_dim_;
Eigen::Matrix<double, cov_dim_, g_dim_*2> Filter;

};


TYPED_TEST_SUITE(ModelTest, MyTypes);


/////////////////////////////////////////////////////////////////////////////////////////////////////////////


TYPED_TEST(ModelTest, AllFunctions) {

// Test init function and set parameters function
this->track.Init(this->params,this->sources.size());
ASSERT_EQ(this->track.Q_, this->params.process_noise_covariance_);
ASSERT_EQ(this->track.err_cov_, TypeParam::Model::Mat::Identity());
ASSERT_EQ(this->track.newest_measurement_time_stamp, 0);
ASSERT_EQ(this->track.model_likelihood_, 0.5);
ASSERT_EQ(this->track.label_, -1);
ASSERT_EQ(this->track.innov_cov_set_.size(), this->sources.size());
ASSERT_EQ(this->track.innovation_covariances_.size(), this->sources.size());
ASSERT_EQ(this->track.new_assoc_meas_.size(), this->sources.size());
ASSERT_EQ(this->track.model_likelihood_update_info_.size(), this->sources.size());




// Test propagate state
this->state = this->state0;
this->state = this->track.PropagateState(this->state, this->dt*this->num_iters);
ASSERT_LE( (this->state.g_.data_- this->states.back().g_.data_).norm(), 1e-12);
ASSERT_EQ(this->state.u_.data_, this->states.back().u_.data_);


// Test the system function Jacobians
// std::cout << "this->F: " << std::endl << this->F << std::endl << std::endl;
// std::cout << "track->F: " << std::endl << this->track.GetLinTransFuncMatState(this->state0,this->dt) << std::endl;
// std::cout << "err: " << std::endl << this->track.GetLinTransFuncMatState(this->state0,this->dt) - this->F << std::endl;
// std::cout << "F: " << std::endl << this->track.GetLinTransFuncMatState(this->state0,this->dt) << std::endl;
// std::cout << "Ft: " << std::endl << this->F << std::endl << std::endl;
ASSERT_LE( ( this->track.GetLinTransFuncMatState(this->state0,this->dt) - this->F).norm(), 1e-12);
ASSERT_LE( ( this->track.GetLinTransFuncMatNoise(this->state0,this->dt) - this->G).norm(), 1e-12);

// Test the Propagate Model Function
std::fill(this->track.innov_cov_set_.begin(), this->track.innov_cov_set_.end(), true); // make sure that this is reset to false.
this->track.state_ = this->state0;
this->track.PropagateModel(this->dt);
for (auto val : this->track.innov_cov_set_) {
    ASSERT_FALSE(val);
}

// ASSERT_EQ(this->track.missed_detection_time_, this->dt);
ASSERT_LE( (this->track.F_- this->F).norm(), 1e-12  );
ASSERT_LE( (this->track.G_- this->G).norm(), 1e-12  );
ASSERT_LE( (this->track.state_.g_.data_ - this->states[1].g_.data_).norm(), 1e-12  );
ASSERT_LE( (this->track.state_.u_.data_ - this->states[1].u_.data_).norm(), 1e-12  );


ASSERT_LE( (this->track.err_cov_ - this->F*this->F.transpose() - this->G*this->track.Q_*this->dt*this->G.transpose()).norm(), 1e-12  );

// Test the innovation covariance function.
////////////////////////////////////////////////
Eigen::MatrixXd H0, H1, V0, V1, S0, S1, S0_e, S1_e; // Observation function jacobians w.r.t. the state and noise for both sources and innovation covariances.
H0 = this->track.GetLinObsMatState(this->sources,this->track.state_,0);
V0 = this->track.GetLinObsMatSensorNoise(this->sources,this->track.state_,0);
H1 = this->track.GetLinObsMatState(this->sources,this->track.state_,1);
V1 = this->track.GetLinObsMatSensorNoise(this->sources,this->track.state_,1);
S0 = H0*this->track.err_cov_*H0.transpose() + V0*this->sources[0].params_.meas_cov_*V0.transpose();
S1 = H1*this->track.err_cov_*H1.transpose() + V1*this->sources[1].params_.meas_cov_*V1.transpose();
S0_e = this->track.GetInnovationCovariance(this->sources,0);
ASSERT_EQ(this->track.innov_cov_set_.size(),this->sources.size());
ASSERT_EQ(this->track.innovation_covariances_.size(),this->sources.size());
ASSERT_TRUE(this->track.innov_cov_set_[0]);  // Make sure the innovation covariance for the first source was set and none other.
ASSERT_GT(this->track.innovation_covariances_[0].rows(),0);
for(int ii = 1; ii < this->sources.size(); ++ii) {
    ASSERT_FALSE(this->track.innov_cov_set_[ii]);
    ASSERT_EQ(this->track.innovation_covariances_[ii].rows(),0);
}
ASSERT_LE( (S0- S0_e).norm(), 1e-12  );

S1_e = this->track.GetInnovationCovariance(this->sources,1);
ASSERT_TRUE(this->track.innov_cov_set_[1]);  // Make sure the innovation covariance for the first source was set and none other.
ASSERT_GT(this->track.innovation_covariances_[1].rows(),0);
ASSERT_LE( (S1- S1_e).norm(), 1e-12  );

this->track.UpdateModel(this->sources, this->params); // Make sure that the update model resets the innov_cov_set to false
for (auto val : this->track.innov_cov_set_) {
    ASSERT_FALSE(val);
}


// Test the update model function
this->track.err_cov_.setIdentity();
this->track.state_ = this->track.state_.Identity();





// std::cerr << "state0: " << this->track.state_.g_.data_ << std::endl;
// std::cerr <<  this->track.state_.u_.data_ << std::endl;
double model_likelihood_prev;
ModelLikelihoodUpdateInfo info;


for (unsigned long int ii=0; ii < this->num_iters; ++ii) {

    this->track.PropagateModel(this->dt);
    this->track.new_assoc_meas_ = this->new_meas[ii];
    model_likelihood_prev = this->track.model_likelihood_;

    this->track.UpdateModel(this->sources, this->params);


    

}
ASSERT_DOUBLE_EQ(this->track.newest_measurement_time_stamp,this->num_iters*this->dt - this->dt);

// Test covariance
ASSERT_LE(  ((this->track.err_cov_ + this->track.err_cov_.transpose())/2.0 - this->track.err_cov_ ).norm(), 1e-12); // symmetric
Eigen::VectorXcd eigen_values = this->track.err_cov_.eigenvalues();
for (int ii =0; ii < eigen_values.rows(); ++ii){                       // positive definite
    ASSERT_GT( std::real(eigen_values(ii)), 0);
}

if (TypeParam::MeasType1 == MeasurementTypes::SEN_POS|| TypeParam::MeasType1 == MeasurementTypes::SEN_POS_VEL) {
    
    // Eigen::MatrixXd error = (*this->track.sources_)[1].params_.meas_cov_.inverse()*(*this->track.sources_)[1].OMinus( this->new_meas.back().back().back(),(*this->track.sources_)[1].GetEstMeas(this->track.state_));
    // std::cerr << "norm error: " << error.norm() << std::endl;
    // std::cerr << " cov inverse: " << this->meas_cov2.inverse() << std::endl;
    // std::cerr << " blah er: " << (*this->track.sources_)[1].OMinus( this->new_meas.back().back().back(),(*this->track.sources_)[1].GetEstMeas(this->track.state_))<< std::endl;
    Eigen::MatrixXd error = this->meas_cov2.inverse().sqrt()*this->sources[1].OMinus( this->new_meas.back().back().back(),this->sources[1].GetEstMeas(this->track.state_));
    EXPECT_LE( error.norm(),  5) << "This result might not always be true since updating the track is not deterministic";


} else {
    EXPECT_LE( (this->Filter.transpose()*this->track.err_cov_.inverse().sqrt()*this->Filter*this->track.state_.OMinus(this->states.back())).norm(),  5) << "This result might not always be true since updating the track is not deterministic";
    // std::cout << "norm error " << (this->Filter.transpose()*this->track.err_cov_.inverse().sqrt()*this->Filter*this->track.state_.OMinus(this->states.back())).norm() << std::endl;

}

// 

// std::cout << "g_est: " << std::endl << this->track.state_.g_.data_ << std::endl;
// std::cout << "g_true: " << std::endl << this->states.back().g_.data_<< std::endl;
// std::cout << "u_est: " << std::endl << this->track.state_.u_.data_ << std::endl;
// std::cout << "u_true: " << std::endl << this->states.back().u_.data_ << std::endl;
// std::cout << "cov: " << std::endl << this->track.err_cov_ << std::endl;
// std::cout << "norm error " << (this->track.state_.OMinus(this->states.back())).norm() << std::endl;


// Verify the Consensus Set
ASSERT_EQ(this->track.cs_.Size(), this->num_iters);

ConsensusSet<Meas<double>> set;
set.consensus_set_.begin();

for( std::list<std::vector<Meas<double>>>::iterator it = this->track.cs_.consensus_set_.begin(); it!= this->track.cs_.consensus_set_.end(); ++it) {
    ASSERT_EQ( (*it).size(), 4);
}

this->track.PruneConsensusSet(0);
ASSERT_EQ(this->track.cs_.Size(), this->num_iters-1);




}






} // namespace rransac
