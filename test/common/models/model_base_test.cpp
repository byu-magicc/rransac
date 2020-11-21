#include <gtest/gtest.h>
#include <typeinfo>
#include <Eigen/Core>


#include "common/models/model_base.h"
#include "common/models/model_RN.h"
#include "common/sources/source_base.h"
#include "common/sources/source_RN.h"
#include "common/sources/source_SEN_pos_vel.h"
#include "common/sources/source_SEN_pose_twist.h"
#include "parameters.h"
#include "common/transformations/transformation_base.h"
#include "common/transformations/trans_homography.h"
#include "common/transformations/transformation_null.h"


namespace rransac
{

using namespace lie_groups;

typedef ModelRN<R2_r2, TransformNULL<R2_r2>> Model1;
typedef ModelRN<R3_r3, TransformNULL<R3_r3>> Model2;
// typedef ModelBase<SE2_se2, SE2_se2::StateType, SourceSENPosVel<SE2_se2>,    TransformNULL<SE2_se2>> Model3;
// typedef ModelBase<SE2_se2, SE2_se2::StateType, SourceSENPoseTwist<SE2_se2>, TransformNULL<SE2_se2>> Model4;
// typedef ModelBase<SE3_se3, SE3_se3::StateType, SourceSENPosVel<SE3_se3>,    TransformNULL<SE3_se3>> Model5;
// typedef ModelBase<SE3_se3, SE3_se3::StateType, SourceSENPoseTwist<SE3_se3>, TransformNULL<SE3_se3>> Model6;

template<class M, MeasurementTypes MT1, MeasurementTypes MT2>
struct ModelHelper {
    typedef M Model;
    static constexpr MeasurementTypes MeasType1 = MT1;
    static constexpr MeasurementTypes MeasType2 = MT2;
};

typedef ModelHelper<Model1, MeasurementTypes::RN_POS, MeasurementTypes::RN_POS_VEL> ModelHelper1;
typedef ModelHelper<Model2, MeasurementTypes::RN_POS, MeasurementTypes::RN_POS_VEL> ModelHelper2;
// typedef ModelHelper<Model3, MeasurementTypes::SEN_POS, MeasurementTypes::SEN_POS_VEL> ModelHelper3;
// typedef ModelHelper<Model4, MeasurementTypes::SEN_POSE, MeasurementTypes::SEN_POSE_TWIST> ModelHelper4;
// typedef ModelHelper<Model5, MeasurementTypes::SEN_POS, MeasurementTypes::SEN_POS_VEL> ModelHelper5;
// typedef ModelHelper<Model6, MeasurementTypes::SEN_POSE, MeasurementTypes::SEN_POSE_TWIST> ModelHelper6;


// using MyTypes = ::testing::Types<ModelHelper1, ModelHelper2, ModelHelper3, ModelHelper4, ModelHelper5, ModelHelper6 >;
using MyTypes = ::testing::Types< ModelHelper1>;

/////////////////////////////////////////////////////////////////////////////////////////////


template<class Model>
class ModelTest : public ::testing::Test {

protected:

static constexpr unsigned int meas_dim = Model::Model::MB_Source::meas_dim;
static constexpr unsigned int state_dim = Model::Model::MB_State::g_type_::dim_*2;

void SetUp() override {

// std::cerr << "setup" << std::endl;
// Setup the parameters
double meas_cov_scale = 0.01;
double system_cov_scale = 0.01;
meas_cov1.setIdentity();
meas_cov1 *= meas_cov_scale;
meas_std1.setIdentity();
meas_std1 *= sqrt(meas_cov_scale);
source_params1.meas_cov_fixed_ = true;
// source_params1.meas_cov_fixed_ = false;
source_params1.meas_cov_ = meas_cov1;
source_params1.type_ = Model::MeasType1;
source_params1.source_index_ = 0;

source_params2.meas_cov_fixed_ = false;
meas_cov2.setIdentity();
meas_cov2 *= meas_cov_scale;
// source_params2.meas_cov_fixed_ = true;
// source_params2.meas_cov_ = meas_cov2;
source_params2.type_ = Model::MeasType2;
source_params2.source_index_ = 1;
meas_std2.setIdentity();
meas_std2 *= sqrt(meas_cov_scale);

// std::cerr << "here0" << std::endl;
typename Model::Model::MB_Source source1;
typename Model::Model::MB_Source source2;
source1.Init(source_params1);
source2.Init(source_params2);
// std::cerr << "here00" << std::endl;

sources.push_back(source1);
sources.push_back(source2);
// std::cerr << source1.gsd_ptr_ << std::endl;
// std::cerr << sources.front().gsd_ptr_ << std::endl;
params.process_noise_covariance_ = Eigen::Matrix<double,state_dim,state_dim>::Identity()*system_cov_scale;

// std::cerr << "here 1" << std::endl;

// Generate Trajectory and measurements
Meas m1;
Meas m2;
std::vector<std::vector<Meas>> gen_meas(2);
std::vector<Meas> meas1(num_meas);
std::vector<Meas> meas2(num_meas);

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
        m2.meas_cov = meas_cov2;
        m1.weight = 1.0/(num_meas+1);
        m2.weight = 1.0/(num_meas+1);
        meas1[jj] = m1;
        meas2[jj] = m2;
    }
    gen_meas[0] = meas1;
    gen_meas[1] = meas2;
    new_meas.push_back(gen_meas);

}

// std::cerr << "here 3" << std::endl;
// Construct Jacobians with state0 and dt
F.block(0,0,g_dim_, g_dim_) = typename Model::Model::MB_State::g_type_(Model::Model::MB_State::u_type_::Exp(state0.u_.data_*dt)).Adjoint();
F.block(0,g_dim_,g_dim_,g_dim_) = (state0.u_*dt).Jr()*dt;
F.block(g_dim_,0,g_dim_,g_dim_).setZero();
F.block(g_dim_,g_dim_,g_dim_,g_dim_).setIdentity(); 

G.block(0,0,g_dim_, g_dim_) = (state0.u_*dt).Jr()*dt;
G.block(0,g_dim_,g_dim_,g_dim_) = (state0.u_*dt).Jr()*dt*dt/2;
G.block(g_dim_,0,g_dim_,g_dim_).setZero(); 
G.block(g_dim_,g_dim_,g_dim_,g_dim_)= Eigen::Matrix<double,Model::Model::g_dim_,Model::Model::g_dim_>::Identity()*dt;

}

//-----------------------------------------------------------------------------------------------

SourceParameters source_params1;
SourceParameters source_params2;
Eigen::Matrix<double,meas_dim, meas_dim> meas_cov1;
Eigen::Matrix<double,meas_dim*2, meas_dim*2> meas_cov2;
Eigen::Matrix<double,meas_dim, meas_dim> meas_std1;
Eigen::Matrix<double,meas_dim*2, meas_dim*2> meas_std2;
Parameters params;

typename Model::Model::MB_State state;
typename Model::Model::MB_State state0;
std::vector<typename Model::Model::MB_State> states;

std::vector<typename Model::Model::MB_Source> sources;

std::vector<std::vector<std::vector<Meas>>> new_meas; // time, source, measurements

int num_meas = 2;
int num_iters = 10;
double dt = 0.1;

typename Model::Model track;
typename Model::Model::Mat F;
typename Model::Model::Mat G;
static constexpr unsigned int g_dim_ = Model::Model::g_dim_;

};


TYPED_TEST_SUITE(ModelTest, MyTypes);


/////////////////////////////////////////////////////////////////////////////////////////////////////////////


TYPED_TEST(ModelTest, AllFunctions) {

// Test init function and set parameters function
this->track.Init(this->sources, this->params);
ASSERT_EQ((*this->track.sources_)[0].params_.source_index_, this->sources[0].params_.source_index_);
ASSERT_EQ((*this->track.sources_)[1].params_.source_index_, 1);
ASSERT_EQ(this->track.Q_, this->params.process_noise_covariance_);
ASSERT_EQ(this->track.err_cov_, TypeParam::Model::Mat::Identity());

// Test propagate state
this->state = this->state0;
this->state = this->track.PropagateState(this->state, this->dt*this->num_iters);
ASSERT_LE( (this->state.g_.data_- this->states.back().g_.data_).norm(), 1e-12);
ASSERT_EQ(this->state.u_.data_, this->states.back().u_.data_);


// Test the system function Jacobians
ASSERT_EQ(this->track.GetLinTransFuncMatState(this->state0,this->dt), this->F);
ASSERT_EQ(this->track.GetLinTransFuncMatNoise(this->state0,this->dt), this->G);

// Test the Propagate Model Function
this->track.state_ = this->state0;
this->track.PropagateModel(this->dt);
// std::cout << "err_cov: " << this->track.err_cov_ << std::endl << std::endl;
// std::cout << "G_: " << this->track.G_ << std::endl << std::endl;
ASSERT_LE( (this->track.F_- this->F).norm(), 1e-12  );
ASSERT_LE( (this->track.G_- this->G).norm(), 1e-12  );
ASSERT_LE( (this->track.state_.g_.data_ - this->states[1].g_.data_).norm(), 1e-12  );
ASSERT_LE( (this->track.state_.u_.data_ - this->states[1].u_.data_).norm(), 1e-12  );
ASSERT_LE( (this->track.err_cov_ - this->F*this->F.transpose() - this->G*this->track.Q_*this->G.transpose()).norm(), 1e-12  );

// Test the update model function
this->track.err_cov_.setIdentity();
this->track.state_ = this->track.state_.Identity();

std::cerr << "state0: " << this->track.state_.g_.data_ << std::endl;
std::cerr <<  this->track.state_.u_.data_ << std::endl;

for (unsigned long int ii=0; ii < this->num_iters; ++ii) {

    this->track.PropagateModel(this->dt);
    this->track.new_assoc_meas_ = this->new_meas[ii];
    this->track.UpdateModel(this->params);

}

std::cout << "g_est: " << std::endl << this->track.state_.g_.data_ << std::endl;
std::cout << "g_true: " << std::endl << this->states.back().g_.data_<< std::endl;
std::cout << "u_est: " << std::endl << this->track.state_.u_.data_ << std::endl;
std::cout << "u_true: " << std::endl << this->states.back().u_.data_ << std::endl;
std::cout << "cov: " << std::endl << this->track.err_cov_ << std::endl;


}




} // namespace rransac
