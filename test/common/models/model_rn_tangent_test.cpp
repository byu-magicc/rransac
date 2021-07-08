#include <gtest/gtest.h>
#include <typeinfo>
#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

#include "rransac/common/models/model_base.h"
#include "rransac/common/models/model_RN.h"
#include "rransac/common/sources/source_base.h"
#include "rransac/common/sources/source_RN.h"
#include "rransac/parameters.h"
#include "rransac/common/transformations/transformation_null.h"
#include "rransac/common/sources/source_container.h"
#include <unsupported/Eigen/MatrixFunctions>

namespace rransac {

using namespace lie_groups;



// Create state types
typedef State<Rn,double,1,1> StateR1_1;
typedef State<Rn,double,1,2> StateR1_2;
typedef State<Rn,double,1,3> StateR1_3;

typedef State<Rn,double,2,1> StateR2_1;
typedef State<Rn,double,2,2> StateR2_2;
typedef State<Rn,double,2,3> StateR2_3;

typedef State<Rn,double,3,1> StateR3_1;
typedef State<Rn,double,3,2> StateR3_2;
typedef State<Rn,double,3,3> StateR3_3;

typedef State<Rn,double,4,1> StateR4_1;
typedef State<Rn,double,4,2> StateR4_2;
typedef State<Rn,double,4,3> StateR4_3;

// Create Source Types
typedef SourceRN<StateR1_1,MeasurementTypes::RN_POS_VEL,TransformNULL> SourceR1_1;
typedef SourceRN<StateR1_2,MeasurementTypes::RN_POS_VEL,TransformNULL> SourceR1_2;
typedef SourceRN<StateR1_3,MeasurementTypes::RN_POS_VEL,TransformNULL> SourceR1_3;

typedef SourceRN<StateR2_1,MeasurementTypes::RN_POS_VEL,TransformNULL> SourceR2_1;
typedef SourceRN<StateR2_2,MeasurementTypes::RN_POS_VEL,TransformNULL> SourceR2_2;
typedef SourceRN<StateR2_3,MeasurementTypes::RN_POS_VEL,TransformNULL> SourceR2_3;

typedef SourceRN<StateR3_1,MeasurementTypes::RN_POS_VEL,TransformNULL> SourceR3_1;
typedef SourceRN<StateR3_2,MeasurementTypes::RN_POS_VEL,TransformNULL> SourceR3_2;
typedef SourceRN<StateR3_3,MeasurementTypes::RN_POS_VEL,TransformNULL> SourceR3_3;


// Create Source Container Types
typedef SourceContainer<SourceR1_1> SourceContainerR1_1;
typedef SourceContainer<SourceR1_2> SourceContainerR1_2;
typedef SourceContainer<SourceR1_3> SourceContainerR1_3;

typedef SourceContainer<SourceR2_1> SourceContainerR2_1;
typedef SourceContainer<SourceR2_2> SourceContainerR2_2;
typedef SourceContainer<SourceR2_3> SourceContainerR2_3;

typedef SourceContainer<SourceR3_1> SourceContainerR3_1;
typedef SourceContainer<SourceR3_2> SourceContainerR3_2;
typedef SourceContainer<SourceR3_3> SourceContainerR3_3;


// Create Model Types

typedef ModelRN<SourceContainerR1_1> ModelR1_1;
typedef ModelRN<SourceContainerR1_2> ModelR1_2;
typedef ModelRN<SourceContainerR1_3> ModelR1_3;

typedef ModelRN<SourceContainerR2_1> ModelR2_1;
typedef ModelRN<SourceContainerR2_2> ModelR2_2;
typedef ModelRN<SourceContainerR2_3> ModelR2_3;

typedef ModelRN<SourceContainerR3_1> ModelR3_1;
typedef ModelRN<SourceContainerR3_2> ModelR3_2;
typedef ModelRN<SourceContainerR3_3> ModelR3_3;



using MyTypes = ::testing::Types<ModelR1_1,ModelR1_2,ModelR1_3,ModelR2_1,ModelR2_2,ModelR2_3,ModelR3_1,ModelR3_2,ModelR3_3>;
// using MyTypes = ::testing::Types<ModelR2_1>;


template<class Model>
class ModelTest : public ::testing::Test {

public: 



typedef typename Model::State State;
typedef typename Model::DataType DataType;
typedef typename Model::MatModelCov MatModelCov;
typedef typename Model::TransformDataType TransformDataType;

static constexpr unsigned int g_dim_ = State::Group::dim_;                      /**< The dimension of the pose of the state, i.e. the dimension of the group portion of the state. */
static constexpr unsigned int num_tangent_spaces_ = State::NumTangentSpaces;    /**< The number of tangent spaces in the group. Ex 1 tangent space is velocity and 2 is velocity and acceleration. */
static constexpr unsigned int u_dim_ = g_dim_*num_tangent_spaces_;              /**< The dimension of the total tangent space. */
static constexpr unsigned int cov_dim_ = g_dim_+u_dim_;                         /**< The dimension of the error covariance. */
// typedef Eigen::Matrix<DataType,cov_dim_,cov_dim_> Mat;                          /**< The object type of the error covariance, Jacobians, and others. */
typedef Eigen::Matrix<DataType,cov_dim_,1> VecCov;                                  /**< The object type of the error state. */

static constexpr unsigned int meas_pose_dim = Model::SourceContainer::Source0::meas_pose_dim_;
static constexpr unsigned int total_meas_dim_ = Model::SourceContainer::Source0::Base::total_meas_dim_;
static constexpr unsigned int state_dim = Model::State::dim_;
static constexpr unsigned int cov_dim = Model::cov_dim_;

typedef typename Model::SourceContainer::Source0::MatMeasCov   MeasSpaceCovMat;
typedef typename Model::Base::Measurement Measurement;
// static constexpr unsigned int a_vel_dim = Model::Model::cov_dim_ - Model::Model::State::g_type_::dim_-1;
// static constexpr unsigned int t_vel_dim = state_dim - Model::Model::State::g_type_::dim_ - a_vel_dim;

protected:

void SetUp() override {

// std::cerr << "setup" << std::endl;
// Setup the parameters
double meas_cov_scale = 0.001;
double system_cov_scale = 0.1;
meas_cov1.setIdentity();
meas_cov1 *= meas_cov_scale;
meas_std1.setIdentity();
meas_std1 *= sqrt(meas_cov_scale);
// source_params1.meas_cov_fixed_ = false;
source_params1.meas_cov_ = meas_cov1;
source_params1.type_ = MeasurementTypes::RN_POS_VEL;
source_params1.source_index_ = 0;
source_params1.spacial_density_of_false_meas_ = 0.8;
source_params1.probability_of_detection_ = 0.8;
source_params1.gate_probability_ = 0.8;



// std::cerr << "here0" << std::endl;
source_container.AddSource(source_params1);



params.process_noise_covariance_ = Eigen::Matrix<double,cov_dim,cov_dim>::Identity()*system_cov_scale;

// std::cerr << "here 1" << std::endl;

// Generate Trajectory and measurements
Measurement m1;
std::vector<std::vector<Measurement>> gen_meas(1);
std::vector<Measurement> meas1(num_meas);
std::vector<Measurement> meas2(num_meas);

state = Model::GetRandomState();
state.g_.data_*=10;
state.u_.data_*=5;



state0 = state;                           // Initial state
states.push_back(state);  
bool transform_state = false;
typename Model::SourceContainer::Source0::TransformDataType EmptyMat;              

// std::cout << "state 0: " << std::endl;
// std::cout << "state pose: " << std::endl << state0.g_.data_ << std::endl;
// std::cout << "state twist: " << std::endl << state0.u_.data_ << std::endl;

for (unsigned int ii = 0; ii < num_iters; ++ii) {

    // Propagate state
    Model::PropagateState(state,dt);

    states.push_back(state); 

    // Get measurements
    for (unsigned int jj=0; jj < num_meas; ++jj) {

        m1 = source_container.GenerateRandomMeasurement(source_params1.source_index_, meas_std1, state, transform_state, EmptyMat);
        m1.weight = 1.0/(num_meas+1);
        m1.time_stamp = ii*dt;
        meas1[jj] = m1;
    }
    gen_meas[0] = meas1;
    new_meas.push_back(gen_meas);

}

statef = state;

// std::cout << "state propagated: " << std::endl;
// std::cout << "state pose: " << std::endl << state.g_.data_ << std::endl;
// std::cout << "state twist: " << std::endl << state.u_.data_ << std::endl;


// Construct Jacobians


A.setZero();
A.block(0,g_dim_,u_dim_,u_dim_).setIdentity();

F = (A*dt).exp();




G.setIdentity();

}




SourceParameters source_params1;
MeasSpaceCovMat meas_cov1;
MeasSpaceCovMat meas_std1;
Parameters params;

State state;
State statef;
State state0;
std::vector<State> states;

typename Model::SourceContainer source_container;

std::vector<std::vector<std::vector<Measurement>>> new_meas; // time, source, measurements

int num_meas = 2;
int num_iters = 20;
double dt = 0.1;

Model track;
MatModelCov F;
MatModelCov G;

MatModelCov A;




};

TYPED_TEST_SUITE(ModelTest, MyTypes);

TYPED_TEST(ModelTest, JacobiansPropagateOMinusOPlus) {

    typedef ModelTest<TypeParam> TestType;
    typedef typename TestType::Measurement Measurement;

    // Test init function and set parameters function
    this->track.Init(this->params);
    ASSERT_EQ(this->track.Q_, this->params.process_noise_covariance_);
    ASSERT_EQ(this->track.err_cov_, TypeParam::MatModelCov::Identity());
    ASSERT_EQ(this->track.newest_measurement_time_stamp, 0);
    ASSERT_EQ(this->track.model_likelihood_, 0.5);
    ASSERT_EQ(this->track.label_, -1);
    ASSERT_EQ(this->track.innov_cov_set_.size(), 1);
    ASSERT_EQ(this->track.innovation_covariances_.size(), 1);
    ASSERT_EQ(this->track.new_assoc_meas_.size(), 1);
    ASSERT_EQ(this->track.model_likelihood_update_info_.size(), 1);


    // Test the system function Jacobians
    // std::cout << "F: " << std::endl << this->F << std::endl;
    // std::cout << "F2: " << std::endl << this->track.GetLinTransFuncMatState(this->state0,this->dt) << std::endl;
    ASSERT_LE( ( this->track.GetLinTransFuncMatState(this->state0,this->dt) - this->F).norm(), 1e-12);
    ASSERT_LE( ( this->track.GetLinTransFuncMatNoise(this->state0,this->dt) - this->G).norm(), 1e-12);

    // Test Propagate state

    typename TestType::VecCov state;
    state.block(0,0,TypeParam::g_dim_,1) = this->state0.g_.data_;
    state.block(TypeParam::g_dim_,0,TypeParam::u_dim_,1) = this->state0.u_.data_;
    state = (this->A*this->dt*this->num_iters).exp()*state;
    
    ASSERT_LT( (this->statef.g_.data_- state.block(0,0,TypeParam::g_dim_,1)).norm(), 1e-10);
    ASSERT_LT( (this->statef.u_.data_- state.block(TypeParam::g_dim_,0,TypeParam::u_dim_,1)).norm(), 1e-10);




    // Test the Propagate Model Function
    std::fill(this->track.innov_cov_set_.begin(), this->track.innov_cov_set_.end(), true); // make sure that this is reset to false.
    this->track.state_ = this->state0;
    Eigen::MatrixXd err_cov0 = this->track.err_cov_;
    this->track.PropagateModel(this->dt);
    for (auto val : this->track.innov_cov_set_) {
        ASSERT_FALSE(val);
    }

    // ASSERT_EQ(this->track.missed_detection_time_, this->dt);
    ASSERT_LE( (this->track.F_- this->F).norm(), 1e-12  );
    ASSERT_LE( (this->track.G_- this->G).norm(), 1e-12  );
    ASSERT_LE( (this->track.state_.g_.data_ - this->states[1].g_.data_).norm(), 1e-12  );
    ASSERT_LE( (this->track.state_.u_.data_ - this->states[1].u_.data_).norm(), 1e-12  );

    ASSERT_LE( (this->track.err_cov_ - this->F*err_cov0*this->F.transpose() - this->G*this->track.Q_*this->dt*this->G.transpose()).norm(), 1e-12  );



    // Test the innovation covariance function.
    ////////////////////////////////////////////////
    bool transform_state = false;
    typename TestType::TransformDataType EmptyMat;
    Eigen::MatrixXd H0, V0, S0, S0_e; // Observation function jacobians w.r.t. the state and noise for both sources and innovation covariances.
    H0 = this->source_container.GetLinObsMatState(0, this->track.state_, transform_state, EmptyMat);
    V0 = this->source_container.GetLinObsMatSensorNoise(0, this->track.state_, transform_state, EmptyMat);
    S0 = H0*this->track.err_cov_*H0.transpose() + V0*this->source_container.GetParams(0).meas_cov_*V0.transpose();
    S0_e = this->track.GetInnovationCovariance(this->source_container,0, transform_state, EmptyMat);
    ASSERT_EQ(this->track.innov_cov_set_.size(),1);
    ASSERT_EQ(this->track.innovation_covariances_.size(),1);
    ASSERT_TRUE(this->track.innov_cov_set_[0]);  // Make sure the innovation covariance for the first source was set and none other.
    ASSERT_GT(this->track.innovation_covariances_[0].rows(),0);
    ASSERT_LE( (S0- S0_e).norm(), 1e-12  );

    this->track.UpdateModel(this->source_container, this->params); // Make sure that the update model resets the innov_cov_set to false
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

        this->track.UpdateModel(this->source_container, this->params);


        

    }
    ASSERT_DOUBLE_EQ(this->track.newest_measurement_time_stamp,this->num_iters*this->dt - this->dt);

    // Test covariance
    ASSERT_LE(  ((this->track.err_cov_ + this->track.err_cov_.transpose())/2.0 - this->track.err_cov_ ).norm(), 1e-12); // symmetric
    Eigen::VectorXcd eigen_values = this->track.err_cov_.eigenvalues();
    for (int ii =0; ii < eigen_values.rows(); ++ii){                       // positive definite
        ASSERT_GT( std::real(eigen_values(ii)), 0);
    }


    EXPECT_LE( (this->track.err_cov_.inverse().sqrt()*this->track.state_.OMinus(this->states.back())).norm(),  5) << "This result might not always be true since updating the track is not deterministic";

    // Verify the Consensus Set
    ASSERT_EQ(this->track.cs_.Size(), this->num_iters);

    ConsensusSet<Measurement> set;
    set.consensus_set_.begin();

    for( typename std::list<std::vector<Measurement>>::iterator it = this->track.cs_.consensus_set_.begin(); it!= this->track.cs_.consensus_set_.end(); ++it) {
        ASSERT_EQ( (*it).size(), 2);
    }

    this->track.PruneConsensusSet(0);
    ASSERT_EQ(this->track.cs_.Size(), this->num_iters-1);

}


} // rransac