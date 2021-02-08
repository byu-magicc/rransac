#include <gtest/gtest.h>
#include "rransac/common/sources/source_SEN_pos_vel.h"


namespace rransac
{
    
using namespace lie_groups;

// using MyTypes = ::testing::Types<R2_r2,R3_r3,SO2_so2,SO3_so3,SE2_se2,SE3_se3>;
using MyTypes = ::testing::Types<SE2_se2,SE3_se3>;

template<typename T>
class SEN_POS_VEL_Test :public testing::Test {
    public:
    typedef T type;

};

TYPED_TEST_SUITE(SEN_POS_VEL_Test, MyTypes);

TYPED_TEST(SEN_POS_VEL_Test, INIT) {

typedef Eigen::Matrix<double,TypeParam::Group::dim_pos_,1> Mat_p;
typedef Eigen::Matrix<double,TypeParam::Group::dim_pos_*2,1> Mat_pv;
typedef Eigen::Matrix<double,TypeParam::Group::dim_pos_,  TypeParam::Group::dim_pos_> Mat_p_c;
typedef Eigen::Matrix<double,TypeParam::Group::dim_pos_*2,TypeParam::Group::dim_pos_*2> Mat_pv_c;
typedef SourceSENPosVel<TypeParam> Source;

SourceParameters params;
params.spacial_density_of_false_meas_ = 0.1;
params.gate_probability_ = 0.8;
params.probability_of_detection_ = 0.9;

Source source;

// Valid state type
if (typeid(TypeParam).name() == typeid(SE2_se2).name() || typeid(TypeParam).name() == typeid(SE3_se3).name()) {

// Valid measurement types
params.type_ = MeasurementTypes::SEN_POS;
params.meas_cov_ = Eigen::Matrix<double,Source::meas_space_dim_,Source::meas_space_dim_>::Identity();
ASSERT_NO_THROW(source.Init(params));
params.meas_cov_ = Eigen::Matrix<double,Source::meas_space_dim_*2,Source::meas_space_dim_*2>::Identity();
params.type_ = MeasurementTypes::SEN_POS_VEL;
ASSERT_NO_THROW(source.Init(params));

// Invalid measurements
params.type_ = MeasurementTypes::RN_POS;
ASSERT_ANY_THROW(source.Init(params));
params.type_ = MeasurementTypes::RN_POS_VEL;
ASSERT_ANY_THROW(source.Init(params));
params.type_ = MeasurementTypes::SEN_POSE;
ASSERT_ANY_THROW(source.Init(params));
params.type_ = MeasurementTypes::SEN_POSE_TWIST;
ASSERT_ANY_THROW(source.Init(params));


}

// Invalid state type
else {

    // Valid measurement types
    params.type_ = MeasurementTypes::SEN_POS;
    ASSERT_ANY_THROW(source.Init(params));

}


}


///////////////////////////////////////
TYPED_TEST(SEN_POS_VEL_Test, Funcs) {

typedef TypeParam S;
typedef Eigen::Matrix<double,TypeParam::Group::dim_pos_,1> Mat_p;
typedef Eigen::Matrix<double,TypeParam::Group::dim_pos_*2,1> Mat_pv;
typedef Eigen::Matrix<double,TypeParam::Group::dim_pos_,  TypeParam::Group::dim_pos_> Mat_p_c;
typedef Eigen::Matrix<double,TypeParam::Group::dim_pos_*2,TypeParam::Group::dim_pos_*2> Mat_pv_c;
typedef SourceSENPosVel<TypeParam> Source;


SourceParameters params;
params.spacial_density_of_false_meas_ = 0.1;
params.gate_probability_ = 0.8;
params.probability_of_detection_ = 0.9;

Source source;

// Construct a state
TypeParam state = TypeParam::Random();
state.g_.R_.block(0,0,state.u_.p_.rows(),1) = state.u_.p_.normalized(); 

if(state.g_.dim_ == 3) {
    state.g_.R_.block(0,1,state.u_.p_.rows(),1) << - state.g_.data_(1,0), state.g_.data_(0,0);
} else {

    state.g_.R_.block(0,1,state.u_.p_.rows(),1) << -state.g_.data_(1,0), state.g_.data_(0,0), state.g_.data_(2,0);
    // std::cout << state.g_.R_.block(0,2,state.u_.p_.rows(),1)<< std::endl << std::endl;
    // std::cout << S::Group::RotAlgebra::Wedge(state.g_.data_.block(0,0,state.u_.p_.rows(),1)) << std::endl << std::endl;
    // std::cout << state.g_.data_.block(0,1,state.u_.p_.rows(),1)<< std::endl << std::endl;
    
    state.g_.R_.block(0,2,state.u_.p_.rows(),1) = S::Group::RotAlgebra::Wedge(state.g_.data_.block(0,0,state.u_.p_.rows(),1))*state.g_.data_.block(0,1,state.u_.p_.rows(),1);

}

double px = state.u_.p_.norm();
state.u_.p_.setZero();
state.u_.p_(0,0) = px;

// Construct the expected measurement
Meas<double> m;
m.pose = state.g_.t_;
m.twist = state.g_.R_*state.u_.p_;

// Construct the Jacobians
Eigen::Matrix<double,TypeParam::Group::dim_pos_,TypeParam::Group::dim_pos_> V_pos;
V_pos.setIdentity();
Eigen::MatrixXd V_pos_vel = Eigen::Matrix<double,S::Group::dim_pos_ + S::Algebra::dim_t_vel_,S::Group::dim_pos_+ S::Algebra::dim_t_vel_>::Identity();

Eigen::MatrixXd H_pos = Eigen::Matrix<double, S::Group::dim_pos_, SourceSENPosVel<TypeParam>::cov_dim_>::Zero();
H_pos.block(0,0,S::Group::dim_pos_, S::Group::dim_pos_) = state.g_.R_;

Eigen::MatrixXd H_pos_vel = Eigen::Matrix<double,S::Group::dim_pos_+S::Algebra::dim_t_vel_, SourceSENPosVel<TypeParam>::cov_dim_>::Zero();
H_pos_vel.block(0,0,S::Group::dim_pos_, S::Group::dim_pos_) = state.g_.R_;
H_pos_vel.block(S::Group::dim_pos_, S::Group::dim_, S::Algebra::dim_t_vel_,1) = state.g_.R_.block(0,0,S::Algebra::dim_t_vel_,1);

if ( S::Group::dim_pos_== 2) { 
H_pos_vel.block(S::Group::dim_pos_,S::Group::dim_pos_, S::Algebra::dim_t_vel_, S::Group::RotAlgebra::dim_) = state.g_.R_ * S::Group::RotAlgebra::Wedge(Eigen::Matrix<double,S::Group::RotAlgebra::dim_,1>::Ones()) * state.u_.p_;
} else {
    H_pos_vel.block(S::Group::dim_pos_,S::Group::dim_pos_, S::Algebra::dim_t_vel_, S::Group::RotAlgebra::dim_) = -state.g_.R_ * S::Group::RotAlgebra::Wedge(state.u_.p_.block(0,0,S::Group::RotAlgebra::dim_,1));
}

// std::cout << V_pos << std::endl << std::endl;
// std::cout << V_pos_vel << std::endl << std::endl;
// std::cout << H_pos << std::endl << std::endl;
// std::cout << H_pos_vel << std::endl << std::endl;
// std::cout << state.g_.R_ << std::endl << std::endl;
// std::cout << state.g_.R_.transpose()*state.g_.R_ << std::endl << std::endl;




// Tests
params.type_ = MeasurementTypes::SEN_POS;
params.meas_cov_ = Eigen::Matrix<double,Source::meas_space_dim_,Source::meas_space_dim_>::Identity();
ASSERT_NO_THROW(source.Init(params));

ASSERT_EQ(source.GetLinObsMatState(state),H_pos);
ASSERT_EQ(source.GetLinObsMatState(state,params.type_),H_pos);
ASSERT_EQ(source.GetLinObsMatSensorNoise(state),V_pos);
ASSERT_EQ(source.GetLinObsMatSensorNoise(state,params.type_),V_pos);
ASSERT_EQ(source.GetEstMeas(state).pose,m.pose);
ASSERT_EQ(source.GetEstMeas(state,params.type_).pose,m.pose);
// ASSERT_EQ(source.GetEstMeas(state).twist,m.twist);

params.type_ = MeasurementTypes::SEN_POS_VEL;
params.meas_cov_ = Eigen::Matrix<double,Source::meas_space_dim_*2,Source::meas_space_dim_*2>::Identity();
ASSERT_NO_THROW(source.Init(params));

ASSERT_EQ(source.GetLinObsMatState(state),H_pos_vel);
ASSERT_EQ(source.GetLinObsMatState(state,params.type_),H_pos_vel);
ASSERT_EQ(source.GetLinObsMatSensorNoise(state),V_pos_vel);
ASSERT_EQ(source.GetLinObsMatSensorNoise(state,params.type_),V_pos_vel);
ASSERT_EQ(source.GetEstMeas(state).pose,m.pose);
ASSERT_EQ(source.GetEstMeas(state,params.type_).pose,m.pose);
ASSERT_EQ(source.GetEstMeas(state).twist,m.twist);
ASSERT_EQ(source.GetEstMeas(state,params.type_).twist,m.twist);


// Test ToEuclidean


// Test OMinus
SourceParameters params1, params2;

params1.spacial_density_of_false_meas_ = 0.1;
params1.gate_probability_ = 0.8;
params1.probability_of_detection_ = 0.9;

params2.spacial_density_of_false_meas_ = 0.1;
params2.gate_probability_ = 0.8;
params2.probability_of_detection_ = 0.9;


params1.type_ = MeasurementTypes::SEN_POS;
params1.meas_cov_ = Eigen::Matrix<double,Source::meas_space_dim_,Source::meas_space_dim_>::Identity();
params2.type_ = MeasurementTypes::SEN_POS_VEL;
params2.meas_cov_ = Eigen::Matrix<double,Source::meas_space_dim_*2,Source::meas_space_dim_*2>::Identity();

Source source1, source2;
source1.Init(params1);
source2.Init(params2);

Meas<double> m3, m4;
m3.type = params1.type_;
m4.type = params2.type_;
m3.pose = Eigen::Matrix<double,TypeParam::Group::dim_pos_,1>::Random();
m3.twist = Eigen::Matrix<double,TypeParam::Group::dim_pos_,1>::Random();
m4.pose = Eigen::Matrix<double,TypeParam::Group::dim_pos_,1>::Random();
m4.twist = Eigen::Matrix<double,TypeParam::Group::dim_pos_,1>::Random();

Eigen::Matrix<double,TypeParam::Group::dim_pos_*2,1> error2;
error2.block(0,0,TypeParam::Group::dim_pos_,1) = m3.pose - m4.pose;
error2.block(TypeParam::Group::dim_pos_,0,TypeParam::Group::dim_pos_,1) = m3.twist - m4.twist;

m4.type = params1.type_;
ASSERT_EQ( source1.OMinus(m3,m4), m3.pose - m4.pose);
m3.type = params2.type_;
m4.type = params2.type_;
ASSERT_EQ( source2.OMinus(m3,m4), error2);


// Test Random Measurements
const int num_rand = 10000;
double std_scalar = 0.1;
Mat_p_c std1 = Mat_p_c::Identity()*std_scalar;
Mat_pv_c std2 = Mat_pv_c::Identity()*std_scalar;
std::vector<Meas<double>> rand_meas1(num_rand);
std::vector<Meas<double>> rand_meas2(num_rand);
std::vector<Mat_p> error_1(num_rand);
std::vector<Mat_pv> error_2(num_rand);
Mat_p error_mean1;
Mat_pv error_mean2;
error_mean1.setZero();
error_mean2.setZero();

// std::cout << "std1: " << std::endl << std1 << std::endl;

// Generate the random measurement, calculate the error between the state and measurement, and get the mean of the error
for (unsigned long int ii = 0; ii < num_rand; ++ii) {

    rand_meas1[ii].type = source1.params_.type_;
    rand_meas1[ii] = source1.GenerateRandomMeasurement(state,std1);
    rand_meas2[ii].type = source2.params_.type_;
    rand_meas2[ii] = source2.GenerateRandomMeasurement(state,std2);

    error_1[ii] = source1.OMinus(rand_meas1[ii], source1.GetEstMeas(state));
    error_2[ii] = source2.OMinus(rand_meas2[ii], source2.GetEstMeas(state));
    error_mean1+=error_1[ii];
    error_mean2+=error_2[ii];

}
error_mean1 /=num_rand;
error_mean2 /=num_rand;

// Calculate the covariance
Mat_p_c cov1;
Mat_pv_c cov2;
cov1.setZero();
cov2.setZero();
for (unsigned long int ii = 0; ii < num_rand; ++ii){

    cov1 += (error_mean1 - error_1[ii])*(error_mean1 - error_1[ii]).transpose();
    cov2 += (error_mean2 - error_2[ii])*(error_mean2 - error_2[ii]).transpose();

}
cov1 /= num_rand;
cov2 /= num_rand;


// std::cout << "error_mean1: " << std::endl << error_mean1 << std::endl;
// std::cout << "error_mean2: " << std::endl << error_mean2 << std::endl;
// std::cout << "cov1: " << std::endl << cov1 << std::endl;
// std::cout << "cov2: " << std::endl << cov2 << std::endl;


ASSERT_LE( (error_mean1).norm(), 0.1);
ASSERT_LE( (cov1 - std1*std1).norm(), 0.1);
ASSERT_LE( (error_mean2).norm(), 0.1);
ASSERT_LE( (cov2 - std2*std2).norm(), 0.1);



}




} // namespace rransac
