#include <gtest/gtest.h>
#include "common/sources/source_SEN_pos_vel.h"


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

typedef Eigen::Matrix<double,TypeParam::g_type_::dim_pos_,1> Mat_p;
typedef Eigen::Matrix<double,TypeParam::g_type_::dim_pos_*2,1> Mat_pv;
typedef Eigen::Matrix<double,TypeParam::g_type_::dim_pos_,  TypeParam::g_type_::dim_pos_> Mat_p_c;
typedef Eigen::Matrix<double,TypeParam::g_type_::dim_pos_*2,TypeParam::g_type_::dim_pos_*2> Mat_pv_c;

SourceParameters params;
SourceSENPosVel<TypeParam> source;

// Valid state type
if (typeid(TypeParam).name() == typeid(SE2_se2).name() || typeid(TypeParam).name() == typeid(SE3_se3).name()) {

// Valid measurement types
params.type_ = MeasurementTypes::SEN_POS;
ASSERT_NO_THROW(source.Init(params));
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
typedef Eigen::Matrix<double,TypeParam::g_type_::dim_pos_,1> Mat_p;
typedef Eigen::Matrix<double,TypeParam::g_type_::dim_pos_*2,1> Mat_pv;
typedef Eigen::Matrix<double,TypeParam::g_type_::dim_pos_,  TypeParam::g_type_::dim_pos_> Mat_p_c;
typedef Eigen::Matrix<double,TypeParam::g_type_::dim_pos_*2,TypeParam::g_type_::dim_pos_*2> Mat_pv_c;

SourceParameters params;
SourceSENPosVel<TypeParam> source;

// Construct a state
TypeParam state = TypeParam::Random();
state.g_.R_.block(0,0,state.u_.p_.rows(),1) = state.u_.p_.normalized(); 

if(state.g_.dim_ == 3) {
    state.g_.R_.block(0,1,state.u_.p_.rows(),1) << - state.g_.data_(1,0), state.g_.data_(0,0);
} else {

    state.g_.R_.block(0,1,state.u_.p_.rows(),1) << -state.g_.data_(1,0), state.g_.data_(0,0), state.g_.data_(2,0);
    // std::cout << state.g_.R_.block(0,2,state.u_.p_.rows(),1)<< std::endl << std::endl;
    // std::cout << S::g_type_::rot_algebra::Wedge(state.g_.data_.block(0,0,state.u_.p_.rows(),1)) << std::endl << std::endl;
    // std::cout << state.g_.data_.block(0,1,state.u_.p_.rows(),1)<< std::endl << std::endl;
    
    state.g_.R_.block(0,2,state.u_.p_.rows(),1) = S::g_type_::rot_algebra::Wedge(state.g_.data_.block(0,0,state.u_.p_.rows(),1))*state.g_.data_.block(0,1,state.u_.p_.rows(),1);

}

double px = state.u_.p_.norm();
state.u_.p_.setZero();
state.u_.p_(0,0) = px;

// Construct the expected measurement
Meas m;
m.pose = state.g_.t_;
m.twist = state.g_.R_*state.u_.p_;

// Construct the Jacobians
Eigen::Matrix<double,TypeParam::g_type_::dim_pos_,TypeParam::g_type_::dim_pos_> V_pos;
V_pos.setIdentity();
Eigen::MatrixXd V_pos_vel = Eigen::Matrix<double,S::g_type_::dim_pos_ + S::u_type_::dim_t_vel_,S::g_type_::dim_pos_+ S::u_type_::dim_t_vel_>::Identity();

Eigen::MatrixXd H_pos = Eigen::Matrix<double, S::g_type_::dim_pos_, SourceSENPosVel<TypeParam>::cov_dim_>::Zero();
H_pos.block(0,0,S::g_type_::dim_pos_, S::g_type_::dim_pos_) = state.g_.R_;

Eigen::MatrixXd H_pos_vel = Eigen::Matrix<double,S::g_type_::dim_pos_+S::u_type_::dim_t_vel_, SourceSENPosVel<TypeParam>::cov_dim_>::Zero();
H_pos_vel.block(0,0,S::g_type_::dim_pos_, S::g_type_::dim_pos_) = state.g_.R_;
H_pos_vel.block(S::g_type_::dim_pos_, S::g_type_::dim_, S::u_type_::dim_t_vel_,1) = state.g_.R_.block(0,0,S::u_type_::dim_t_vel_,1);

if ( S::g_type_::dim_pos_== 2) { 
H_pos_vel.block(S::g_type_::dim_pos_,S::g_type_::dim_pos_, S::u_type_::dim_t_vel_, S::g_type_::rot_algebra::dim_) = state.g_.R_ * S::g_type_::rot_algebra::Wedge(Eigen::Matrix<double,S::g_type_::rot_algebra::dim_,1>::Ones()) * state.u_.p_;
} else {
    H_pos_vel.block(S::g_type_::dim_pos_,S::g_type_::dim_pos_, S::u_type_::dim_t_vel_, S::g_type_::rot_algebra::dim_) = -state.g_.R_ * S::g_type_::rot_algebra::Wedge(state.u_.p_.block(0,0,S::g_type_::rot_algebra::dim_,1));
}

// std::cout << V_pos << std::endl << std::endl;
// std::cout << V_pos_vel << std::endl << std::endl;
// std::cout << H_pos << std::endl << std::endl;
// std::cout << H_pos_vel << std::endl << std::endl;
// std::cout << state.g_.R_ << std::endl << std::endl;
// std::cout << state.g_.R_.transpose()*state.g_.R_ << std::endl << std::endl;




// Tests
params.type_ = MeasurementTypes::SEN_POS;
ASSERT_NO_THROW(source.Init(params));

ASSERT_EQ(source.GetLinObsMatState(state),H_pos);
ASSERT_EQ(source.GetLinObsMatSensorNoise(state),V_pos);
ASSERT_EQ(source.GetEstMeas(state).pose,m.pose);
// ASSERT_EQ(source.GetEstMeas(state).twist,m.twist);

params.type_ = MeasurementTypes::SEN_POS_VEL;
ASSERT_NO_THROW(source.Init(params));

ASSERT_EQ(source.GetLinObsMatState(state),H_pos_vel);
ASSERT_EQ(source.GetLinObsMatSensorNoise(state),V_pos_vel);
ASSERT_EQ(source.GetEstMeas(state).pose,m.pose);
ASSERT_EQ(source.GetEstMeas(state).twist,m.twist);


// Test ToEuclidean

m.pose_euclidean = source.ToEuclidean(m);
ASSERT_EQ(m.pose,m.pose_euclidean);

// Test OMinus
SourceParameters params1, params2;
params1.type_ = MeasurementTypes::SEN_POS;
params2.type_ = MeasurementTypes::SEN_POS_VEL;

SourceSENPosVel<TypeParam> source1, source2;
source1.Init(params1);
source2.Init(params2);

Meas m3, m4;
m3.pose = Eigen::Matrix<double,TypeParam::g_type_::dim_pos_,1>::Random();
m3.twist = Eigen::Matrix<double,TypeParam::g_type_::dim_pos_,1>::Random();
m4.pose = Eigen::Matrix<double,TypeParam::g_type_::dim_pos_,1>::Random();
m4.twist = Eigen::Matrix<double,TypeParam::g_type_::dim_pos_,1>::Random();

Eigen::Matrix<double,TypeParam::g_type_::dim_pos_*2,1> error2;
error2.block(0,0,TypeParam::g_type_::dim_pos_,1) = m3.pose - m4.pose;
error2.block(TypeParam::g_type_::dim_pos_,0,TypeParam::g_type_::dim_pos_,1) = m3.twist - m4.twist;

ASSERT_EQ( source1.OMinus(m3,m4), m3.pose - m4.pose);
ASSERT_EQ( source2.OMinus(m3,m4), error2);


// Test Random Measurements
const int num_rand = 10000;
double std_scalar = 0.1;
Mat_p_c std1 = Mat_p_c::Identity()*std_scalar;
Mat_pv_c std2 = Mat_pv_c::Identity()*std_scalar;
std::vector<Meas> rand_meas1(num_rand);
std::vector<Meas> rand_meas2(num_rand);
std::vector<Mat_p> error_1(num_rand);
std::vector<Mat_pv> error_2(num_rand);
Mat_p error_mean1;
Mat_pv error_mean2;
error_mean1.setZero();
error_mean2.setZero();

// std::cout << "std1: " << std::endl << std1 << std::endl;

// Generate the random measurement, calculate the error between the state and measurement, and get the mean of the error
for (unsigned long int ii = 0; ii < num_rand; ++ii) {

    rand_meas1[ii] = source1.GenerateRandomMeasurement(state,std1);
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
