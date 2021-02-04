#include <gtest/gtest.h>
#include "common/sources/source_SEN_pose_twist.h"



namespace rransac {

using namespace lie_groups;


using MyTypes = ::testing::Types<SE2_se2,SE3_se3>;

template<typename T>
class SEN_POSE_TWIST_Test :public testing::Test {
    public:
    typedef T type;
};

TYPED_TEST_SUITE(SEN_POSE_TWIST_Test, MyTypes);

TYPED_TEST(SEN_POSE_TWIST_Test, INIT) {

typedef Eigen::Matrix<double,6,6> Mat6d;
typedef Eigen::Matrix<double,12,12> Mat12d;
typedef SourceSENPoseTwist<TypeParam> Source;
typedef Eigen::Matrix<double,Source::meas_space_dim_,Source::meas_space_dim_> MatMeas1;
typedef Eigen::Matrix<double,Source::meas_space_dim_*2,Source::meas_space_dim_*2> MatMeas2;

SourceParameters params;
params.spacial_density_of_false_meas_ = 0.1;
params.gate_probability_ = 0.8;
params.probability_of_detection_ = 0.9;


Source source;

// Valid state type
if (typeid(TypeParam).name() == typeid(SE2_se2).name() || typeid(TypeParam).name() == typeid(SE3_se3).name()) {

// Valid measurement types
params.type_ = MeasurementTypes::SEN_POSE;
params.meas_cov_ = MatMeas1::Identity();
ASSERT_NO_THROW(source.Init(params));
params.type_ = MeasurementTypes::SEN_POSE_TWIST;
params.meas_cov_ = MatMeas2::Identity();
ASSERT_NO_THROW(source.Init(params));

// Invalid measurements
params.type_ = MeasurementTypes::RN_POS;
params.meas_cov_ = MatMeas1::Identity();
ASSERT_ANY_THROW(source.Init(params));
params.type_ = MeasurementTypes::RN_POS_VEL;
params.meas_cov_ = MatMeas2::Identity();
ASSERT_ANY_THROW(source.Init(params));
params.type_ = MeasurementTypes::SEN_POS;
params.meas_cov_ = MatMeas1::Identity();
ASSERT_ANY_THROW(source.Init(params));
params.type_ = MeasurementTypes::SEN_POS_VEL;
params.meas_cov_ = MatMeas2::Identity();
ASSERT_ANY_THROW(source.Init(params));


}

// Invalid state type
else {

    // Valid measurement types
    params.type_ = MeasurementTypes::SEN_POSE;
    params.meas_cov_ = MatMeas1::Identity();
    ASSERT_ANY_THROW(source.Init(params));

}


}

///////////////////////////////////////
TYPED_TEST(SEN_POSE_TWIST_Test, Funcs) {

typedef Eigen::Matrix<double,6,6> Mat6d;
typedef Eigen::Matrix<double,12,12> Mat12d;
typedef SourceSENPoseTwist<TypeParam> Source;
typedef Eigen::Matrix<double,Source::meas_space_dim_,Source::meas_space_dim_> MatMeas1;
typedef Eigen::Matrix<double,Source::meas_space_dim_*2,Source::meas_space_dim_*2> MatMeas2;

typedef TypeParam S;
typedef Eigen::Matrix<double,TypeParam::Group::dim_,1> Mat_p;
typedef Eigen::Matrix<double,TypeParam::Group::dim_*2,1> Mat_pt;
typedef Eigen::Matrix<double,TypeParam::Group::dim_,  TypeParam::Group::dim_> Mat_p_c;
typedef Eigen::Matrix<double,TypeParam::Group::dim_*2,TypeParam::Group::dim_*2> Mat_pt_c;


SourceParameters params;
params.spacial_density_of_false_meas_ = 0.1;
params.gate_probability_ = 0.8;
params.probability_of_detection_ = 0.9;

Source source;
TypeParam state = TypeParam::Random();


// Construct the expected measurement
Meas<double> m;
m.pose = state.g_.data_;
m.twist = state.u_.data_;

// Construct the Jacobians
Eigen::MatrixXd V_pose = Eigen::Matrix<double,S::Group::dim_,S::Group::dim_>::Identity();
Eigen::MatrixXd V_pose_twist = Eigen::Matrix<double,S::dim_,S::dim_>::Identity();

Eigen::MatrixXd H_pose = Eigen::Matrix<double, S::Group::dim_, S::dim_>::Zero();
H_pose.block(0,0,S::Group::dim_,S::Group::dim_).setIdentity();

Eigen::MatrixXd H_pose_twist = Eigen::Matrix<double,S::dim_,S::dim_>::Identity();

// std::cout << V_pose << std::endl << std::endl;
// std::cout << V_pose_twist << std::endl << std::endl;
// std::cout << H_pose << std::endl << std::endl;
// std::cout << H_pose_twist << std::endl << std::endl;


// Tests
params.type_ = MeasurementTypes::SEN_POSE;
params.meas_cov_ = MatMeas1::Identity();
ASSERT_NO_THROW(source.Init(params));

ASSERT_EQ(source.GetLinObsMatState(state),H_pose);
ASSERT_EQ(source.GetLinObsMatState(state,params.type_),H_pose);
ASSERT_EQ(source.GetLinObsMatSensorNoise(state),V_pose);
ASSERT_EQ(source.GetLinObsMatSensorNoise(state,params.type_),V_pose);
ASSERT_EQ(source.GetEstMeas(state).pose,m.pose);
ASSERT_EQ(source.GetEstMeas(state,params.type_).pose,m.pose);

params.type_ = MeasurementTypes::SEN_POSE_TWIST;
m.type = SEN_POSE_TWIST;
params.meas_cov_ = MatMeas2::Identity();
ASSERT_NO_THROW(source.Init(params));

ASSERT_EQ(source.GetLinObsMatState(state),H_pose_twist);
ASSERT_EQ(source.GetLinObsMatState(state,params.type_),H_pose_twist);
ASSERT_EQ(source.GetLinObsMatSensorNoise(state),V_pose_twist);
ASSERT_EQ(source.GetLinObsMatSensorNoise(state,params.type_),V_pose_twist);
ASSERT_EQ(source.GetEstMeas(state).pose,m.pose);
ASSERT_EQ(source.GetEstMeas(state,params.type_).pose,m.pose);
ASSERT_EQ(source.GetEstMeas(state).twist,m.twist);
ASSERT_EQ(source.GetEstMeas(state,params.type_).twist,m.twist);





// Test OMinus
SourceParameters params1, params2;

params1.spacial_density_of_false_meas_ = 0.1;
params1.meas_cov_ = MatMeas1::Identity();
params1.gate_probability_ = 0.8;
params1.probability_of_detection_ = 0.9;

params2.spacial_density_of_false_meas_ = 0.1;
params2.meas_cov_ = MatMeas2::Identity();
params2.gate_probability_ = 0.8;
params2.probability_of_detection_ = 0.9;


params1.type_ = MeasurementTypes::SEN_POSE;
params2.type_ = MeasurementTypes::SEN_POSE_TWIST;

Source source1, source2;
source1.Init(params1);
source2.Init(params2);

TypeParam state_tmp = TypeParam::Random();

Meas<double> m3, m4;
m3.pose = state_tmp.g_.data_;
m3.twist = Eigen::Matrix<double,TypeParam::Group::dim_,1>::Random();
state_tmp = TypeParam::Random();
m4.pose = state_tmp.g_.data_;
m4.twist = Eigen::Matrix<double,TypeParam::Group::dim_,1>::Random();

Eigen::Matrix<double,TypeParam::Group::dim_*2,1> error2;
error2.block(0,0,TypeParam::Group::dim_,1) = TypeParam::Group::OMinus(m3.pose,m4.pose);
error2.block(TypeParam::Group::dim_,0,TypeParam::Group::dim_,1) = m3.twist - m4.twist;

m3.type = MeasurementTypes::SEN_POSE;
m4.type = MeasurementTypes::SEN_POSE;
ASSERT_LE( (source1.OMinus(m3,m4) - TypeParam::Group::OMinus(m3.pose,m4.pose)).norm(), 1e-8  );
m3.type = MeasurementTypes::SEN_POSE_TWIST;
m4.type = MeasurementTypes::SEN_POSE_TWIST;
ASSERT_LE( (source2.OMinus(m3,m4) - error2).norm(), 1e-8) ;



// Test Random Measurements
const int num_rand = 10000;
double std_scalar = 0.1;
Mat_p_c std1 = Mat_p_c::Identity()*std_scalar;
Mat_pt_c std2 = Mat_pt_c::Identity()*std_scalar;
std::vector<Meas<double>> rand_meas1(num_rand);
std::vector<Meas<double>> rand_meas2(num_rand);
std::vector<Mat_p> error_1(num_rand);
std::vector<Mat_pt> error_2(num_rand);
Mat_p error_mean1;
Mat_pt error_mean2;
error_mean1.setZero();
error_mean2.setZero();

// std::cout << "std1: " << std::endl << std1 << std::endl;

// Generate the random measurement, calculate the error between the state and measurement, and get the mean of the error
for (unsigned long int ii = 0; ii < num_rand; ++ii) {

    rand_meas1[ii].type = MeasurementTypes::SEN_POSE;
    rand_meas1[ii] = source1.GenerateRandomMeasurement(state,std1);
    rand_meas2[ii].type = MeasurementTypes::SEN_POSE_TWIST;
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
Mat_pt_c cov2;
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