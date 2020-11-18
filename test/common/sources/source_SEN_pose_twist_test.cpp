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

SourceParameters params;
SourceSENPoseTwist<TypeParam> source;

// Valid state type
if (typeid(TypeParam).name() == typeid(SE2_se2).name() || typeid(TypeParam).name() == typeid(SE3_se3).name()) {

// Valid measurement types
params.type_ = MeasurementTypes::SEN_POSE;
ASSERT_NO_THROW(source.Init(params));
params.type_ = MeasurementTypes::SEN_POSE_TWIST;
ASSERT_NO_THROW(source.Init(params));

// Invalid measurements
params.type_ = MeasurementTypes::RN_POS;
ASSERT_ANY_THROW(source.Init(params));
params.type_ = MeasurementTypes::RN_POS_VEL;
ASSERT_ANY_THROW(source.Init(params));
params.type_ = MeasurementTypes::SEN_POS;
ASSERT_ANY_THROW(source.Init(params));
params.type_ = MeasurementTypes::SEN_POS_VEL;
ASSERT_ANY_THROW(source.Init(params));


}

// Invalid state type
else {

    // Valid measurement types
    params.type_ = MeasurementTypes::SEN_POSE;
    ASSERT_ANY_THROW(source.Init(params));

}


}

///////////////////////////////////////
TYPED_TEST(SEN_POSE_TWIST_Test, Funcs) {

typedef TypeParam S;

SourceParameters params;
SourceSENPoseTwist<TypeParam> source;
TypeParam state = TypeParam::Random();

// Construct the expected measurement
Meas m;
m.pose = state.g_.data_;
m.twist = state.u_.data_;

// Construct the Jacobians
Eigen::MatrixXd V_pose = Eigen::Matrix<double,S::g_type_::dim_,S::g_type_::dim_>::Identity();
Eigen::MatrixXd V_pose_twist = Eigen::Matrix<double,S::dim_,S::dim_>::Identity();

Eigen::MatrixXd H_pose = Eigen::Matrix<double, S::g_type_::dim_, S::dim_>::Zero();
H_pose.block(0,0,S::g_type_::dim_,S::g_type_::dim_).setIdentity();

Eigen::MatrixXd H_pose_twist = Eigen::Matrix<double,S::dim_,S::dim_>::Identity();

// std::cout << V_pose << std::endl << std::endl;
// std::cout << V_pose_twist << std::endl << std::endl;
// std::cout << H_pose << std::endl << std::endl;
// std::cout << H_pose_twist << std::endl << std::endl;


// Tests
params.type_ = MeasurementTypes::SEN_POSE;
ASSERT_NO_THROW(source.Init(params));

ASSERT_EQ(source.GetLinObsMatState(state),H_pose);
ASSERT_EQ(source.GetLinObsMatSensorNoise(state),V_pose);
ASSERT_EQ(source.GetEstMeas(state).pose,m.pose);
ASSERT_EQ(source.GetEstMeas(state).twist,m.twist);

params.type_ = MeasurementTypes::SEN_POSE_TWIST;
ASSERT_NO_THROW(source.Init(params));

ASSERT_EQ(source.GetLinObsMatState(state),H_pose_twist);
ASSERT_EQ(source.GetLinObsMatSensorNoise(state),V_pose_twist);
ASSERT_EQ(source.GetEstMeas(state).pose,m.pose);
ASSERT_EQ(source.GetEstMeas(state).twist,m.twist);


// Test to Euclidean function by undoing it and testing the result against the original pose.
m.pose_euclidean = source.ToEuclidean(m);

Eigen::Matrix<double, TypeParam::g_type_::size1_, TypeParam::g_type_::size2_> pose;
pose.setZero();
pose(TypeParam::g_type_::size1_-1,TypeParam::g_type_::size2_-1) = 1;
Eigen::Matrix<double,  TypeParam::g_type_::dim_pos_,1> t = m.pose_euclidean.block(0,0,TypeParam::g_type_::dim_pos_,1);
Eigen::Matrix<double,  TypeParam::g_type_::dim_pos_,TypeParam::g_type_::dim_pos_> R;
Eigen::Matrix<double,  TypeParam::g_type_::dim_pos_,TypeParam::g_type_::dim_pos_> I = Eigen::Matrix<double,  TypeParam::g_type_::dim_pos_,TypeParam::g_type_::dim_pos_>::Identity();
R = TypeParam::g_type_::rot_algebra::Wedge(m.pose_euclidean.block(TypeParam::g_type_::dim_pos_,0,TypeParam::g_type_::dim_rot_,1));
R = (I+R/2.0)*(I-R/2.0).inverse();



pose.block(0,0,TypeParam::g_type_::dim_pos_,TypeParam::g_type_::dim_pos_) = R;
pose.block(0,TypeParam::g_type_::dim_pos_,TypeParam::g_type_::dim_pos_,1) = t;

ASSERT_LE( (m.pose-pose).norm() , 1e-8);


// Test OMinus
SourceParameters params1, params2;
params1.type_ = MeasurementTypes::SEN_POSE;
params2.type_ = MeasurementTypes::SEN_POSE_TWIST;

SourceSENPoseTwist<TypeParam> source1, source2;
source1.Init(params1);
source2.Init(params2);

TypeParam state_tmp = TypeParam::Random();

Meas m3, m4;
m3.pose = state_tmp.g_.data_;
m3.twist = Eigen::Matrix<double,TypeParam::g_type_::dim_,1>::Random();
state_tmp = TypeParam::Random();
m4.pose = state_tmp.g_.data_;
m4.twist = Eigen::Matrix<double,TypeParam::g_type_::dim_,1>::Random();

Eigen::Matrix<double,TypeParam::g_type_::dim_*2,1> error2;
error2.block(0,0,TypeParam::g_type_::dim_,1) = TypeParam::g_type_::OMinus(m3.pose,m4.pose);
error2.block(TypeParam::g_type_::dim_,0,TypeParam::g_type_::dim_,1) = m3.twist - m4.twist;

ASSERT_LE( (source1.OMinus(m3,m4) - TypeParam::g_type_::OMinus(m3.pose,m4.pose)).norm(), 1e-8  );
ASSERT_LE( (source2.OMinus(m3,m4) - error2).norm(), 1e-8) ;

}




} // namespace rransac