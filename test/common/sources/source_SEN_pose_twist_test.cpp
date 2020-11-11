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

}




} // namespace rransac