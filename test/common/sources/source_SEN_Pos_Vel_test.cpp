#include <gtest/gtest.h>
#include "common/sources/source_SEN_Pos_Vel.h"


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

SourceParameters params;
SourceSENPosVel<TypeParam> source;
TypeParam state = TypeParam::Random();

// Construct the expected measurement
Meas m;
m.pose = state.g_.t_;
m.twist = state.u_.p_;

// Construct the Jacobians
Eigen::Matrix<double,TypeParam::g_type_::dim_pos_,TypeParam::g_type_::dim_pos_> V_pos;
V_pos.setIdentity();
Eigen::MatrixXd V_pos_vel = Eigen::Matrix<double,S::g_type_::dim_pos_ + S::u_type_::dim_t_vel_,S::g_type_::dim_pos_+ S::u_type_::dim_t_vel_>::Identity();

Eigen::MatrixXd H_pos = Eigen::Matrix<double, S::g_type_::dim_pos_, S::dim_>::Zero();
H_pos.block(0,0,S::g_type_::dim_pos_, S::g_type_::dim_pos_) = state.g_.R_;

Eigen::MatrixXd H_pos_vel = Eigen::Matrix<double,S::g_type_::dim_pos_+S::u_type_::dim_t_vel_, S::dim_>::Zero();
H_pos_vel.block(S::g_type_::dim_pos_,S::g_type_::dim_,S::u_type_::dim_t_vel_,S::u_type_::dim_).setIdentity();
H_pos_vel.block(0,0,S::g_type_::dim_pos_, S::g_type_::dim_pos_) = state.g_.R_;

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
ASSERT_EQ(source.GetEstMeas(state).twist,m.twist);

params.type_ = MeasurementTypes::SEN_POS_VEL;
ASSERT_NO_THROW(source.Init(params));

ASSERT_EQ(source.GetLinObsMatState(state),H_pos_vel);
ASSERT_EQ(source.GetLinObsMatSensorNoise(state),V_pos_vel);
ASSERT_EQ(source.GetEstMeas(state).pose,m.pose);
ASSERT_EQ(source.GetEstMeas(state).twist,m.twist);

}




} // namespace rransac
