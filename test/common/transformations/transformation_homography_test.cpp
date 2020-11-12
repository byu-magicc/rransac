#include <gtest/gtest.h>
#include "common/transformations/trans_homography.h"
#include <chrono> 


namespace rransac
{

using namespace lie_groups;  

// template<class S>
// class TransHomograpy : public TransformHomography<S> {

// public:



// Eigen::Matrix2d GetH1(){return this->H1_;}
// Eigen::Matrix<double,2,1> Geth2(){return this->h2_;}
// Eigen::Matrix<double,1,2> Geth3(){return this->h3_T_;}
// Eigen::Matrix<double,1,1> Geth4(){return this->h4_;}


// };

template<class S>
class HomographyTest : public ::testing::Test {

protected:



};

using MyTypes = ::testing::Types<R2_r2, SE2_se2>;
TYPED_TEST_SUITE(HomographyTest, MyTypes);


TYPED_TEST(HomographyTest, SetData) {

Eigen::Matrix3d H;
Eigen::Matrix2d H1;
Eigen::Matrix<double,2,1> h2;
Eigen::Matrix<double,1,2> h3_T;
Eigen::Matrix<double,1,1> h4;
H << 1,2,3,4,5,6,7,8,9;
H1 << 1,2,4,5;
h2 <<3,6;
h3_T << 7,8;
h4 << 9;

TransformHomography<TypeParam> trans;


trans.SetData(H);

ASSERT_EQ(trans.GetData(), H);
ASSERT_EQ(trans.H1_, H1);
ASSERT_EQ(trans.h2_, h2);
ASSERT_EQ(trans.h3_T_, h3_T);
ASSERT_EQ(trans.h4_, h4);

}


// Measurement transformation

TEST(TransformHomographyTest, MeasurementTransformationR2) {


TransformHomography<R2_r2> trans;

// Construct Homography.
SO3 R = SO3::Random();
Eigen::Matrix3d H{R.data_};

trans.SetData(H);

// Construct pixel measurement
Meas m1,m2,m3;
m1.pose = Eigen::Matrix<double,2,1>::Random();
m1.twist = Eigen::Matrix<double,2,1>::Random();
m1.type = MeasurementTypes::RN_POS_VEL;
m2 = m1;

// Verify velocity transform. i.e. ConstructVelTransform function;
double dt = 1e-8;
Eigen::Matrix<double,2,1> d1,d2;
d1 << dt,0;
d2 << 0,dt;
Eigen::Matrix2d vel_trans_numerical;
vel_trans_numerical.block(0,0,2,1) = (trans.TransformPosition(m1.pose + d1) - trans.TransformPosition(m1.pose))/dt;
vel_trans_numerical.block(0,1,2,1) = (trans.TransformPosition(m1.pose + d2) - trans.TransformPosition(m1.pose))/dt;
Eigen::Matrix2d vel_trans_analytical = trans.ConstructVelTransform(m1.pose);
ASSERT_LT( (vel_trans_analytical-vel_trans_numerical).norm(),1e-8); // small error


// Verify a part of the covariance transform testing the function ConstructCovTrans12
Eigen::Matrix2d cov_12_numerical;
cov_12_numerical.block(0,0,2,1) = (trans.ConstructVelTransform(m1.pose + d1)*m1.twist - trans.ConstructVelTransform(m1.pose)*m1.twist)/dt;
cov_12_numerical.block(0,1,2,1) = (trans.ConstructVelTransform(m1.pose + d2)*m1.twist - trans.ConstructVelTransform(m1.pose)*m1.twist)/dt;
Eigen::Matrix2d cov_12_analytical = trans.ConstructCovTrans12(m1.pose,m1.twist);
ASSERT_LT( (cov_12_analytical-cov_12_numerical).norm(),1e-8); // small error



// std::cout << "num: " << std::endl << cov_12_numerical << std::endl;
// std::cout << "anal: " << std::endl << cov_12_analytical << std::endl;






// Calculate own measurement
Eigen::Matrix<double,3,1> pos1,pos2, vel, vel_trans;
pos1.block(0,0,2,1) = m1.pose;
pos1(2) = 1;
vel.block(0,0,2,1) = m1.twist;
vel(2) = 0;
pos2 = pos1 + vel;



// Transform them using the homography
pos1 = H*pos1;
pos1 = pos1/pos1(2);
pos2 = H*pos2;
pos2 = pos2/pos2(2);
vel_trans = pos2-pos1;
m3.pose = pos1.block(0,0,2,1);
m3.twist = vel_trans.block(0,0,2,1);

trans.TransformMeasurement(m2);

ASSERT_EQ(m2.pose,m3.pose);
ASSERT_LT( (m2.twist.normalized()-m3.twist.normalized()).norm(),1e-10); // small error


// Construct state and error covariance to test track transformation
R2_r2 state1 = R2_r2::Random();
R2_r2 state_trans; // the transformed state
Eigen::Matrix4d error_cov, error_cov_transformed;
error_cov.setIdentity();
Eigen::Matrix4d cov_trans;
cov_trans.setZero();
cov_trans.block(0,0,2,2) = trans.ConstructVelTransform(state1.g_.data_);
cov_trans.block(2,2,2,2) = trans.ConstructVelTransform(state1.g_.data_);
cov_trans.block(2,0,2,2) = trans.ConstructCovTrans12(state1.g_.data_,state1.u_.data_);

// Create transformed state and error covariance
state_trans.g_.data_ = trans.TransformPosition(state1.g_.data_);
state_trans.u_.data_ = trans.ConstructVelTransform(state1.g_.data_)*state1.u_.data_;
error_cov_transformed = cov_trans * error_cov *cov_trans.transpose();

// Transform state and covariance
trans.TransformTrack(state1,error_cov);

std::cout << "cov: " << std::endl << cov_trans << std::endl;
ASSERT_EQ(state1.g_.data_, state_trans.g_.data_);
ASSERT_EQ(state1.u_.data_, state_trans.u_.data_);
ASSERT_EQ(error_cov, error_cov_transformed);





}


} // namespace rransac
