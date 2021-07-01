#include <gtest/gtest.h>
#include <chrono> 
#include <iostream>

#include "rransac/common/transformations/trans_homography.h"
#include "rransac/common/measurement/measurement_base.h"


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

using MyTypes = ::testing::Types<TransformHomography<R2_r2>,TransformHomography<SE2_se2>>;
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

TypeParam trans;

trans.Init();
trans.SetData(H);

ASSERT_EQ(trans.GetData(), H);
ASSERT_EQ(trans.H1_, H1);
ASSERT_EQ(trans.h2_, h2);
ASSERT_EQ(trans.h3_T_, h3_T);
ASSERT_EQ(trans.h4_, h4);

}


// Test the transformation for R2_r2

TEST(TransformHomographyTest, R2_r2_TransformationR2) {


TransformHomography<R2_r2> trans;

// Construct Homography.
SO3<double> R = SO3<double>::Random();
Eigen::Matrix3d H{R.data_};

trans.Init();
trans.SetData(H);

// Construct pixel measurement
Meas<double> m1,m2,m3;
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
Eigen::Matrix2d vel_trans_analytical = trans.ConstructTranslationalVelTransform(m1.pose);
ASSERT_LT( (vel_trans_analytical-vel_trans_numerical).norm(),1e-8); // small error


// Verify a part of the covariance transform testing the function ConstructCovTrans12
Eigen::Matrix2d cov_12_numerical;
cov_12_numerical.block(0,0,2,1) = (trans.ConstructTranslationalVelTransform(m1.pose + d1)*m1.twist - trans.ConstructTranslationalVelTransform(m1.pose)*m1.twist)/dt;
cov_12_numerical.block(0,1,2,1) = (trans.ConstructTranslationalVelTransform(m1.pose + d2)*m1.twist - trans.ConstructTranslationalVelTransform(m1.pose)*m1.twist)/dt;
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
error_cov.setRandom();
Eigen::Matrix4d cov_trans;
cov_trans.setZero();
cov_trans.block(0,0,2,2) = trans.ConstructTranslationalVelTransform(state1.g_.data_);
cov_trans.block(2,2,2,2) = trans.ConstructTranslationalVelTransform(state1.g_.data_);
cov_trans.block(2,0,2,2) = trans.ConstructCovTrans12(state1.g_.data_,state1.u_.data_);

// Create transformed state and error covariance
state_trans.g_.data_ = trans.TransformPosition(state1.g_.data_);
state_trans.u_.data_ = trans.ConstructTranslationalVelTransform(state1.g_.data_)*state1.u_.data_;
error_cov_transformed = cov_trans * error_cov *cov_trans.transpose();

// Transform state and covariance
trans.TransformTrack(state1,error_cov);

ASSERT_EQ(state1.g_.data_, state_trans.g_.data_);
ASSERT_EQ(state1.u_.data_, state_trans.u_.data_);
ASSERT_EQ(error_cov, error_cov_transformed);

// // Test



// Verify transformation data
Eigen::MatrixXd transform_data;
transform_data = Eigen::Matrix3d::Identity();

ASSERT_TRUE(trans.IsAcceptableTransformData(transform_data));

transform_data = Eigen::Matrix2d::Identity();

ASSERT_FALSE(trans.IsAcceptableTransformData(transform_data));

}

// Test the transformation for SE2_se2

TEST(TransformHomographyTest, SE2_se2_Transformation) {

///////////////////////////////////////
// Test the transform rotation function
///////////////////////////////////////

typedef Eigen::Matrix<double,5,5> MatCov;

TransformHomography<SE2_se2> trans;



// Construct a state and make sure the velocity isn't zero.
SE2_se2 state = SE2_se2::Random();
while (state.u_.data_.norm() < 1e-5) {
    state = SE2_se2::Random();
}

// Construct Homography.
SO2<double> R_psi = SO2<double>::Random();
Eigen::Matrix<double,3,1> t,n;
t.setRandom();
n << 0,0,1;

// std::cout << "n: " << std::endl << n << std::endl;
// std::cout << "t: " << std::endl << t << std::endl;
// std::cout << "tn: " << std::endl << t*n.transpose() << std::endl;

Eigen::Matrix3d H;
H.setZero();
H.block(0,0,2,2) = R_psi.data_;
H.block(0,2,3,1) = -t;
H(2,2) +=1;

Eigen::Matrix2d H1 = H.block(0,0,2,2);

// std::cout << "H: " << std::endl << H << std::endl << std::endl;
trans.Init();
trans.SetData(H);

// Transform the velocity and rotation and construct the proper rotation.
Eigen::Matrix<double,2,1> vel_transformed = trans.ConstructTranslationalVelTransform(state.g_.t_)*state.g_.R_*state.u_.p_; 
Eigen::Matrix2d R_transformed = trans.TransformRotation(vel_transformed);
Eigen::Matrix2d R_proper;
R_proper.block(0,0,2,1) = vel_transformed.normalized();
R_proper(0,1) = - R_proper(1,0);
R_proper(1,1) = R_proper(0,0);

ASSERT_EQ(R_proper, R_transformed);

// Test measurement transformation
Meas<double> m, m_transformed;
m.pose = Eigen::Matrix<double,2,1>::Random();
m.twist = Eigen::Matrix<double,2,1>::Random();
m.time_stamp = 0;
m.type = MeasurementTypes::SEN_POS_VEL;

m_transformed = m;
trans.TransformMeasurement(m_transformed);
ASSERT_NE(m.pose, m_transformed.pose);
ASSERT_NE(m.twist, m_transformed.twist);



/////////////////////////////////////////////////////
// Test the covariance portion of the transform track;
/////////////////////////////////////////////////////

double dt = 1e-8;
Eigen::Matrix<double,3,1> dt_t1;
Eigen::Matrix<double,3,1> dt_t2;
Eigen::Matrix<double,3,1> dt_r;
Eigen::Matrix<double,3,1> dt_p1;
Eigen::Matrix<double,3,1> dt_p2;
Eigen::Matrix<double,3,1> dt_w;

dt_t1 << dt,0,0;
dt_t2 << 0,dt,0;
dt_r  << 0,0,dt;
dt_p1 = dt_t1;
dt_p2 = dt_t2;
dt_w  = dt_r;


SE2_se2 state_base = SE2_se2::Random();                           // Get a random state and modify it to meet specifications.
Eigen::Matrix<double,2,1> v = state_base.u_.p_;
state_base.u_.p_ << v.norm(), 0;
v.normalize();
state_base.g_.R_ << v(0), -v(1), v(1), v(0);
MatCov cov;                                    // Construct a covariance matrix. Not used
MatCov cov_transformed_analytical, cov_original_extracted;
Eigen::Matrix<double,6,6> cov_original, cov_transformed_numerical;
cov_original.setZero();
cov_original_extracted.setRandom();
cov_original.block(0,0,4,4) = cov_original_extracted.block(0,0,4,4);
cov_original.block(5,0,1,4) = cov_original_extracted.block(4,0,1,4);
cov_original.block(0,5,4,1) = cov_original_extracted.block(0,4,4,1);
cov_original(5,5) = cov_original_extracted(4,4);

cov_transformed_analytical = cov_original_extracted;

cov_transformed_numerical = cov_original;
Eigen::Matrix<double,6,6> cov_transform_numerical;


// Construct perturbed states
SE2_se2 state_t1, state_t2, state_r, state_p1, state_p2, state_w;
state_t1 = state_base;
state_t2 = state_base;
state_r  = state_base;
state_p1 = state_base;
state_p2 = state_base;
state_w  = state_base;
state_t1.g_.OPlusEq(dt_t1);
state_t2.g_.OPlusEq(dt_t2);
state_r.g_.OPlusEq(dt_r);
state_p1.u_.data_+= dt_p1;
state_p2.u_.data_+= dt_p2;
state_w.u_.data_ += dt_w ;

trans.TransformTrack(state_t1, cov);
trans.TransformTrack(state_t2, cov);
trans.TransformTrack(state_r, cov);
trans.TransformTrack(state_p1, cov);
trans.TransformTrack(state_p2, cov);
trans.TransformTrack(state_w, cov);
trans.TransformTrack(state_base, cov_transformed_analytical);

cov_transform_numerical.block(0,0,6,1) = SE2_se2::OMinus(state_base,state_t1)/dt;
cov_transform_numerical.block(0,1,6,1) = SE2_se2::OMinus(state_base,state_t2)/dt;
cov_transform_numerical.block(0,2,6,1) = SE2_se2::OMinus(state_base,state_r)/dt;
cov_transform_numerical.block(0,3,6,1) = SE2_se2::OMinus(state_base,state_p1)/dt;
cov_transform_numerical.block(0,4,6,1) = SE2_se2::OMinus(state_base,state_p2)/dt;
cov_transform_numerical.block(0,5,6,1) = SE2_se2::OMinus(state_base,state_w)/dt;

MatCov cov_transform_numerical_extracted;
cov_transform_numerical_extracted.block(0,0,4,4) = cov_transform_numerical.block(0,0,4,4);
cov_transform_numerical_extracted.block(4,0,1,4) = cov_transform_numerical.block(5,0,1,4);
cov_transform_numerical_extracted.block(0,4,4,1) = cov_transform_numerical.block(0,5,4,1);
cov_transform_numerical_extracted(4,4) = cov_transform_numerical(5,5); 

cov_transformed_numerical = cov_transform_numerical*cov_original*cov_transform_numerical.transpose();

MatCov cov_extracted;
cov_extracted.block(0,0,4,4) = cov_transformed_numerical.block(0,0,4,4);
cov_extracted.block(4,0,1,4) = cov_transformed_numerical.block(5,0,1,4);
cov_extracted.block(0,4,4,1) = cov_transformed_numerical.block(0,5,4,1);
cov_extracted(4,4) = cov_transformed_numerical(5,5); 

// We undo the analytical transformation using the numerical transformations to test that the numerical is close to the analytical.
MatCov cov_undone = cov_transform_numerical_extracted.inverse()* cov_transformed_analytical * (cov_transform_numerical_extracted.transpose()).inverse();

ASSERT_LE( (cov_extracted - cov_transformed_analytical).norm(), 1e-7 );
ASSERT_LE( (cov_undone - cov_original_extracted).norm(), 1e-7 );

// Test the reprojection of the homography rotation 
// R_bad is a rotation that left SO2 and R good is the projection of R_bad back onto SO_2
// these matrices were generated from matlab
Eigen::Matrix2d R_bad, R_good, R_track;
R_bad << -0.112231140832595,	-0.990882177336245, 0.992832207585224,	-0.114921220083900;
R_good << -0.113765174035777,	-0.993507667396991, 0.993507667396991,	-0.113765174035778;
H.block(0,0,2,2) = R_bad;
trans.SetData(H);

R_track = R_good*state_t1.g_.R_;
trans.TransformTrack(state_t1, cov);

ASSERT_LT( (state_t1.g_.R_.transpose()*state_t1.g_.R_ - Eigen::Matrix2d::Identity() ).norm(),1e-12  );
ASSERT_LT( (R_track -state_t1.g_.R_ ).norm(),1e-12  );

}

} // namespace rransac
