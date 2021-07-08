#include <gtest/gtest.h>
#include <chrono> 
#include <iostream>


#include "rransac/common/transformations/trans_radar_R2_R3_with_SE2_SE3.h"
#include "lie_groups/state.h"


namespace rransac {

using namespace lie_groups;



// Create state types
typedef State<Rn,double,2,1> StateR2_1;
typedef State<Rn,double,2,2> StateR2_2;
typedef State<Rn,double,2,3> StateR2_3;

typedef State<Rn,double,3,1> StateR3_1;
typedef State<Rn,double,3,2> StateR3_2;
typedef State<Rn,double,3,3> StateR3_3;

using MyTypes = ::testing::Types<StateR2_1,StateR2_2,StateR2_3,StateR3_1,StateR3_2,StateR3_3>;
// using MyTypes = ::testing::Types<StateR2_1,StateR2_2,StateR2_3>;

template<class _State>
class TransRadarR2R3WithSE2SE3Test : public ::testing::Test {

public:

typedef TransRadarR2R3WithSE2SE3<_State> Transformation;
typedef typename Transformation::State State;                              /**< The State type being used. */
typedef typename Transformation::DataType DataType;                        /**< The scalar data type. */
typedef typename Transformation::TransformDataType TransformDataType;      /**< The transform data type being used. It is either an element of SE2 for R2 or SE3 for R3. */
typedef typename Transformation::MatCov MatCov;                            /**< The covariance type of the track, and the transform jacobian type. */
typedef typename Transformation::Measurement Measurement;                  /**< The measurement type. */
static constexpr unsigned int meas_dim_ = State::Group::dim_;
typedef Eigen::Matrix<DataType,meas_dim_,1> VecMeas;


Transformation transformation_;

//------------------------------------------

TransformDataType GetRandomTransformation(){
    Eigen::MatrixXd random_transformation;
    if (Transformation::polar_coordinates_ ) {
        random_transformation = lie_groups::SE2<DataType>::Random();
    } else {
        random_transformation = lie_groups::SE3<DataType>::Random(); 
    } 
    return random_transformation;
}


protected:

void SetUp() override {

}


};


TYPED_TEST_SUITE(TransRadarR2R3WithSE2SE3Test, MyTypes);

TYPED_TEST(TransRadarR2R3WithSE2SE3Test, NonImplementedFunctions){

typedef TransRadarR2R3WithSE2SE3Test<TypeParam> TestType;
typename TestType::State state;
typename TestType::MatCov cov;
typename TestType::Measurement meas;
typename TestType::TransformDataType transform_data;

ASSERT_ANY_THROW(this->transformation_.SetData(transform_data));
ASSERT_ANY_THROW(this->transformation_.TransformMeasurement(meas));
ASSERT_ANY_THROW(this->transformation_.TransformTrack(state,cov));


}

//------------------------------------------------------------------

TYPED_TEST(TransRadarR2R3WithSE2SE3Test, TransformMeasurement){

typedef TransRadarR2R3WithSE2SE3Test<TypeParam> TestType;
typename TestType::State state;
typename TestType::MatCov cov;
typename TestType::Measurement meas;
typename TestType::TransformDataType transform_data;
typename TestType::VecMeas meas_data, meas_data_transformed, pos;
meas_data.setRandom();
meas_data= meas_data.cwiseAbs()*M_PI;
meas.pose = meas_data;
transform_data = this->GetRandomTransformation();
ASSERT_TRUE(this->transformation_.IsAcceptableTransformData(transform_data));

if(TestType::meas_dim_ == 2) {

    double r = meas_data(0);
    double azimuth = meas_data(1);
    pos << r*cos(azimuth), r*sin(azimuth);
    pos = transform_data.block(0,0,2,2)*pos + transform_data.block(0,2,2,1);
    meas_data_transformed(0) = pos.norm();
    meas_data_transformed(1) = atan2(pos(1),pos(0));


} else {
    double r = meas_data(0);
    double azimuth = meas_data(1);
    double zenith = meas_data(2);
    pos << r*cos(azimuth)*sin(zenith), r*sin(azimuth)*sin(zenith), r*cos(zenith);
    pos = transform_data.block(0,0,3,3)*pos+transform_data.block(0,3,3,1);
    meas_data_transformed(0) = pos.norm();
    meas_data_transformed(1) = atan2(pos(1),pos(0));
    meas_data_transformed(2) = atan2( pos.block(0,0,2,1).norm(), pos(2));
}

// std::cout << "meas_data " << std::endl << meas_data << std::endl;
// std::cout << "transform data: " << std::endl << transform_data << std::endl;
// std::cout << "pos: " << std::endl << pos << std::endl;
// std::cout << "meas transformed: " << std::endl << meas_data_transformed << std::endl;

ASSERT_LT ( (meas_data_transformed - this->transformation_.TransformMeasurement(meas,transform_data).pose).norm(),1e-10 );

// Ensure that they are the same when using the identity transform
transform_data.setIdentity();
ASSERT_LT ( (meas.pose - this->transformation_.TransformMeasurement(meas,transform_data).pose).norm(),1e-10 );


}

//------------------------------------------------------------------

TYPED_TEST(TransRadarR2R3WithSE2SE3Test, TransformTrack){

typedef TransRadarR2R3WithSE2SE3Test<TypeParam> TestType;
typename TestType::State state, state_original;
typename TestType::MatCov cov, cov_original;
typename TestType::Measurement meas;
typename TestType::TransformDataType transform_data;
typename TestType::VecMeas meas_data, meas_data_transformed, pos;

state_original.g_.data_.setRandom();
state_original.u_.data_.setRandom();
cov_original.setRandom();

transform_data = this->transformation_.GetRandomTransform(10);
ASSERT_TRUE(this->transformation_.IsAcceptableTransformData(transform_data));

transform_data.setIdentity();
ASSERT_TRUE(this->transformation_.IsAcceptableTransformData(transform_data));


state = state_original;
cov = cov_original;

// Identity transformation
this->transformation_.TransformTrack(state,cov,transform_data);
ASSERT_EQ(state.g_.data_, state_original.g_.data_);
ASSERT_EQ(state.u_.data_, state_original.u_.data_);
ASSERT_EQ(cov,cov_original);

// translating the data wont affect the velocity or the covariance.
transform_data.block(0,TestType::State::Group::dim_,TestType::State::Group::dim_,1).setRandom();
this->transformation_.TransformTrack(state,cov,transform_data);
ASSERT_EQ(state.g_.data_, state_original.g_.data_+transform_data.block(0,TestType::State::Group::dim_,TestType::State::Group::dim_,1));
ASSERT_EQ(state.u_.data_, state_original.u_.data_);
ASSERT_EQ(cov,cov_original);


transform_data = this->GetRandomTransformation();
const Eigen::Matrix<typename TestType::DataType, TestType::State::Group::dim_,TestType::State::Group::dim_>& rotation = transform_data.block(0,0,TestType::State::Group::dim_,TestType::State::Group::dim_);
const Eigen::Matrix<typename TestType::DataType, TestType::State::Group::dim_,1> position_offset = transform_data.block(0,TestType::State::Group::dim_,TestType::State::Group::dim_,1); 

Eigen::Matrix<typename TestType::DataType, TestType::State::Group::dim_,TestType::State::Group::dim_> identity;
identity.setIdentity();
ASSERT_LT (  (rotation*rotation.transpose() - identity).norm(),1e-10);

// construct the covariance transformation
typename TestType::MatCov TJ;
TJ.setZero();
for (int ii = 0; ii < TestType::Transformation::num_tangent_spaces_+1; ++ii) {
    TJ.block(TestType::Transformation::group_dim_*ii,TestType::Transformation::group_dim_*ii,TestType::Transformation::group_dim_,TestType::Transformation::group_dim_) = rotation;
}

state = state_original;
cov = cov_original;
this->transformation_.TransformTrack(state,cov,transform_data);
ASSERT_LT(  (TJ*TJ.transpose() - TestType::Transformation::MatCov::Identity()).norm(),1e-10); // Ensure that it is unitary
ASSERT_LT( (cov_original - TJ.transpose()*cov*TJ).norm(),1e-10);


ASSERT_LT( (TJ.block(TestType::Transformation::group_dim_,TestType::Transformation::group_dim_,TestType::Transformation::group_dim_*TestType::Transformation::num_tangent_spaces_, TestType::Transformation::group_dim_*TestType::Transformation::num_tangent_spaces_)   *state_original.u_.data_ - state.u_.data_).norm(),1e-10);


}




} // rransac