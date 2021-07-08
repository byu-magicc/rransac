#include <gtest/gtest.h>
#include <typeinfo>
#include <Eigen/Core>


#include "rransac/common/sources/source_base.h"
#include "rransac/common/sources/source_R2_R3_radar.h"
#include "rransac/parameters.h"
#include "lie_groups/state.h"
#include "rransac/common/transformations/transformation_null.h"
#include "rransac/common/transformations/trans_radar_R2_R3_with_SE2_SE3.h"
#include "rransac/parameters.h"
#include "lie_groups/lie_groups/SE3.h"
#include "lie_groups/lie_groups/SE2.h"

namespace rransac {

using namespace lie_groups;



// Create state types
typedef State<Rn,double,2,1> StateR2_1;
typedef State<Rn,double,2,2> StateR2_2;
typedef State<Rn,double,2,3> StateR2_3;

typedef State<Rn,double,3,1> StateR3_1;
typedef State<Rn,double,3,2> StateR3_2;
typedef State<Rn,double,3,3> StateR3_3;

// Create Source Types
typedef SourceRadarR2_R3<StateR2_1,MeasurementTypes::R2_R3_RADAR,TransRadarR2R3WithSE2SE3> SourceR2_1;
typedef SourceRadarR2_R3<StateR2_2,MeasurementTypes::R2_R3_RADAR,TransRadarR2R3WithSE2SE3> SourceR2_2;
typedef SourceRadarR2_R3<StateR2_3,MeasurementTypes::R2_R3_RADAR,TransRadarR2R3WithSE2SE3> SourceR2_3;
typedef SourceRadarR2_R3<StateR2_1,MeasurementTypes::R2_R3_RADAR_DEPTH_DERIV,TransRadarR2R3WithSE2SE3> SourceR2_1_Deriv;
typedef SourceRadarR2_R3<StateR2_2,MeasurementTypes::R2_R3_RADAR_DEPTH_DERIV,TransRadarR2R3WithSE2SE3> SourceR2_2_Deriv;
typedef SourceRadarR2_R3<StateR2_3,MeasurementTypes::R2_R3_RADAR_DEPTH_DERIV,TransRadarR2R3WithSE2SE3> SourceR2_3_Deriv;
typedef SourceRadarR2_R3<StateR3_1,MeasurementTypes::R2_R3_RADAR,TransRadarR2R3WithSE2SE3> SourceR3_1;
typedef SourceRadarR2_R3<StateR3_2,MeasurementTypes::R2_R3_RADAR,TransRadarR2R3WithSE2SE3> SourceR3_2;
typedef SourceRadarR2_R3<StateR3_3,MeasurementTypes::R2_R3_RADAR,TransRadarR2R3WithSE2SE3> SourceR3_3;
typedef SourceRadarR2_R3<StateR3_1,MeasurementTypes::R2_R3_RADAR_DEPTH_DERIV,TransRadarR2R3WithSE2SE3> SourceR3_1_Deriv;
typedef SourceRadarR2_R3<StateR3_2,MeasurementTypes::R2_R3_RADAR_DEPTH_DERIV,TransRadarR2R3WithSE2SE3> SourceR3_2_Deriv;
typedef SourceRadarR2_R3<StateR3_3,MeasurementTypes::R2_R3_RADAR_DEPTH_DERIV,TransRadarR2R3WithSE2SE3> SourceR3_3_Deriv;


using MyTypes = ::testing::Types<SourceR2_1,SourceR2_2,SourceR2_3,SourceR2_1_Deriv,SourceR2_2_Deriv,SourceR2_3_Deriv,SourceR3_1,SourceR3_2,SourceR3_3,SourceR3_1_Deriv,SourceR3_2_Deriv,SourceR3_3_Deriv>;
// using MyTypes = ::testing::Types<SourceR3_2_Deriv>;

template<class _Source>
class SourceR2R3RadarTest : public ::testing::Test {


public:

typedef _Source Source;
typedef typename Source::State State;
typedef typename Source::DataType DataType;
typedef utilities::MatXT<DataType> MatXd;
static constexpr unsigned int state_dim_      = Source::state_dim_;                        /**< The dimension of the state. */
static constexpr unsigned int meas_pose_rows_ = Source::meas_pose_rows_;                   /**< The number of rows in the pose measurement. */
static constexpr unsigned int meas_pose_cols_ = Source::meas_pose_cols_;                   /**< The number of columns in the pose measurement. */
static constexpr unsigned int meas_twist_rows_= Source::meas_twist_rows_;                  /**< The number of rows in the twist measurement. */
static constexpr unsigned int meas_twist_cols_= Source::meas_twist_cols_;                  /**< The number of columns in the twist measurement. */
static constexpr unsigned int total_meas_dim_ = meas_pose_rows_+meas_twist_rows_;          /**< The total measurement dimension. */
static constexpr unsigned int group_dim_      = Source::State::Group::dim_;                /**< The dimensions of the group part of the state. */
static constexpr unsigned int algebra_dim_    = Source::State::Algebra::total_num_dim_;    /**< The dimensions of the algebra part of the state. */
static constexpr bool polar_coordinates_      = Source::polar_coordinates_;                /**< This variable is true if the measurement is in polar coordinates. */
static constexpr bool spherical_coordinates_  = Source::spherical_coordinates_;            /**< This value is true if the measurement is in spherical coordinates. */
static constexpr bool has_vel_                = Source::Base::has_vel_;                    /**< This value is true if the measurement contains velocity component. */
typedef Eigen::Matrix<DataType,state_dim_,1> StateVec;
typedef Eigen::Matrix<DataType,total_meas_dim_,1> MeasVec;
typedef Eigen::Matrix<DataType,total_meas_dim_,state_dim_> JacobianHMat; /**< The data type for the jacobians.*/
typedef Eigen::Matrix<DataType,total_meas_dim_,total_meas_dim_> JacobianVMat; /**< The data type for the jacobians.*/
typedef typename Source::Measurement Measurement;

double dt_ = 1e-8;
std::vector<StateVec> e_;
JacobianHMat H_;
JacobianVMat V_;
State state_;
bool transform_state_ = false;
typename Source::TransformDataType EmptyMat_;
Source source_;


//-----------------------------------------------------------

State GetRandomState() {
    State state = State::Random();
    state.g_.data_*=10;
    state.u_.data_*=10;
    return state;
}

//-------------------------------
StateVec ConvertStateToVec(const State& state) {
    StateVec v;
    v.block(0,0,group_dim_,1) = state.g_.data_;
    v.block(group_dim_,0,algebra_dim_,1) = state.u_.data_;
    return v;
}

//-------------------------------
MeasVec ConvertMeasToVec(const Measurement& meas) {
    MeasVec v;
    v.block(0,0,meas_pose_rows_,1) = meas.pose;
    if (Source::Base::has_vel_){
        v.block(meas_pose_rows_,0,1,1) = meas.twist;
    }
    return v;
}

//-----------------------------------------------------------

Eigen::Matrix<DataType,Eigen::Dynamic,Eigen::Dynamic> GetRandomTransform() {

    if (polar_coordinates_) {
        return SE2<DataType>::Random();
    } else {
        return SE3<DataType>::Random();
    }


}

//-----------------------------------------------------------

JacobianHMat ComputeNumericalH(const State& state) {
    JacobianHMat H;
    MeasVec v;
    for (int ii = 0; ii < state_dim_; ++ii){
        v = (ConvertMeasToVec(Source::GetEstMeas(state.OPlus(e_[ii]*dt_),transform_state_,EmptyMat_)) - ConvertMeasToVec(Source::GetEstMeas(state,transform_state_,EmptyMat_)))/dt_;
        H.block(0,ii,total_meas_dim_,1) = v;
    }

    return H;
}





protected:

void SetUp() override {

// Setup the source
SourceParameters params;
params.meas_cov_ = JacobianVMat::Identity();
params.spacial_density_of_false_meas_ = 0.1;
params.type_ = Source::measurement_type_;
params.gate_probability_ = 0.8;
params.probability_of_detection_ = 0.9;
params.source_index_ = 0;
this->source_.Init(params);



// construct the state deviations used to verify the Jacobians
for (int ii = 0; ii < state_dim_; ++ii) {
    StateVec v;
    v.setZero();
    v(ii,0) = 1;

    e_.push_back(v);;
}

// Construct the Jacobian V
V_.setIdentity();


}



};



TYPED_TEST_SUITE(SourceR2R3RadarTest, MyTypes);

TYPED_TEST(SourceR2R3RadarTest, JacobianTestEstMeasTest){


typedef SourceR2R3RadarTest<TypeParam> TestType;

/////////////////////////
// Test the Get Est Meas;
/////////////////////////
int num_trials = 1;
double tolerance = 1e-6;
typename TestType::Measurement meas_test_case, meas;
typename TestType::State state;
for (int ii = 0; ii < num_trials; ++ii) {

    state=this->GetRandomState();

    meas = this->source_.GetEstMeas(state,this->transform_state_,this->EmptyMat_);

    ASSERT_DOUBLE_EQ(meas.pose(0,0),state.g_.data_.norm()); // Test the range

    double r = meas.pose(0,0);
    double azimuth = meas.pose(1,0);

    if(TestType::polar_coordinates_){

        double x = r*cos(azimuth);
        double y = r*sin(azimuth);
        ASSERT_LT( fabs(x-state.g_.data_(0,0)), 1e-10);
        ASSERT_LT( fabs(y-state.g_.data_(1,0)), 1e-10);

   
    } else if(TestType::spherical_coordinates_) {

        double zenith = meas.pose(2,0);
        double x = r*cos(azimuth)*sin(zenith);
        double y = r*sin(azimuth)*sin(zenith);
        double z = r*cos(zenith);
        ASSERT_LT( fabs(x-state.g_.data_(0,0)),1e-10);
        ASSERT_LT( fabs(y-state.g_.data_(1,0)),1e-10);
        ASSERT_LT( fabs(z-state.g_.data_(2,0)),1e-10);

    }

    // Test the change in depth
    if(TestType::has_vel_) {
        double r0 = r;
        double r1 = (state.g_.data_ + state.u_.data_.block(0,0,TestType::group_dim_,1)*this->dt_).norm();
        double rd_numerical = (r1-r0)/this->dt_;
        ASSERT_LT( fabs(rd_numerical - meas.twist(0,0)), tolerance);
    }
    


}




///////////////////////
// Test the Jacobian H
///////////////////////
typename TestType::JacobianHMat H_numerical, H_exact;
for (int ii = 0; ii < num_trials; ++ii) {

    state = this->GetRandomState();
    H_numerical = this->ComputeNumericalH(state);
    H_exact = this->source_.GetLinObsMatState(state,this->transform_state_,this->EmptyMat_);

    // std::cout << "H_numerical: " << std::endl << H_numerical << std::endl;
    // std::cout << "H_exact: " << std::endl << H_exact << std::endl;

    ASSERT_LT( (H_numerical - H_exact).norm(),tolerance);

}

///////////////////////
// Test the Jacobian V
///////////////////////

ASSERT_EQ(this->V_, this->source_.GetLinObsMatSensorNoise(state,this->transform_state_,this->EmptyMat_));


////////////////////////////////////
// Test OMinus
////////////////////////////////////
typename TestType::Measurement meas1, meas2;
state = this->GetRandomState();
typename TestType::StateVec state_vec = this->ConvertStateToVec(state);
typename TestType::StateVec random_state_vec1, random_state_vec2;
typename TestType::MeasVec meas_diff_truth, meas_diff_source;
random_state_vec1.setRandom();
random_state_vec2.setRandom();
meas1 = this->source_.GetEstMeas(state.OPlus(random_state_vec1),this->transform_state_,this->EmptyMat_);
meas2 = this->source_.GetEstMeas(state.OPlus(random_state_vec2),this->transform_state_,this->EmptyMat_);

meas_diff_truth.block(0,0,TestType::meas_pose_rows_,1) = meas1.pose - meas2.pose; 
if(TestType::has_vel_){
meas_diff_truth.block(TestType::meas_pose_rows_,0,1,1) = meas1.twist - meas2.twist; 
    
}

meas_diff_source = this->source_.OMinus(meas1,meas2);
ASSERT_LT( (meas_diff_truth - meas_diff_source).norm(), 1e-10);





////////////////////////////////////
// Test Generate Random Measurement
////////////////////////////////////

state = this->GetRandomState();
const int num_rand = 10000;
double std_scalar = 0.1;
typename TestType::JacobianVMat meas_std = TestType::JacobianVMat::Identity()*std_scalar;

std::vector<typename TestType::MeasVec> error(num_rand);
typename TestType::MeasVec error_mean;
error_mean.setZero();

// std::cout << "std1: " << std::endl << std1 << std::endl;

// Generate the random measurement, calculate the error between the state and measurement, and get the mean of the error
for (unsigned long int ii = 0; ii < num_rand; ++ii) {

    meas = this->source_.GenerateRandomMeasurement(meas_std,state,this->transform_state_,this->EmptyMat_);
    error[ii] = this->source_.OMinus(meas, this->source_.GetEstMeas(state, this->transform_state_, this->EmptyMat_));
    error_mean+=error[ii];

}
error_mean /=num_rand;

// Calculate the covariance
typename TestType::JacobianVMat cov;
cov.setZero();
for (unsigned long int ii = 0; ii < num_rand; ++ii){

    cov += (error_mean - error[ii])*(error_mean - error[ii]).transpose();

}
cov /= num_rand;





ASSERT_LE( (error_mean - TestType::MeasVec::Zero()).norm(), 0.1);
ASSERT_LE( (cov - meas_std*meas_std).norm(), 0.1);


////////////////////////////////////
// Test Distance
////////////////////////////////////
Parameters sys_params;
state = this->GetRandomState();
typename TestType::Measurement meas1_original, meas2_original;
typename TestType::Measurement meas1_transformed, meas2_transformed;
typename TestType::Source::Transformation::TransformDataType transform_data1, transform_data2;
transform_data1 = this->GetRandomTransform();
transform_data2 = this->GetRandomTransform();
meas1_original = this->source_.GetEstMeas(state, false, transform_data1);
meas1_original.transform_meas = false;
meas1_transformed = this->source_.GetEstMeas(state, true, transform_data1);
meas1_transformed.transform_meas = true;
meas1_transformed.transform_data_m_t = transform_data1.inverse();
meas2_original = this->source_.GetEstMeas(state, false, transform_data2);
meas2_transformed = this->source_.GetEstMeas(state, true, transform_data2);
meas2_transformed.transform_meas = true;
meas2_transformed.transform_data_m_t = transform_data2.inverse();

ASSERT_LT( fabs(this->source_.GetSpatialDistance(meas1_transformed,meas2_transformed,sys_params)),1e-6);
ASSERT_LT( fabs(this->source_.GetSpatialDistance(meas2_original,meas2_transformed,sys_params)   ),1e-6);
ASSERT_LT( fabs(this->source_.GetSpatialDistance(meas1_transformed,meas1_original,sys_params)   ),1e-6);
ASSERT_LT( fabs(this->source_.GetSpatialDistance(meas1_original,meas2_original,sys_params)      ),1e-6);


}




} //namespace