#include <gtest/gtest.h>
#include <typeinfo>
#include <Eigen/Core>

#include "rransac/common/sources/source_base.h"
#include "rransac/common/sources/source_SE2_SE3_radar.h"
#include "rransac/common/sources/source_container.h"
#include "rransac/common/measurement/measurement_base.h"
#include "rransac/common/transformations/transformation_null.h"
#include "rransac/common/transformations/trans_radar_SE2_SE3_with_SE2_SE3.h"
#include "rransac/parameters.h"
#include "lie_groups/state.h"
#include "rransac/common/utilities.h"

using namespace rransac;
using namespace lie_groups;


typedef SourceRadarSE2SE3<SE2_se2,MeasurementTypes::SE2_SE3_RADAR,TransRadarSE2SE3WithSE2SE3> SourceSE2Radar;
typedef SourceRadarSE2SE3<SE3_se3,MeasurementTypes::SE2_SE3_RADAR,TransRadarSE2SE3WithSE2SE3> SourceSE3Radar;

using MyTypes = ::testing::Types<SourceSE2Radar,SourceSE3Radar>;


template<class _Source>
class SourceRadarSE2SE3Test : public ::testing::Test {

public:

typedef _Source Source;
typedef typename Source::Transformation Transformation;
typedef typename Source::State State;
typedef typename Source::DataType DataType;
typedef rransac::utilities::MatXT<DataType> MatXd;
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
std::vector<typename State::Vec_SC> e_;
JacobianHMat H_;
JacobianVMat V_;
State state_;
bool transform_state_ = false;
typename Source::TransformDataType EmptyMat_;
Source source_;


JacobianHMat ComputeNumericalH(const State& state) {

    JacobianHMat H;
    State state_perturbed = state;
    H.setZero();
    Measurement m = Source::GetEstMeas(state, transform_state_, EmptyMat_);
    Measurement m_perturbed;

    for (size_t ii = 0; ii < state_dim_; ++ii) {

        state_perturbed = state.OPlus(e_[ii]*dt_);
        m_perturbed = Source::GetEstMeas(state_perturbed, transform_state_, EmptyMat_);
        H.block(0,ii,total_meas_dim_,1) = Source::OMinus(m_perturbed,m)/dt_;

    }

    return H;

}


protected:

void SetUp() override {

SourceParameters source_params;
source_params.meas_cov_ = JacobianVMat::Identity();
source_params.spacial_density_of_false_meas_ = 0.1;
source_params.type_ = Source::measurement_type_;
source_params.gate_probability_ = 0.8;
source_params.probability_of_detection_ = 0.9;
source_params.source_index_ = 0;
this->source_.Init(source_params);

// construct the state deviations used to verify the Jacobians
for (int ii = 0; ii < state_dim_; ++ii) {
    typename State::Vec_SC v;
    int index = ii;
    v.setZero();

    if(Source::polar_coordinates_) {

        if(ii > 3) {
            index+=1;
        }

    } else {
        if (ii > 5) {
            index +=2;
        }
    }

    v(index,0) = 1;

    e_.push_back(v);;
}

// Construct the Jacobian V
V_.setIdentity();


}



};

//------------------------------------------------------------------------------------------------------------------------------------------------------------------

TYPED_TEST_SUITE(SourceRadarSE2SE3Test, MyTypes);

TYPED_TEST(SourceRadarSE2SE3Test, JacobianTestEstMeasTest){


typedef SourceRadarSE2SE3Test<TypeParam> TestType;

/////////////////////////
// Test the Get Est Meas;
/////////////////////////
size_t num_trials = 10;
double tolerance = 1e-6;
typename TestType::State state;
typename TestType::Measurement meas_test_case, meas;

for(size_t ii = 0; ii < num_trials; ++ii) {

    state = TestType::State::Random(10);

    meas = this->source_.GetEstMeas(state,this->transform_state_,this->EmptyMat_);

    ASSERT_DOUBLE_EQ(meas.pose(0,0),state.g_.t_.norm()); // Test the distance

    double r = meas.pose(0,0);
    double azimuth = meas.pose(1,0);

    if(TestType::polar_coordinates_){

        double x = r*cos(azimuth);
        double y = r*sin(azimuth);
        ASSERT_LT( fabs(x-state.g_.t_(0)), 1e-10);
        ASSERT_LT( fabs(y-state.g_.t_(1)), 1e-10);

   
    } else if(TestType::spherical_coordinates_) {

        double zenith = meas.pose(2,0);
        double x = r*cos(azimuth)*sin(zenith);
        double y = r*sin(azimuth)*sin(zenith);
        double z = r*cos(zenith);
        ASSERT_LT( fabs(x-state.g_.t_(0,0)),1e-10);
        ASSERT_LT( fabs(y-state.g_.t_(1,0)),1e-10);
        ASSERT_LT( fabs(z-state.g_.t_(2,0)),1e-10);

    }


}


////////////////////////////////////
// Test the Jacobian H and OMinus
////////////////////////////////////
typename TestType::JacobianHMat H_numerical, H_exact;
for (size_t ii = 0; ii < num_trials; ++ii) {

    state = TestType::State::Random(10);
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
// Test Generate Random Measurement
////////////////////////////////////

state = TestType::State::Random(10);
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
state = TestType::State::Random(10);
typename TestType::Measurement meas1_original, meas2_original;
typename TestType::Measurement meas1_transformed, meas2_transformed;
typename TestType::Transformation::TransformDataType transform_data1, transform_data2;
transform_data1 = TestType::Transformation::GetRandomTransform();
transform_data2 = TestType::Transformation::GetRandomTransform();
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