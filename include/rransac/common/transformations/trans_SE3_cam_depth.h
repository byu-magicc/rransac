#ifndef RRANSAC_COMMON_TRANSFORATIONS_SE3_CAM_DEPTH
#define RRANSAC_COMMON_TRANSFORATIONS_SE3_CAM_DEPTH
#pragma once


#include <Eigen/Core>
#include "rransac/common/measurement/measurement_base.h"
#include "rransac/common/transformations/transformation_base.h"
#include "lie_groups/state.h"

namespace rransac
{


/** \class TransformSE3CamDepth
 * This transformation is used when target tracking is done on SE3 and the target is observed by a camera that returns the target's location
 * on the normalized image sphere and the relative depth of the target to the camera. 
 * The transformation data is an element of SE3 and transforms objects from the previous frame to the current frame. 
 * We assume that the velocity is expressed in the
 * body frame, so when the tracking frame moves, we only need to transform the pose of the target and not the velocity.
*/

template<class tState>
class TransformSE3CamDepth : public TransformBase<Eigen::Matrix<typename tState::DataType,4,4>, tState, Eigen::Matrix<typename tState::DataType, 10,10>, false, TransformSE3CamDepth<tState>> {

public:

typedef tState State;
typedef typename tState::DataType DataType; /**< The scalar object for the data. Ex. float, double, etc. */
typedef Eigen::Matrix<DataType,4,4> MatData;
typedef Eigen::Matrix<DataType,10,10> MatCov;
typedef Eigen::Matrix<DataType,Eigen::Dynamic, Eigen::Dynamic> MatXd;  


static_assert(std::is_same<tState,lie_groups::State< lie_groups::SE3,typename tState::DataType,tState::N>>::value, "TransformSE3CamDepth: The state is not supported");

/**
 * Used to initialize the object, but it doesn't need to initialize anyting.
 */ 
void DerivedInit() {};

/** 
 * Sets the data and the block components of the homography.
 * @param data The data required to transform the measurements, states, and error covariance. This should be a transformation on SE(3) from the previous tracking
 * frame to the current tracking frame.
 */ 
void DerivedSetData(const MatData data) {
    this->data_ = data;
}

/** 
 * Transforms the measurement from the previous tracking frame to the current one.
 * @param[in,out] meas The measurement to be transformed.
 */ 
void DerivedTransformMeasurement(Meas<DataType>& meas) const ;

/** 
 * Transforms the track provided that the state is SE3_se3.
 * Since the velocity is expressed in the body frame, it doesn't need to be transformed.
 * Since the error covariance is on the error state, it doesn't need to be transformed. 
 * @param[in,out] state The track's state to be transformed.
 * @param[in,out] cov   The track's error covariance.
 */ 
void DerivedTransformTrack(tState& state, MatCov& cov) const {
    state.g_.data_ = this->data_*state.g_.data_;
}


/** 
 * Transforms the state and error covariance using user provided transform data.
 * The transform data is an element of SE(3) and transforms the state from the current
 * tracking frame to a measurement frame.
 * @param[in] state The track's state to be transformed.
 * @param[in] cov   The track's error covariance to be transformed.
 * @param[in] transform_data The data used to transform the state and error covariance
 */ 
static void DerivedTransformTrack(State& state, MatCov& cov, const MatData& transform_data) {
   state.g_.data_ = transform_data*state.g_.data_;
}

/** 
 * Transforms the state using user provided transform data.
 * The transform data is an element of SE(3) and transforms the state from the current
 * tracking frame to a measurement frame.
 * @param[in] state The track's state to be transformed.
 * @param[in] transform_data The data used to transform the state and error covariance
 */ 
static State DerivedTransformState(const State& state, const MatData& transform_data) {
    State transformed_state;
    transformed_state.u_.data_= state.u_.data_;
    transformed_state.g_.data_ = transform_data*state.g_.data_;
    return transformed_state;
}

/** 
 * Returns the Jacobian of the transformation
 * The transform data is an element of SE(3) and transforms the state from the current
 * tracking frame to a measurement frame.
 * @param[in] state The state of the target after it has been transformed using transform_data
 * @param[in] transform_data The data used in the transformation
 */ 
static MatXd DerivedGetTransformationJacobian(const State& state, const MatData& transform_data) {
   return Eigen::Matrix<DataType,4,4>::Identity();
}


/**
 * Verifies that the transform data provided by the user is in the requested from.
 * @param transform_data The tranformation data to be tested. 
 */
static bool DerivedIsAcceptableTransformData(const Eigen::MatrixXd& transform_data) {

    bool correct;

    // transform data should be a 3x3 matrix
    if(transform_data.cols() != 4 || transform_data.rows() != 4) {
        correct = false;
    } else if (transform_data.determinant() !=1) { // The determinant should be 1
        correct = false;
    } else if(transform_data.block(0,0,3,3)*transform_data.block(0,0,3,3).transpose() != Eigen::Matrix3d::Identity()) { // Make sure it is a rotation matrix
        correct = false;
    } else if (transform_data.block(3,0,1,2).norm() != 0) { // These cells should be zero
        correct = false;
    } else {
        correct = true;
    }



    return correct;
} 


};

template<class tState>
void TransformSE3CamDepth<tState>::DerivedTransformMeasurement(Meas<DataType>& meas) const {

    if(!meas.transform_state) {
        meas.transform_state = true;
        meas.transform_data_m_t = this->data_;
        meas.transform_data_t_m = this->data_.inverse();
    } else {
        meas.transform_data_m_t = this->data_ * meas.transform_data_m_t;
        meas.transform_data_t_m = meas.transform_data_t_m * this->data_.inverse();
    }

    // DataType d_oa_a = meas.pose(0,0);
    // Eigen::Matrix<DataType,3,1> s_oa_a = meas.pose.block(1,0,3,1);
    // Eigen::Matrix<DataType,3,1> sd_oa_a = meas.twist;
    // Eigen::Matrix<DataType,3,3> I = Eigen::Matrix<DataType,3,3>::Identity();

    // Eigen::Matrix<DataType,3,3> R_ab = this->data_.block(0,0,3,3);
    // Eigen::Matrix<DataType,3,1> t_ab_b = this->data_.block(0,3,3,1);

    // Eigen::Matrix<DataType,3,1> t_oa_a = d_oa_a*s_oa_a;
    // Eigen::Matrix<DataType,3,1> t_ob_b = R_ab*t_oa_a + t_ab_b;
    // Eigen::Matrix<DataType,3,1> v_oa_a = (I*pow(d_oa_a,2) - (t_oa_a*t_oa_a.transpose())).inverse()*pow(d_oa_a,3)*sd_oa_a;

    // std::cout << "I: " << std::endl << I << std::endl;
    // std::cout << "d_oa_a: " << std::endl << d_oa_a << std::endl;
    // std::cout << "t_oa_a: " << std::endl << t_oa_a << std::endl;
    // std::cout << "sd_oa_a:  " << std::endl << sd_oa_a << std::endl;


    // DataType d_ob_b = sqrt(  (t_ob_b.transpose()*t_ob_b)(0,0)   );
    // Eigen::Matrix<DataType,3,1> s_ob_b = t_ob_b/d_ob_b;
    // Eigen::Matrix<DataType,3,1> sd_ob_b = R_ab*v_oa_a/d_ob_b - (t_ob_b*t_ob_b.transpose()*R_ab*v_oa_a)/pow(d_ob_b,3);

    // meas.pose(0,0) = d_ob_b;
    // meas.pose.block(1,0,3,1) = s_ob_b;
    // meas.twist = sd_ob_b;

    // std::cout << "v: " << std::endl << v_oa_a << std::endl;


}





} // namespace rransac
#endif // RRANSAC_COMMON_TRANSFORATIONS_SE3_CAM_DEPTH