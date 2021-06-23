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
class TransformSE3CamDepth : public TransformBase<Eigen::Matrix<typename tState::DataType,4,4>, tState, Eigen::Matrix<typename tState::DataType, 10,10>, TransformSE3CamDepth<tState>> {

public:

typedef typename tState::DataType DataType; /**< The scalar object for the data. Ex. float, double, etc. */
typedef Eigen::Matrix<double,4,4> MatData;
typedef Eigen::Matrix<double,10,10> MatCov;

static_assert(std::is_same<tState,lie_groups::State< lie_groups::SE3,typename tState::DataType,tState::N>>::value, "TransformSE3CamDepth: The state is not supported");

/**
 * Used to initialize the object, but it doesn't need to initialize anyting.
 */ 
void DerivedInit() {};

/** 
 * Sets the data and the block components of the homography.
 * @param data The data required to transform the measurements, states, and error covariance
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




};

template<class tState>
void TransformSE3CamDepth<tState>::DerivedTransformMeasurement(Meas<DataType>& meas) const {

    if(!meas.state_transform_data) {
        meas.state_transform_data = true;
        meas.trans_data = this->data_;
    } else {
        meas.trans_data = this->data_ * meas.trans_data;
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