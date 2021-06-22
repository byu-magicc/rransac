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
class TransformSE3CamDepth : public TransformBase<Eigen::Matrix<typename tState::DataType,3,3>, tState, Eigen::Matrix<typename tState::DataType,TransformHomographyCovDim(tState::g_type_::dim_),TransformHomographyCovDim(tState::g_type_::dim_)>, TransformHomography<tState>> {

public:

typedef typename tState::DataType DataType; /**< The scalar object for the data. Ex. float, double, etc. */
typedef Eigen::Matrix<double,4,4> MatData;

std::static_assert(std::is_same<tState,lie_groups::State<tState::G,tState::DataType,tState::N>::value, "TransformSE3CamDepth: The state is not supported");

/**
 * Used to initialize the object, but it doesn't need to initialize anyting.
 */ 
void DerivedInit() {};

/** 
 * Sets the data and the block components of the homography.
 * @param data The data required to transform the measurements, states, and error covariance
 */ 
void DerivedSetData(const Mat3d data) {
    this->data_ = data;
}

/** 
 * Transforms the measurement from the previous tracking frame to the current one.
 * @param[in] meas The measurement to be transformed.
 */ 
void DerivedTransformMeasurement(Meas<DataType>& meas) const {

    if (meas.twist.rows() != 0) {
        Mat2d&& tmp = ConstructTranslationalVelTransform(meas.pose);
        meas.twist = tmp*meas.twist;
    }

    meas.pose = TransformPosition(meas.pose);
}

/** 
 * Transforms the track provided that the state is SE3_se3.
 * Since the velocity is expressed in the body frame, it doesn't need to be transformed.
 * Since the error covariance is on the error state, it doesn't need to be transformed. 
 * @param[in] state The track's state to be transformed.
 * @param[in] cov   The track's error covariance.
 */ 
void DerivedTransformTrack(tState& state, MatCov& cov) const {
    state.g_.data_ = this->data_*state_.g_.data_;
}




};





} // namespace rransac
#endif // RRANSAC_COMMON_TRANSFORATIONS_SE3_CAM_DEPTH