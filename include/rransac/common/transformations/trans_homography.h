#ifndef RRANSAC_COMMON_TRANS_HOMOGRAPHY_R2_H_
#define RRANSAC_COMMON_TRANS_HOMOGRAPHY_R2_H_

#include <Eigen/Core>
#include "common/measurement/measurement_base.h"
#include "lie_groups/liegroups/state.h"
#include "common/transformations/transformation_base.h"

namespace rransac
{
/** \class TransHomography
 * Transforms the measurements, states and error covariance using the homography. The homography is Eigen::Matrix3d
 * the state is R2 or SE2 and the measurement is of type _______
*/

template<class S>
class TransformHomography : public TransformBase<Eigen::Matrix3d, Meas, S> {

public:

/** 
 * Transforms the measurement using data_ from the previous surveillance frame to the current one.
 * @param meas The measurement to be transformed.
 */ 
void TransformMeasurement(Meas& meas) override;

/** 
 * Transforms the state using data_ from the previous surveillance frame to the current one.
 * @param state The state to be transformed.
 */ 
template<class S>
void TransformState(S& state) override {
    throw std::runtime_error("TransformHomography: The state is not supported.");
}

/** 
 * Transforms the error covariance using data_ from the previous surveillance frame to the current one.
 * @param cov The error covariance to be transformed.
 */ 
template<class S>
void TransformErrorCov(Eigen::MatrixXd& cov) override {
 throw std::runtime_error("TransformHomography: The state is not supported.");
}


};

//-------------------------------------------------------------------------------------------

template<>
void TransformHomography<lie_gorups::R2_r2>::TransformState(lie_gorups::R2_r2& state){

}

template<>
void TransformHomography<lie_gorups::SE2_se2>::TransformState(lie_gorups::SE2_se2& state){
    
}

}

#endif // RRANSAC_COMMON_TRANS_HOMOGRAPHY_R2_H_
