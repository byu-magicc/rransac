#ifndef RRANSAC_COMMON_TRANSFORM_NULL_H_
#define RRANSAC_COMMON_TRANSFORM_NULL_H_
#pragma once


#include <Eigen/Core>

#include "rransac/common/measurement/measurement_base.h"
#include "rransac/common/transformations/transformation_base.h"

namespace rransac
{
/** \class TransformNULL
 * This transform class is used when the measurements and the track do not need to be transformed. 
 * None of it's member functions does anything so it should be optimized out. @see TransformBase
*/

template<class tState>
class TransformNULL : public TransformBase<Eigen::Matrix<typename tState::DataType,Eigen::Dynamic, Eigen::Dynamic>, tState, Eigen::Matrix<typename tState::DataType,tState::dim_,tState::dim_>, true, TransformNULL<tState>> {

public:

typedef typename tState::DataType DataType;  /**< The scalar object for the data. Ex. float, double, etc. */
typedef Eigen::Matrix<DataType,3,3> Mat3d;
typedef Mat3d MatData;
typedef Eigen::Matrix<DataType,Eigen::Dynamic, Eigen::Dynamic> MatXd;


void DerivedInit() {}

/** 
 * Doesn't set the data. 
 * @param data The data required to transform the measurements, states, and error covariance
 */ 
void DerivedSetData(const Mat3d data) {}

/** 
 * Doesn't transform the measurements.
 * @param meas The measurement to be transformed.
 */ 
void DerivedTransformMeasurement(Meas<DataType>& meas) const {}

/** 
 * Doesn't transform the track.
 * @param[in] state The track's state to be transformed.
 * @param[in] cov   The track's error covariance to be transformed.
 */ 
void DerivedTransformTrack(tState& state, Eigen::Matrix<DataType,tState::dim_,tState::dim_>& cov) const {}

/** 
 * Transforms the state and error covariance using user provided transform data.
 * This function is useful when you need to transform a track to a measurement frame. 
 * @param[in] state The track's state to be transformed.
 * @param[in] cov   The track's error covariance to be transformed.
 * @param[in] transform_data The data used to transform the state and error covariance
 */ 
static void DerivedTransformTrack(tState& state, MatData& cov, const Mat3d& transform_data) {}

/** 
 * Transforms the state using user provided transform data.
 * This function is useful when you need to transform a track to a measurement frame. 
 * @param[in] state The track's state to be transformed.
 * @param[in] transform_data The data used to transform the state and error covariance
 */ 
static tState DerivedTransformState(const tState& state, const Mat3d& transform_data) { return state; }

/** 
 * Returns the Jacobian of the transformation
 * @param[in] state The state of the target after it has been transformed using transform_data
 * @param[in] transform_data The data used in the transformation
 */ 
static MatXd DerivedGetTransformationJacobian(const tState& state, const Mat3d& transform_data) {
   return Eigen::Matrix<DataType,1,1>::Zero();
}


/**
 * Verifies that the transform data provided by the user is in the requested from.
 * @param transform_data The tranformation data to be tested. 
 */
static bool DerivedIsAcceptableTransformData(const Eigen::MatrixXd& transform_data) {
    return true;
} 


};





}// namespace rransac

#endif // RRANSAC_COMMON_TRANSFORM_NULL_H_
