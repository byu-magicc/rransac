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

template<class _State>
class TransformNULL : public TransformBase< TransformDerivedTraits<_State,Eigen::Matrix<typename _State::DataType,Eigen::Dynamic, Eigen::Dynamic>,Eigen::Matrix<typename _State::DataType,Eigen::Dynamic, Eigen::Dynamic>,true>, TransformNULL> {

public:

typedef TransformBase< TransformDerivedTraits<_State,Eigen::Matrix<typename _State::DataType,Eigen::Dynamic, Eigen::Dynamic>,Eigen::Matrix<typename _State::DataType,Eigen::Dynamic, Eigen::Dynamic>,true>, TransformNULL> Base;
typedef typename Base::State State;                                      /**< The State type being used. */
typedef typename Base::DataType DataType;                                /**< The scalar data type. */
typedef typename Base::TransformDataType TransformDataType;              /**< The transform data type being used. It is either an element of SE2 for R2 or SE3 for R3. */
typedef typename Base::MatCov MatCov;                                    /**< The covariance type of the track, and the transform jacobian type. */
typedef typename Base::Measurement Measurement;                          /**< The measurement type. */

void DerivedInit() {}

/** 
 * Doesn't set the data. 
 * @param data The data required to transform the measurements, states, and error covariance
 */ 
void DerivedSetData(const TransformDataType data) {}

/** 
 * Doesn't transform the measurements.
 * @param meas The measurement to be transformed.
 */ 
void DerivedTransformMeasurement(const Measurement& meas) const {}

/** 
 * Doesn't transform the measurements.
 * @param meas The measurement to be transformed.
 */ 
static void DerivedTransformMeasurement(const Measurement& meas, const TransformDataType& transform_data) {}

/** 
 * Doesn't transform the track.
 * @param[in] state The track's state to be transformed.
 * @param[in] cov   The track's error covariance to be transformed.
 */ 
void DerivedTransformTrack(const State& state, const MatCov& cov) const {}

/** 
 * Transforms the state and error covariance using user provided transform data.
 * This function is useful when you need to transform a track to a measurement frame. 
 * @param[in] state The track's state to be transformed.
 * @param[in] cov   The track's error covariance to be transformed.
 * @param[in] transform_data The data used to transform the state and error covariance
 */ 
static void DerivedTransformTrack(const State& state, const MatCov& cov, const TransformDataType& transform_data) {}

/** 
 * Transforms the state using user provided transform data.
 * This function is useful when you need to transform a track to a measurement frame. 
 * @param[in] state The track's state to be transformed.
 * @param[in] transform_data The data used to transform the state and error covariance
 */ 
static State DerivedTransformState(const State& state, const TransformDataType& transform_data) { return state; }

/** 
 * Returns the Jacobian of the transformation
 * @param[in] state The state of the target after it has been transformed using transform_data
 * @param[in] transform_data The data used in the transformation
 */ 
static MatCov DerivedGetTransformationJacobian(const State& state, const TransformDataType& transform_data) {
   
     return MatCov::Identity();
}


/**
 * Verifies that the transform data provided by the user is in the requested from.
 * @param transform_data The tranformation data to be tested. 
 */
static bool DerivedIsAcceptableTransformData(const TransformDataType& transform_data) {
    return true;
} 

/**
 * Generates random transform data. The function can use the parameter scalar in order to 
 * generate a larger distribution of random transformations.
 * @param scalar A scalar used to generate a larger distribution. 
 */ 
static TransformDataType DerivedGetRandomTransform(const DataType scalar = static_cast<DataType>(1.0)){
    return TransformDataType::Zero();
}


};





}// namespace rransac

#endif // RRANSAC_COMMON_TRANSFORM_NULL_H_
