#ifndef RRANSAC_COMMON_TRANSFORM_NULL_H_
#define RRANSAC_COMMON_TRANSFORM_NULL_H_
#pragma once


#include <Eigen/Core>

#include "common/measurement/measurement_base.h"
#include "common/transformations/transformation_base.h"

namespace rransac
{
/** \class TransformNULL
 * This transform class is used when the measurements and the track do not need to be transformed. 
 * None of it's member functions does anything so it should be optimized out. @see TransformBase
*/

template<class tState>
class TransformNULL : public TransformBase<Eigen::Matrix<typename tState::DataType,Eigen::Dynamic, Eigen::Dynamic>, tState, Eigen::Matrix<typename tState::DataType,tState::dim_,tState::dim_>, TransformNULL<tState>> {

public:

typedef typename tState::DataType DataType;  /**< The scalar object for the data. Ex. float, double, etc. */
typedef Eigen::Matrix<DataType,3,3> Mat3d;
typedef Mat3d MatData;


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
void DerivedTransformMeasurement(Meas<double>& meas) const {}

/** 
 * Doesn't transform the track.
 * @param[in] state The track's state to be transformed.
 * @param[in] cov   The track's error covariance to be transformed.
 */ 
void DerivedTransformTrack(tState& state, Eigen::Matrix<DataType,tState::dim_,tState::dim_>& cov) const {}



};





}// namespace rransac

#endif // RRANSAC_COMMON_TRANSFORM_NULL_H_
