#ifndef RRANSAC_COMMON_TRANSFORM_NULL_H_
#define RRANSAC_COMMON_TRANSFORM_NULL_H_

#include <Eigen/Core>
#include "common/measurement/measurement_base.h"
#include "common/transformations/transformation_base.h"

namespace rransac
{
/** \class TransformNULL
 * This transform class is used when the measurements and the track do not need to be transformed. 
*/

template<class tState>
class TransformNULL : public TransformBase<Eigen::Matrix<typename tState::DataType,Eigen::Dynamic, Eigen::Dynamic>, tState, Eigen::Matrix<typename tState::DataType,tState::dim_,tState::dim_>, TransformNULL<tState>> {

public:

typedef typename tState::DataType DataType;
typedef Eigen::Matrix<DataType,3,3> Mat3d;
typedef Mat3d MatData;


void DerivedInit() {
    this->transform_null_ = true;
}

/** 
 * The parent class sets the transformation data member variable. This derived class has the 
 * opportunity to perform other calculations using the data. 
 * @param data The data required to transform the measurements, states, and error covariance
 */ 
void DerivedSetData(const Mat3d data) {
    throw std::runtime_error("TransformNULL::SetData Not Implemented.");
}

/** 
 * Transforms the measurement using data_ from the previous surveillance frame to the current one.
 * @param meas The measurement to be transformed.
 */ 
void DerivedTransformMeasurement(Meas<double>& meas) const {
    throw std::runtime_error("TransformNULL::TransformMeasurement Not Implemented.");
}

/** 
 * Transforms the track using the transform data. i.e. transform the estimated 
 * state and error covariance.
 * @param state The track's state to be transformed.
 * @param cov   The track's error covariance to be transformed.
 */ 
void DerivedTransformTrack(tState& state, Eigen::Matrix<DataType,tState::dim_,tState::dim_>& cov) const {
    throw std::runtime_error("TransformNULL::TransformTrack Not Implemented.");
}



};





}// namespace rransac

#endif // RRANSAC_COMMON_TRANSFORM_NULL_H_
