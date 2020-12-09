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

template<class tState, class tMatCov = Eigen::Matrix<double,tState::dim_,tState::dim_> >
class TransformNULL : public TransformBase<Eigen::MatrixXd, tState, tMatCov, TransformNULL<tState>> {

public:


void DerivedInit() {
    this->transform_null_ = true;
}

/** 
 * The parent class sets the transformation data member variable. This derived class has the 
 * opportunity to perform other calculations using the data. 
 * @param data The data required to transform the measurements, states, and error covariance
 */ 
void DerivedSetData(Eigen::Matrix3d data) {
    throw std::runtime_error("TransformNULL::SetData Not Implemented.");
}

/** 
 * Transforms the measurement using data_ from the previous surveillance frame to the current one.
 * @param meas The measurement to be transformed.
 */ 
void DerivedTransformMeasurement(Meas& meas) {
    throw std::runtime_error("TransformNULL::TransformMeasurement Not Implemented.");
}

/** 
 * Transforms the track using the transform data. i.e. transform the estimated 
 * state and error covariance.
 * @param state The track's state to be transformed.
 * @param cov   The track's error covariance to be transformed.
 */ 
void DerivedTransformTrack(tState& state, Eigen::Matrix<double,tState::dim_,tState::dim_>& cov) {
    throw std::runtime_error("TransformNULL::TransformTrack Not Implemented.");
}



};





}// namespace rransac

#endif // RRANSAC_COMMON_TRANSFORM_NULL_H_
