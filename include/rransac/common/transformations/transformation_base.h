#ifndef RRANSAC_COMMON_TRANSFORMATION_TRANSFORMATION_BASE_H_
#define RRANSAC_COMMON_TRANSFORMATION_TRANSFORMATION_BASE_H_
#pragma once


#include <Eigen/Core>

#include "lie_groups/state.h"
#include "rransac/common/measurement/measurement_base.h"



namespace rransac
{


/** \class TransformBase
 * When the target tracking frame changes, the measurements and tracks need to be transformed from the previous tracking frame to the current
 * tracking frame. This class uses the curiously recurring template pattern design methodology to specify the API for the derived classes. It
 * it responsible for setting the transformation data and transforming measurements and tracks. 
 * 
*/

template <class tData, class tState, class tMatCov, class tDerived>
class TransformBase {

public:

typedef tData Data;                             /**< The object type of the data needed to transform the measurements and tracks. */
typedef tState State;                           /**< The state of the target. @see State. */
typedef tMatCov MatCov;                         /**< The object type of the tracks error covariance. */
typedef tDerived Derived;                       /**< The child class that implements the specific member functions. */
typedef typename tState::DataType DataType;     /**< The scalar object for the data. Ex. float, double, etc. */


/**
 * Used to initialize the transformation object.
 */ 
void Init() {
    static_cast<tDerived*>(this)->DerivedInit();
}

/** 
 * Sets the transformation data member variable. This function will also call the derived classes set data in case
 * other stuff needs to be done. 
 * @param[in] data The data required to transform the measurements, states, and error covariance
 */ 
void SetData(const Data& data) {
    static_cast<tDerived*>(this)->DerivedSetData(data);
}


/** 
 * Returns the transformation data member variable.
 */ 
tData GetData() {
    return data_;
}

/** 
 * Transforms the measurement from the previous tracking frame to the current one.
 * @param[in] meas The measurement to be transformed.
 */ 
void TransformMeasurement(Meas<DataType>& meas) const {
    static_cast<const tDerived*>(this)->DerivedTransformMeasurement(meas);
}

/** 
 * Transforms the track using the transform data. i.e. transform the estimated 
 * state and error covariance.
 * @param[in] state The track's state to be transformed.
 * @param[in] cov   The track's error covariance to be transformed.
 */ 
void TransformTrack(State& state, MatCov& cov) const {
    static_cast<const tDerived*>(this)->DerivedTransformTrack(state,cov);
}


private:
TransformBase()=default;
~TransformBase()=default;
friend Derived;
// private:

Data data_;

};



} // namespace rransac

#endif // RRANSAC_COMMON_TRANSFORMATION_TRANSFORMATION_BASE_H_
