#ifndef RRANSAC_COMMON_TRANSFORMATION_TRANSFORMATION_BASE_H_
#define RRANSAC_COMMON_TRANSFORMATION_TRANSFORMATION_BASE_H_

#include <Eigen/Core>
#include "common/measurement/measurement_base.h"
#include "state.h"


namespace rransac
{
/** \class Transformation
 * When the global frame changes, the measurements and models need to be transformed into the current global frame.
 * This struct provides the necessary data in order to transform the measurements and models. It is provided by the
 * user and used by R-RANSAC.
 * 
 * The transformation being applied is dependent on the data type, state type, and derived/child class
*/

template <class Data, class State, class Derived>
class TransformBase {

public:

/** 
 * Sets the transformation data member variable. This function will also call the derived classes set data in case
 * other stuff needs to be done. 
 * @param data The data required to transform the measurements, states, and error covariance
 */ 
void SetData(Data& data) {
    static_cast<Derived*>(this)->SetData(data);
}


/** 
 * Returns the transformation data member variable.
 */ 
Data GetData() {
    return data_;
}

/** 
 * Transforms the measurement using data_ from the previous surveillance frame to the current one.
 * @param meas The measurement to be transformed.
 */ 
void TransformMeasurement(Meas& meas) {
    static_cast<Derived*>(this)->TransformMeasurement(meas);
}

/** 
 * Transforms the track using the transform data. i.e. transform the estimated 
 * state and error covariance.
 * @param state The track's state to be transformed.
 * @param cov   The track's error covariance to be transformed.
 */ 
void TransformTrack(State& state, Eigen::Matrix<double,State::dim_,State::dim_>& cov) {
// void TransformTrack(State& state, Eigen::MatrixXd& cov) {
    static_cast<Derived*>(this)->TransformTrack(state,cov);
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
