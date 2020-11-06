#ifndef RRANSAC_COMMON_TRANSFORMATION_H_
#define RRANSAC_COMMON_TRANSFORMATION_H_

#include <Eigen/Core>

namespace rransac
{
/** \class Transformation
 * When the global frame changes, the measurements and models need to be transformed into the current global frame.
 * This struct provides the necessary data in order to transform the measurements and models. It is provided by the
 * user and used by R-RANSAC.
 * 
 * The transformation being applied is dependent on the data type D, measurement type M, and state type S
*/

template <class D, class M, class S>
class TransformBase {

public:

/** 
 * Transforms the measurement using data_ from the previous surveillance frame to the current one.
 * @param meas The measurement to be transformed.
 */ 
virtual void TransformMeasurement(M& meas);

/** 
 * Transforms the state using data_ from the previous surveillance frame to the current one.
 * @param state The state to be transformed.
 */ 
virtual void TransformState(S& state);

/** 
 * Transforms the error covariance using data_ from the previous surveillance frame to the current one.
 * @param cov The error covariance to be transformed.
 */ 
virtual void TransformErrorCov(Eigen::MatrixXd& cov);
D data_;

};

#endif // RRANSAC_COMMON_TRANSFORMATION_H_
