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

template <class tData, class tState, class tMatCov, bool tNullTransform, class tDerived>
class TransformBase {

public:

typedef tData Data;                             /**< The object type of the data needed to transform the measurements and tracks. */
typedef tState State;                           /**< The state of the target. @see State. */
typedef tMatCov MatCov;                         /**< The object type of the tracks error covariance. */
typedef tDerived Derived;                       /**< The child class that implements the specific member functions. */
typedef typename tState::DataType DataType;     /**< The scalar object for the data. Ex. float, double, etc. */
typedef Eigen::Matrix<DataType, Eigen::Dynamic,Eigen::Dynamic> MatXd; /**< Dynamic Eigen Matrix */
static constexpr bool is_null_transform_ = tNullTransform; /**< Indicates if the derived transform call is the null transform */

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
tData GetData() const {
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

/** 
 * Transforms the state and error covariance using user provided transform data.
 * This function is useful when you need to transform a track to a measurement frame. 
 * @param[in] state The track's state to be transformed.
 * @param[in] cov   The track's error covariance to be transformed.
 * @param[in] transform_data The data used to transform the state and error covariance
 */ 
static void TransformTrack(State& state, MatCov& cov, const Data& transform_data) {
   tDerived::DerivedTransformTrack(state,cov,transform_data);
}

/** 
 * Transforms the state using user provided transform data.
 * This function is useful when you need to transform a track to a measurement frame. 
 * @param[in] state The track's state to be transformed.
 * @param[in] transform_data The data used to transform the state and error covariance
 */ 
static State TransformState(const State& state, const Data& transform_data) {
    return tDerived::DerivedTransformState(state,transform_data);
}

/** 
 * Returns the Jacobian of the transformation
 * @param[in] state The state of the target after it has been transformed using transform_data
 * @param[in] transform_data The data used in the transformation
 */ 
static MatXd GetTransformationJacobian(const State& state, const Data& transform_data) {
   return tDerived::DerivedGetTransformationJacobian(state,transform_data);
}

/**
 * Verifies that the transform data provided by the user is in the requested from.
 * @param transform_data The tranformation data to be tested. 
 */
static bool IsAcceptableTransformData(const MatXd& transform_data) {
    return tDerived::DerivedIsAcceptableTransformData(transform_data);
} 


private:
TransformBase()=default;
~TransformBase()=default;
friend Derived;
// private:

Data data_;

};


template <class tData, class tState, class tMatCov, class tDerived>
class TransformBase<tData,tState,tMatCov,true,tDerived> {

public:

typedef typename tState::DataType DataType;                           /**< The scalar object for the data. Ex. float, double, etc. */
typedef Eigen::Matrix<DataType, Eigen::Dynamic,Eigen::Dynamic> MatXd; /**< Dynamic Eigen Matrix */
typedef MatXd Data;                                                   /**< The object type of the data needed to transform the measurements and tracks. */
typedef tState State;                                                 /**< The state of the target. @see State. */
typedef MatXd MatCov;                                                 /**< The object type of the tracks error covariance. */
typedef tDerived Derived;                                             /**< The child class that implements the specific member functions. */


static constexpr bool is_null_transform_ = true; /**< Indicates if the derived transform call is the null transform */


/**
 * Used to initialize the transformation object.
 */ 
void Init() {}

/** 
 * Sets the transformation data member variable. This function will also call the derived classes set data in case
 * other stuff needs to be done. 
 * @param[in] data The data required to transform the measurements, states, and error covariance
 */ 
void SetData(const tData& data) {}


/** 
 * Returns the transformation data member variable.
 */ 
tData GetData() const {return data_;}
/** 
 * Transforms the measurement from the previous tracking frame to the current one.
 * @param[in] meas The measurement to be transformed.
 */ 
void TransformMeasurement(const Meas<typename tState::DataType>& meas) const {}

/** 
 * Transforms the track using the transform data. i.e. transform the estimated 
 * state and error covariance.
 * @param[in] state The track's state to be transformed.
 * @param[in] cov   The track's error covariance to be transformed.
 */ 
void TransformTrack(const tState& state, const MatXd& cov) const {}

/** 
 * Transforms the state and error covariance using user provided transform data.
 * This function is useful when you need to transform a track to a measurement frame. 
 * @param[in] state The track's state to be transformed.
 * @param[in] cov   The track's error covariance to be transformed.
 * @param[in] transform_data The data used to transform the state and error covariance
 */ 
static void TransformTrack(const tState& state, const MatXd& cov, const MatXd& transform_data) {}

/** 
 * Transforms the state using user provided transform data.
 * This function is useful when you need to transform a track to a measurement frame. 
 * @param[in] state The track's state to be transformed.
 * @param[in] transform_data The data used to transform the state and error covariance
 */ 
static tState TransformState(const tState& state, const MatXd& transform_data) {return state;}

/** 
 * Returns the Jacobian of the transformation
 * @param[in] state The state of the target after it has been transformed using transform_data
 * @param[in] transform_data The data used in the transformation
 */ 
static MatXd GetTransformationJacobian(const tState& state, const MatXd& transform_data) {return transform_data;}
/**
 * Verifies that the transform data provided by the user is in the requested from.
 * @param transform_data The tranformation data to be tested. 
 */
static bool IsAcceptableTransformData(const MatXd& transform_data) {return true;}

Data data_;


};



} // namespace rransac

#endif // RRANSAC_COMMON_TRANSFORMATION_TRANSFORMATION_BASE_H_
