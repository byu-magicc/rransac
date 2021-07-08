#ifndef RRANSAC_COMMON_TRANSFORMATION_TRANSFORMATION_BASE_H_
#define RRANSAC_COMMON_TRANSFORMATION_TRANSFORMATION_BASE_H_
#pragma once


#include <Eigen/Core>

#include "lie_groups/state.h"
#include "rransac/common/measurement/measurement_base.h"



namespace rransac
{

/**
 * \class TransformDerivedTraits
 * The traits of the derived class necessary for the base transform class.
 */ 
template<typename _State, typename _TransformDataType, typename _MatCov, bool _IsNullTransform>
struct TransformDerivedTraits {
    typedef _State State;                              /**< The state of the target. @see State.*/
    typedef _TransformDataType TransformDataType;      /**< The object type of the data needed to transform the measurements and tracks. */
    typedef _MatCov MatCov;                            /**< The object type of the tracks error covariance. */
    static constexpr bool is_null_transform_ = _IsNullTransform;
};

/** \class TransformBase
 * When the target tracking frame changes, the measurements and tracks need to be transformed from the previous tracking frame to the current
 * tracking frame. This class uses the curiously recurring template pattern design methodology to specify the API for the derived classes. It
 * it responsible for setting the transformation data and transforming measurements and tracks. 
 * 
*/

template <typename _TransformDerivedTraits, template <typename > typename _Derived>
class TransformBase {

public:

typedef typename _TransformDerivedTraits::TransformDataType TransformDataType;  /**< The object type of the data needed to transform the measurements and tracks. */
typedef typename _TransformDerivedTraits::State State;      /**< The state of the target. @see State. */
typedef typename _TransformDerivedTraits::MatCov MatCov;    /**< The object type of the tracks error covariance. */
typedef _TransformDerivedTraits DerivedTraits;              /**< The traits from the derived class. */
typedef _Derived<State> Derived;                            /**< The child class that implements the specific member functions. */
typedef typename State::DataType DataType;                  /**< The scalar object for the data. Ex. float, double, etc. */
typedef Meas<DataType,TransformDataType> Measurement;       /**< The measurement data type. */
static constexpr bool is_null_transform_ = _TransformDerivedTraits::is_null_transform_; /**< Indicates if the derived transform call is the null transform */

template<typename _State>
using TransformTemplate = _Derived<_State>;

/**
 * Used to initialize the transformation object.
 */ 
void Init() {
    static_cast<Derived*>(this)->DerivedInit();
}

/** 
 * Sets the transformation data member variable. This function will also call the derived classes set data in case
 * other stuff needs to be done. 
 * @param[in] data The data required to transform the measurements, states, and error covariance
 */ 
void SetData(const TransformDataType& data) {
    static_cast<Derived*>(this)->DerivedSetData(data);
}


/** 
 * Returns the transformation data member variable.
 */ 
TransformDataType GetData() const {
    return data_;
}

/** 
 * Transforms the measurement from the previous tracking frame to the current one.
 * @param[in] meas The measurement to be transformed.
 */ 
void TransformMeasurement(Measurement& meas) const {
    static_cast<const Derived*>(this)->DerivedTransformMeasurement(meas);
}

/** 
 * Transforms the measurement from according to the transformation provided.
 * @param[in] meas The measurement to be transformed.
 */ 
static Measurement TransformMeasurement(const Measurement& meas, const TransformDataType& data) {
    return Derived::DerivedTransformMeasurement(meas,data);
}

/** 
 * Transforms the track using the transform data. i.e. transform the estimated 
 * state and error covariance.
 * @param[in] state The track's state to be transformed.
 * @param[in] cov   The track's error covariance to be transformed.
 */ 
void TransformTrack(State& state, MatCov& cov) const {
    static_cast<const Derived*>(this)->DerivedTransformTrack(state,cov);
}

/** 
 * Transforms the state and error covariance using user provided transform data.
 * This function is useful when you need to transform a track to a measurement frame. 
 * @param[in] state The track's state to be transformed.
 * @param[in] cov   The track's error covariance to be transformed.
 * @param[in] transform_data The data used to transform the state and error covariance
 */ 
static void TransformTrack(State& state, MatCov& cov, const TransformDataType& transform_data) {
   Derived::DerivedTransformTrack(state,cov,transform_data);
}

/** 
 * Transforms the state using user provided transform data.
 * This function is useful when you need to transform a track to a measurement frame. 
 * @param[in] state The track's state to be transformed.
 * @param[in] transform_data The data used to transform the state and error covariance
 */ 
static State TransformState(const State& state, const TransformDataType& transform_data) {
    return Derived::DerivedTransformState(state,transform_data);
}

/** 
 * Returns the Jacobian of the transformation
 * @param[in] state The state of the target after it has been transformed using transform_data
 * @param[in] transform_data The data used in the transformation
 */ 
static MatCov GetTransformationJacobian(const State& state, const TransformDataType& transform_data) {
   return Derived::DerivedGetTransformationJacobian(state,transform_data);
}

/**
 * Verifies that the transform data provided by the user is in the requested from.
 * @param transform_data The tranformation data to be tested. 
 */
static bool IsAcceptableTransformData(const TransformDataType& transform_data) {
    return Derived::DerivedIsAcceptableTransformData(transform_data);
} 

/**
 * Generates random transform data. The function can use the parameter scalar in order to 
 * generate a larger distribution of random transformations.
 * @param scalar A scalar used to generate a larger distribution. 
 */ 
static TransformDataType GetRandomTransform(const DataType scalar = static_cast<DataType>(1.0)){
    return Derived::DerivedGetRandomTransform(scalar);
}


private:
TransformBase()=default;
~TransformBase()=default;
friend Derived;
// private:

TransformDataType data_;

};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Template specialization for the null transform.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename _State, typename _TransformDataType, typename _MatCov, template<typename > typename _Derived>
class TransformBase<TransformDerivedTraits<_State,_TransformDataType,_MatCov,true>,_Derived> {

public:
typedef TransformDerivedTraits<_State,_TransformDataType,_MatCov,true> DerivedTraits;
typedef typename DerivedTraits::TransformDataType TransformDataType;           /**< The object type of the data needed to transform the measurements and tracks. */
typedef typename DerivedTraits::State State;               /**< The state of the target. @see State. */
typedef typename DerivedTraits::MatCov MatCov;             /**< The object type of the tracks error covariance. */
typedef _Derived<State> Derived;                           /**< The child class that implements the specific member functions. */
typedef typename State::DataType DataType;                 /**< The scalar object for the data. Ex. float, double, etc. */
typedef Meas<DataType,TransformDataType> Measurement;      /**< The measurement data type. */
static constexpr bool is_null_transform_ = true;           /**< Indicates if the derived transform call is the null transform */



/**
 * Used to initialize the transformation object.
 */ 
void Init() {}

/** 
 * Sets the transformation data member variable. This function will also call the derived classes set data in case
 * other stuff needs to be done. 
 * @param[in] data The data required to transform the measurements, states, and error covariance
 */ 
void SetData(const TransformDataType& data) {}


/** 
 * Returns the transformation data member variable.
 */ 
TransformDataType GetData() const {return data_;}
/** 
 * Transforms the measurement from the previous tracking frame to the current one.
 * @param[in] meas The measurement to be transformed.
 */ 
void TransformMeasurement(const Measurement& meas) const {}

/** 
 * Transforms the measurement from according to the transformation provided.
 * @param[in] meas The measurement to be transformed.
 */ 
static Measurement TransformMeasurement(const Measurement& meas, const TransformDataType& data) {return meas;};


/** 
 * Transforms the track using the transform data. i.e. transform the estimated 
 * state and error covariance.
 * @param[in] state The track's state to be transformed.
 * @param[in] cov   The track's error covariance to be transformed.
 */ 
void TransformTrack(const State& state, const MatCov& cov) const {}

/** 
 * Transforms the state and error covariance using user provided transform data.
 * This function is useful when you need to transform a track to a measurement frame. 
 * @param[in] state The track's state to be transformed.
 * @param[in] cov   The track's error covariance to be transformed.
 * @param[in] transform_data The data used to transform the state and error covariance
 */ 
static void TransformTrack(const State& state, const MatCov& cov, const TransformDataType& transform_data) {}

/** 
 * Transforms the state using user provided transform data.
 * This function is useful when you need to transform a track to a measurement frame. 
 * @param[in] state The track's state to be transformed.
 * @param[in] transform_data The data used to transform the state and error covariance
 */ 
static State TransformState(const State& state, const TransformDataType& transform_data) {return state;}

/** 
 * Returns the Jacobian of the transformation of the track
 * @param[in] state The state of the target after it has been transformed using transform_data
 * @param[in] transform_data The data used in the transformation
 */ 
static MatCov GetTransformationJacobian(const State& state, const TransformDataType& transform_data) {return transform_data;}
/**
 * Verifies that the transform data provided by the user is in the requested from.
 * @param transform_data The tranformation data to be tested. 
 */
static bool IsAcceptableTransformData(const TransformDataType& transform_data) {return true;}

/**
 * Generates random transform data. The function can use the parameter scalar in order to 
 * generate a larger distribution of random transformations.
 * @param scalar A scalar used to generate a larger distribution. 
 */ 
static TransformDataType GetRandomTransform(const DataType scalar = static_cast<DataType>(1.0)){
    throw std::runtime_error("TransformBase::GetRandomTransform: This function shouldn't be called. ");
    // return TransformDataType::Zero();
}


TransformDataType data_;


};



} // namespace rransac

#endif // RRANSAC_COMMON_TRANSFORMATION_TRANSFORMATION_BASE_H_
