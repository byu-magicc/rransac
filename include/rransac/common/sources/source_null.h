#ifndef RRANSAC_COMMON_SOURCES_SOURCE_NULL_H_
#define RRANSAC_COMMON_SOURCES_SOURCE_NULL_H_
#pragma once

#include <Eigen/Core>
#include "rransac/common/sources/source_base.h"
#include "lie_groups/state.h"
#include "rransac/common/transformations/transformation_null.h"
#include "rransac/common/utilities.h"

namespace rransac {

using namespace utilities;

template<typename _State=lie_groups::R2_r2, MeasurementTypes _MeasurementType=MeasurementTypes::NUM_TYPES, template <typename > typename _Transformation = TransformNULL>
class SourceNull : public SourceBase<SourceDerivedTraits<_State,MatXT<typename _State::DataType>,MatXT<typename _State::DataType>,_Transformation,1,1,0,0,1,0,false,_MeasurementType,utilities::CompatibleWithModelNull>,SourceNull> {

public:

    typedef SourceBase<SourceDerivedTraits<_State,MatXT<typename _State::DataType>,MatXT<typename _State::DataType>,_Transformation,1,1,0,0,1,0,false,_MeasurementType,utilities::CompatibleWithModelNull>,SourceNull> Base;
    typedef typename Base::State State;                                            /**< The state of the target. @see State. */
    typedef typename Base::DataType DataType;                                      /**< The scalar object for the data. Ex. float, double, etc. */
    typedef typename Base::MatH MatH;                                              /**< The object type of the Jacobian H. */
    typedef typename Base::MatV MatV;                                              /**< The object type of the Jacobians V. */
    typedef typename Base::Transformation Transformation;                          /**< The transformation used to transform the measurements and tracks. */
    typedef typename Base::MatMeasCov MatMeasCov;                                  /**< The data type of the measurement covariance. */
    typedef typename Base::VecMeas VecMeas;                                        /**< The data type of the measurement covariance. */
    static constexpr unsigned int meas_pose_rows_  = Base::meas_pose_rows_;        /**< The number of rows in the pose measurement. */
    static constexpr unsigned int meas_pose_cols_  = Base::meas_pose_cols_;        /**< The number of columns in the pose measurement. */
    static constexpr unsigned int meas_twist_rows_ = Base::meas_twist_rows_;       /**< The number of rows in the twist measurement. */
    static constexpr unsigned int meas_twist_cols_ = Base::meas_twist_cols_;       /**< The number of columns in the twist measurement. */
    static constexpr unsigned int meas_pose_dim_   = Base::meas_pose_dim_;         /**< The measurement pose dimension. */
    static constexpr unsigned int meas_twist_dim_  = Base::meas_twist_dim_;        /**< The measurement twist dimension. */
    static constexpr unsigned int total_meas_dim_  = Base::total_meas_dim_;        /**< The total measurement dimension. */
    static constexpr bool has_vel_ = Base::has_vel_;                               /**< Indicates if the measurement contains velocity.  */
    static constexpr MeasurementTypes measurement_type_ = Base::measurement_type_; /**< The measurement type of the source. */
    typedef typename Base::TransformDataType TransformDataType;                    /**< The error type of the difference between two measurements. */
    typedef typename Base::Measurement Measurement;                                /**< The measurement data type. */

/** 
 * Initializes the measurement source by setting the parameters using SetParameters, calculating the non user specified parameters,
 * and initializing the Jacobians. By calling this constructor, it is assumed that a track is always inside the surveillance region.
 * This is a good assumption when all of the sources have the same surveillance region.
 * @param[in] params The source parameters.
 */ 
void DerivedInit(const SourceParameters& params) {
    throw std::runtime_error("SourceNull::Init Function Not Implemented, and shouldn't be called. ");    
}


/** 
 * Returns the jacobian of the observation function w.r.t. the states. 
 * @param[in] state A state of the target.
*/
static MatH DerivedGetLinObsMatState(const State& state) {
    throw std::runtime_error("SourceNull::DerivedGetLinObsMatState Function Not Implemented, and shouldn't be called. "); 
}   

/** 
 * Returns the jacobian of the observation function w.r.t. the sensor noise.
 * @param[in] state A state of the target.
 */
static MatV DerivedGetLinObsMatSensorNoise(const State& state)  {
    throw std::runtime_error("SourceNull::DerivedGetLinObsMatSensorNoise Function Not Implemented, and shouldn't be called. "); 

} 

/**
 *  Implements the observation function and returns an estimated measurement based on the state. 
 * @param[in] state A state of the target.
 */
static Measurement DerivedGetEstMeas(const State& state)  {
    throw std::runtime_error("SourceNull::DerivedGetEstMeas Function Not Implemented, and shouldn't be called. "); 

} 

/**
 * Performs the OMinus operation between two measurement (m1 ominus m2) of the same type. In other words, this
 * method computes the geodesic distance between two measurements of the same type.
 * @param[in] m1 a measurement
 * @param[in] m2 a measurement
 */
static VecMeas DerivedOMinus(const Measurement& m1, const Measurement& m2) {
    throw std::runtime_error("SourceNull::DerivedOMinus Function Not Implemented, and shouldn't be called. "); 
} 


/**
 * Generates a random measurement from a Gaussian distribution with mean defined by the state and standard deviation defined by meas_std. This
 * method is used primarily in simulations and tests.
 * @param[in] state    The state that serves as the mean of the Gaussian distribution.
 * @param[in] meas_std The measurement standard deviation.
 */ 
Measurement DerivedGenerateRandomMeasurement(const MatMeasCov& meas_std, const State& state) const {
    throw std::runtime_error("SourceNull::DerivedGenerateRandomMeasurement Function Not Implemented, and shouldn't be called. "); 

}

};


//-----------------------------------------------------------------------------------------------------------------------------


template<typename _Source>
struct IsSourceNull {
    static constexpr bool value = false;
};

template<typename _State, MeasurementTypes _MeasurementType, template<typename > typename _Transformation>
struct IsSourceNull<SourceNull<_State,_MeasurementType,_Transformation>>{
      static constexpr bool value = true; 
};




} // rransac

#endif // RRANSAC_COMMON_SOURCES_SOURCE_NULL_H_