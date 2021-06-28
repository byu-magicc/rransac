#ifndef RRANSAC_COMMON_SOURCES_SOURCE_NULL_H_
#define RRANSAC_COMMON_SOURCES_SOURCE_NULL_H_
#pragma once

#include <Eigen/Core>
#include "rransac/common/sources/source_base.h"
#include "lie_groups/state.h"
#include "rransac/common/transformations/transformation_null.h"
#include "rransac/common/utilities.h"

namespace rransac {


template<typename tState=lie_groups::R2_r2, MeasurementTypes tMeasurementType=MeasurementTypes::NUM_TYPES, typename tTransformation = TransformNULL<lie_groups::R2_r2>>
class SourceNull : public SourceBase<tState,tMeasurementType,tTransformation,SourceNull<tState,tMeasurementType,tTransformation>> {

public:

typedef Eigen::Matrix<typename tState::DataType,Eigen::Dynamic,Eigen::Dynamic> MatXd;
typedef double DataType;
typedef utilities::CompatibleWithModelNull ModelCompatibility;              /**<Indicates which model this source is compatible with. */
static constexpr MeasurementTypes measurement_type_ = tMeasurementType;     /**< The measurement type of the source. */


/** 
 * Initializes the measurement source by setting the parameters using SetParameters, calculating the non user specified parameters,
 * and initializing the Jacobians. By calling this constructor, it is assumed that a track is always inside the surveillance region.
 * This is a good assumption when all of the sources have the same surveillance region.
 * @param[in] params The source parameters.
 */ 
void Init(const SourceParameters& params) {
    throw std::runtime_error("SourceNull::Init Function Not Implemented, and shouldn't be called. ");    
}

/** 
 * Returns the jacobian of the observation function w.r.t. the states. 
 * @param[in] state A state of the target.
*/
static MatXd DerivedGetLinObsMatState(const tState& state) {
    throw std::runtime_error("SourceNull::DerivedGetLinObsMatState Function Not Implemented, and shouldn't be called. "); 
}   

/** 
 * Returns the jacobian of the observation function w.r.t. the sensor noise.
 * @param[in] state A state of the target.
 */
static MatXd DerivedGetLinObsMatSensorNoise(const tState& state)  {
    throw std::runtime_error("SourceNull::DerivedGetLinObsMatSensorNoise Function Not Implemented, and shouldn't be called. "); 

} 

/**
 *  Implements the observation function and returns an estimated measurement based on the state. 
 * @param[in] state A state of the target.
 */
static Meas<double> DerivedGetEstMeas(const tState& state)  {
    throw std::runtime_error("SourceNull::DerivedGetEstMeas Function Not Implemented, and shouldn't be called. "); 

} 

/**
 * Performs the OMinus operation between two measurement (m1 ominus m2) of the same type. In other words, this
 * method computes the geodesic distance between two measurements of the same type.
 * @param[in] m1 a measurement
 * @param[in] m2 a measurement
 */
static MatXd DerivedOMinus(const Meas<DataType>& m1, const Meas<DataType>& m2) {
    throw std::runtime_error("SourceNull::DerivedOMinus Function Not Implemented, and shouldn't be called. "); 
} 


/**
 * Generates a random measurement from a Gaussian distribution with mean defined by the state and standard deviation defined by meas_std. This
 * method is used primarily in simulations and tests.
 * @param[in] state    The state that serves as the mean of the Gaussian distribution.
 * @param[in] meas_std The measurement standard deviation.
 */ 
Meas<double> DerivedGenerateRandomMeasurement(const MatXd& meas_std, const tState& state) const {
    throw std::runtime_error("SourceNull::DerivedGenerateRandomMeasurement Function Not Implemented, and shouldn't be called. "); 

}

};



} // rransac

#endif // RRANSAC_COMMON_SOURCES_SOURCE_NULL_H_