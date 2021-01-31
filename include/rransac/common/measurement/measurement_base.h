#ifndef RRANSAC_COMMON_MEASUREMENT_MEASUREMENT_BASE_H_
#define RRANSAC_COMMON_MEASUREMENT_MEASUREMENT_BASE_H_
#pragma once



#include <Eigen/Core>
#include <string>

#include "rransac/parameters.h"

namespace rransac
{


/** \class MeasurementTypes
 * Lists the different types of measurements available. These types are used for indexing so do not change the order of
 * any existing ones except that NUM_TYPES must be the last one. Currently RRANSAC doesn't support
 * a mix of measurement types with different target spaces. For example, you cannot have 
 * measurements that are both SEN_POSE and RN_POS; however, you can have measurements that are
 * RN_POS and RN_POS_VEL. 
 */ 
enum MeasurementTypes {
    RN_POS,                   // The target space is RN and the position is measured
    RN_POS_VEL,               // The target space is RN and the position and velocity is measured
    SEN_POSE,                 // The target space is SEN and the pose is measured and is an element of SEN. The pose is in the surveillance region frame
    SEN_POSE_TWIST,           // The target space is SEN and the pose and twist is measured and are elements of SEN and sen. The pose is in the surveillance region frame and the twist is in the body frame 
    SEN_POS,                  // The target space is SEN and the position is measured. The position is in the surveillance region frame
    SEN_POS_VEL,              // The target space is SEN and the position and velocity is measured. The position and velocity are measured in the surveillance region frame
    NUM_TYPES
    // SO3_ATT
};


/** \struct Meas
 * This struct contains information regarding a measurement. The user uses this object to supply R-RANSAC with measurements.
 * The user is responsible to provide Meas::data, Meas::time_stamp and Meas::source_id. If the measurement covariance is not provided
 * to the Source object, the user must also provide Meas::meas_cov. The other member variables are set by R-RANSAC.
*/
template<typename tDataType=double>
struct Meas
{
    typedef Eigen::Matrix<tDataType,Eigen::Dynamic, Eigen::Dynamic> MatX;
    
    double time_stamp;          /**< The time the measurement was taken. */
    unsigned int source_index;  /**< When a new source is added, it is added to the vector System::sources_. The source indexes the vector to grab the corresponding source. So make sure it's the right index. */
    MeasurementTypes type;      /** < The measurement type @see MeasurementTypes */
    MatX pose;       /**< The part of the measurement corresponding to the pose of the target. (position, attitude, or both). */
    MatX twist;      /**< The part of the measurement corresponding to the derivative of the pose. (velocity, angular rates, or both). */
  
    // These member variables are reserved
    double likelihood;          /**< The likelihood that the measurement came from the phenomenon it was associated with. This value is set during the data
                                      association process.*/
    double weight;              /**< The weight of the measurement when updating the model is was associated with. This value is set during the data association
                                     process. */
    double vol;
};


} // namespace rransac

#endif // RRANSAC_COMMON_MEASUREMENT_MEASUREMENT_BASE_H_
