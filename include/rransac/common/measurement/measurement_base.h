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
    SE3_CAM_DEPTH,            // The target space is SE3 and the targets location on the normalized image sphere and its depth to the camera are observed.
    NUM_TYPES
};


/** \struct Meas
 * This struct contains information regarding a measurement. The user uses this object to supply R-RANSAC with measurements.
 * The user is responsible to provide Meas::pose, Meas::twist (iff applicable), Meas::time_stamp, Meas::source_index, and Meas::type. 
*/
template<typename tDataType=double>
struct Meas
{
    typedef tDataType DataType;                                            /**< The scalar object for the data. Ex. float, double, etc. */
    typedef Eigen::Matrix<tDataType,Eigen::Dynamic, Eigen::Dynamic> MatX;  /**< The object type for the pose and twist of the measurement. */
    
    double time_stamp;          /**< The time the measurement was taken. */
    unsigned int source_index;  /**< When a new source is added, it is added to the vector System::sources_. The source indexes the vector to grab the corresponding source. So make sure it's the right index. */
    MeasurementTypes type;      /**< The measurement type. @see MeasurementTypes */
    MatX pose;                  /**< The part of the measurement corresponding to the pose of the target. (position, attitude, or both). */
    MatX twist;                 /**< The part of the measurement corresponding to the derivative of the pose. (velocity, angular rates, or both). */
    bool state_transform_data=false;  /**< Sometimes the measurement cannot be transformed into the tracking frame, but the tracks can be transformed into the measurement frame. If the tracks
                                     need to be transformed into the measurement frame, set this value to true. */
    MatX trans_data;            /**< The data used by the transform manager to transform the track from the tracking frame to the measurement frame. */                                


    // These member variables are reserved
    double probability=0;          /**< The probability that the measurement came from the track it was associated with. This value is set during the data
                                      association process.*/
    double weight=0;              /**< The weight assigned to the measurement when updating the track. This value is set during the data association process. */
};

//////// Measurement Type Utilities 






} // namespace rransac

#endif // RRANSAC_COMMON_MEASUREMENT_MEASUREMENT_BASE_H_
