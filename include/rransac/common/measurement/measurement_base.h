#ifndef RRANSAC_COMMON_MEASUREMENT_MEASUREMENT_BASE_H_
#define RRANSAC_COMMON_MEASUREMENT_MEASUREMENT_BASE_H_

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
    SEN_POSE,                 // The target space is SEN and the pose is measured and is an element of SEN
    SEN_POSE_TWIST,           // The target space is SEN and the pose and twist is measured and are elements of SEN and seN
    SEN_POS,                  // The target space is SEN and the position is measured
    SEN_POS_VEL,              // The target space is SEN and the position and velocity is measured
    NUM_TYPES
    // SO3_ATT
};


/** \struct Meas
 * This struct contains information regarding a measurement. The user uses this object to supply R-RANSAC with measurements.
 * The user is responsible to provide Meas::data, Meas::time_stamp and Meas::source_id. If the measurement covariance is not provided
 * to the Source object, the user must also provide Meas::meas_cov. The other member variables are set by R-RANSAC.
*/

struct Meas
{
    
    double time_stamp;          /**< The time the measurement was taken. */
    unsigned int source_index;  /**< When a new source is added, it is added to the vector System::sources_. The source indexes the vector to grab the corresponding source. So make sure it's the right index. */
    MeasurementTypes type;      /** < The measurement type @see MeasurementTypes */
    Eigen::MatrixXd pose;       /**< The part of the measurement corresponding to the pose of the target. (position, attitude, or both). */
    Eigen::MatrixXd twist;      /**< The part of the measurement corresponding to the derivative of the pose. (velocity, angular rates, or both). */
    Eigen::MatrixXd meas_cov;   /**< The measurement covariance. Only used if the measurement covariance changes with different measurements; otherwise, the
                                     measurement covariance given to the Source class is used for every measurement. */
  
    // These member variables are reserved
    double likelihood;          /**< The likelihood that the measurement came from the phenomenon it was associated with. This value is set during the data
                                      association process.*/
    double weight;              /**< The weight of the measurement when updating the model is was associated with. This value is set during the data association
                                     process. */
    Eigen::MatrixXd pose_euclidean; /**< The r-star tree uses bounding boxes in Euclidean space to organize the measurements. Some pose measurements are not
                                         in Euclidean space and need to be mapped to a Euclidean space. This variable contains that counterpart and can be set 
                                         using the source_base. */
};


} // namespace rransac

#endif // RRANSAC_COMMON_MEASUREMENT_MEASUREMENT_BASE_H_
