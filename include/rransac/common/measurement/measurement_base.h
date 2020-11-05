#ifndef RRANSAC_COMMON_MEASUREMENT_MEASUREMENT_BASE_H_
#define RRANSAC_COMMON_MEASUREMENT_MEASUREMENT_BASE_H_

#include <Eigen/Core>
#include <string>

#include "rransac/parameters.h"

namespace rransac
{


/** \class MeasurementTypes
 * Lists the different types of measurements available.
 */ 
enum class MeasurementTypes {
    R2_POSE,                  // The Measurement space is R2 and with position data
    R2_POSE_TWIST,            // The Measurement space is R2 and with position and velocity data
    SE2_POSE,                 // The Measurement space is SE2
    SE2_POSE_TWIST            // The Measurement space is SE2 and se2
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

    double likelihood;          /**< The likelihood that the measurement came from the phenomenon it was associated with. This value is set during the data
                                      association process.*/
    double weight;              /**< The weight of the measurement when updating the model is was associated with. This value is set during the data association
                                     process. */

    Eigen::MatrixXd data;       /**< The data that represents the measurement. */
    Eigen::MatrixXd meas_cov;  /**< The measurement covariance. Only used if the measurement covariance changes with different measurements; otherwise, the
                                     measurement covariance given to the Source class is used for every measurement. */

    // const MeasurementTypes type;
    // Meas(MeasurementTypes type) : type{type} {}   

    // virtual ~MeasBase(){}
};





} // namespace rransac

#endif // RRANSAC_COMMON_MEASUREMENT_MEASUREMENT_BASE_H_
