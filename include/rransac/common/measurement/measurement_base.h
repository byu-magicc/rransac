#ifndef RRANSAC_COMMON_MEASUREMENT_MEASUREMENT_BASE_H_
#define RRANSAC_COMMON_MEASUREMENT_MEASUREMENT_BASE_H_

#include <Eigen/Core>

#include "rransac/common/transformation.h"
#include "rransac/common/param.h"

namespace rransac
{
/** \struct Meas
 * This struct contains information regarding a measurement. The user uses this object to supply R-RANSAC with measurements.
 * The user is responsible to provide Meas::data, Meas::time_stamp and Meas::source_id. If the measurement covariance is not provided
 * to the Source object, the user must also provide Meas::meas_cov. The other member variables are set by R-RANSAC.
*/
struct Meas
{
    Eigen::MatrixXd data;       /**< The data that represents the measurement. */
    double time_stamp;          /**< The time the measurement was taken. */
    unsigned int source_id;     /**< A unique identifier that indicates which source the measurement came from. */
    Eigen::MatrixXd meas_cov;   /**< The measurement covariance. Only used if the measurement covariance changes with different measurements; otherwise, the
                                     measurement covariance given to the Source class is used for every measurement. */
    double likelihood;          /**< The likelihood that the measurement came from the phenomenon it was associated with. This value is set during the data
                                      association process.*/
    double weight;              /**< The weight of the measurement when updating the model is was associated with. This value is set during the data association
                                     process. */
};


/** \enum Meas
 * Used to indicate what type of distance between measurements to calculated.
*/
enum DistanceType
{
    kSpatial=0,
    kTemporal=1,
    kTotal=2
};

/** \enum DerivedMeasurement
 * Used to indicate what type of derived measurement class to create.
*/
enum DerivedMeasurement
{
    kMeasurement=0
}


/** \class MeasurementBase
 * The class is responsible for performing any calculation that involves only measurements.
 * All of the member functions in this class should be static.
*/

class MeasurementBase
{
public:

  MeasurementBase();
  ~MeasurementBase();

  /**
   * Transforms a measurement from the previous global frame to the current global frame.
   * @param[in/out] meas The measurement that will be transformed.
   * @param[in] Transformation The transformation object that contains the data necessary to transform the measurement.
   * @param[in] dt The time elapsed since the last transformation was received.
   */
  virtual void TransformMeasurement(Meas& meas,const Transformation& T,const double dt)=0;

  /**
   * Calculates the distance between two measurements depending on the type of distance to calculate.
   * @param[in] meas1 A measurement.
   * @param[in] meas2 A different measurement.
   * @param[in] type The type of distance to be calculated: spatial, temporal, etc.
   * @param[in] params Contains all of the user defined parameters. A user can define a weight when calculating the distances.
   */
  virtual float GetDistance(const Meas& meas1, const Meas& meas2, const DistanceType& type, const Parameters& params)=0;

};
} // namespace rransac

#endif // RRANSAC_COMMON_MEASUREMENT_MEASUREMENT_BASE_H_
