#include "rransac/common/measurement/measurement_base.h"

/** \class Measurement
 * The default derived measurement class.
*/

class Measurement: public MeasurementBase
{
public:
    /**
     * Transforms a measurement from the previous global frame to the current global frame.
     * @param[in/out] meas The measurement that will be transformed.
     * @param[in] Transformation The transformation object that contains the data necessary to transform the measurement.
     * @param[in] dt The time elapsed since the last transformation was received.
     */
    static void TransformMeasurement(Meas& meas,const Transformation& T,const double dt);


    /**
     * Calculates the distance between two measurements depending on the type of distance to calculate.
     * @param[in] meas1 A measurement.
     * @param[in] meas2 A different measurement.
     * @param[in] type The type of distance to be calculated: spatial, temporal, or total 
     * @param[in] params Contains all of the user defined parameters. A user can define a weight when calculating the distances. 
     */
    static float GetDistance(const Meas& meas1, const Meas& meas2, const DistanceType& type, const Parameters& params);

private:

    /**
     * Calculates the spatial distance between two measurements using the 2-norm.
     * @param[in] meas1 A measurement.
     * @param[in] meas2 A different measurement.
     * @param[in] params Contains all of the user defined parameters. A user can define a weight when calculating the distances.
     */
    static float GetSpatialDistance(Meas& meas1, Meas& meas2, const Parameters& params);

    /**
     * Calculates the temporal distance between two measurements using the 2-norm.
     * @param[in] meas1 A measurement.
     * @param[in] meas2 A different measurement.
     * @param[in] params Contains all of the user defined parameters. A user can define a weight when calculating the distances.
     */
    static float GetTemporalDistance(Meas& meas1, Meas& meas2, const Parameters& params);

    /**
     * Calculates the temporal and spatial distance between two measurements using the 2-norm.
     * @param[in] meas1 A measurement.
     * @param[in] meas2 A different measurement.
     * @param[in] params Contains all of the user defined parameters. A user can define a weight when calculating the distances.
     */
    static float GetTotalDistance(Meas& meas1, Meas& meas2, const Parameters& params);
 

};