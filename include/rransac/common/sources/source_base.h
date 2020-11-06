#ifndef RRANSAC_COMMON_SOURCES_SOURCE_BASE_H_
#define RRANSAC_COMMON_SOURCES_SOURCE_BASE_H_

#include <Eigen/Core>
#include "state.h"
#include "common/measurement/measurement_base.h"
#include "parameters.h"

namespace rransac
{







/** \class SourceParameters
 * This struct contains the parameters needed by a measurement source.
 */ 
struct SourceParameters {
    bool meas_cov_fixed_;            /** < Flag used to indicate if the measurement covariance is the same for every measurement. 
                                      If it is, then the measurement covariance needs to be given to the source. If it isn't, 
                                      then the measurement covariance needs to be given with the measurement. @see Meas */
    Eigen::MatrixXd meas_cov_;       /** < The fixed measurement covariance.*/
    float expected_num_false_meas_;  /** < The expected number of false measurements. We assume that the expected number of false 
                                        measurements from a source per sensor scan can be modeled using a Poisson distribution. */

    MeasurementTypes type_;     /** < The source type @see SourceTypes */

    float probability_of_detection_; /**< The probability that the phenomenon of interest is detected by a source during
                                      a single scan. This value must be between 0 and 1.*/

    
    unsigned int source_index_;  /**< When a new source is added, it is added to the vector System::sources_. This is used to verify that the measurement corresponds to the proper source. */


    // SourceParameters(MeasurementTypes type) : type_(type) {}
    SourceParameters()=default;

};

/** \class SourceBase 
 * A source is an algorithm that takes data from a sensor and extracts measurements from it. Due to the nature of
 * the source, there is measurement noise and false measurements. The measurement noise is assumed to be sampled from
 * a white-noise, zero-mean, Gaussian distribution. If the measurement covariance is fixed, then it should be added to 
 * the source; otherwise, it is supplied with each measurement. We assume that the false measurements are uniformly
 * distributed in the local surveillance region, and that the spatial density of false measurements can 
 * be modeled using a Poisson distribution. 
 * 
 * Since there can be many different sources, we allow users to create their own source type using polymorphism. The 
 * base source must contain all of the necessary member functions and the child classes provide the implementation for
 * specific sources. 
 * 
 * When a measurement is received, the measurement indicates which source it came from using the source ID tag @see Meas.
 * We also need to know the dimension of the measurement.
 * 
 */ 


class SourceBase
{

public:

    SourceParameters params_;  /** < The source parameters @see SourceParameters */

    /** Initializes the measurement source */
    template<MeasurementTypes type>
    inline void Init(const SourceParameters& params);       

    /** Returns the jacobian of the observation function w.r.t. the states */
    template <MeasurementTypes type, class S>
    inline Eigen::MatrixXd GetLinObsMatState(const S& state);                              

    /** Returns the jacobian of the observation function w.r.t. the sensor noise */
    template <MeasurementTypes type, class S>
    inline Eigen::MatrixXd GetLinObsMatSensorNoise(const S& state);                         

    /** Computes the estimated measurement given a state */
    template <MeasurementTypes type, class S>
    inline Eigen::MatrixXd GetEstMeas(const S& state); /** Returns an estimated measurement according to the state. */

    /**
     * Calculates the temporal distance between two measurements.
     * @param[in] meas1 A measurement.
     * @param[in] meas2 A different measurement.
     * @param[in] params Contains all of the user defined parameters. A user can define a weight when calculating the distances.
     * \return Returns temporal distance between two measurements
     */
   
    static float GetTemporalDistance(const Meas& meas1, const Meas& meas2, const Parameters& params) { return fabs(meas1.time_stamp - meas1.time_stamp); }

    /**
     * Calculates the spatial distance between two measurements depending on the type of measurement.
     * @param[in] meas1 A measurement.
     * @param[in] meas2 A different measurement.
     * @param[in] params Contains all of the user defined parameters. A user can define a weight when calculating the distances.
     * \return Returns spatial distance between two measurements
     */
    template<MeasurementTypes M1, MeasurementTypes M2>
    static float GetSpatialDistance(const Meas& meas1, const Meas& meas2, const Parameters& params) {throw std::runtime_error("SourceBase: Spacial Distance is not implemented for the given measurements.");}


private:
    Eigen::MatrixXd H_;
    Eigen::MatrixXd V_;

};




} // namespace rransac

#include "common/sources/source_R2_pos.h"
#include "common/sources/source_R2_pos_vel.h"

#endif // RRANSAC_COMMON_SOURCES_SOURCE_BASE_H_