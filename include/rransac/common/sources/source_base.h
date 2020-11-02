#ifndef RRANSAC_COMMON_SOURCES_SOURCE_BASE_H_
#define RRANSAC_COMMON_SOURCES_SOURCE_BASE_H_

#include <Eigen/Core>
#include "state.h"

namespace rransac
{

/** \class Source 
 * This struct contains the parameters needed by a measurement source.
 */ 
struct SourceParameters {
    bool meas_cov_fixed_;            /** < Flag used to indicate if the measurement covariance is the same for every measurement. 
                                      If it is, then the measurement covariance needs to be given to the source. If it isn't, 
                                      then the measurement covariance needs to be given with the measurement. @see Meas */
    Eigen::MatrixXd meas_cov_;       /** < The fixed measurement covariance.*/
    float expected_num_false_meas_;  /** < The expected number of false measurements. We assume that the expected number of false 
                                        measurements from a source per sensor scan can be modeled using a Poisson distribution. */
};

/** \class Source 
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

// Lists the different types of sources available
enum class SourceTypes {
    R2_POS,                // The Target's manifold is R2 and only position is observed
    R2_POS_VEL,            // The Target's manifold is R2 and position and velocity is observed
    SE2_POS,               // The Target's manifold is SE2 and only position is observed
    SE2_POS_ATT            // The target's manifold is SE2 and the position and attitude is observed
};

template <class G, class U>
class SourceBase
{

public:

    SourceParameters params_;  /** < The source parameters @see SourceParameters */

    unsigned int source_id_;    /** < An identifier that is unique for every source. When a measurement is given, the measurement
                                      indicates which source it came from using the source ID. @see Meas */

    // unsigned int meas_dim_;     /** < The dimension of the measurement. You can also consider this as the number of generalized coordinates measured.
    //                               For example, if the measurement containted the position of a 2D object in the
    //                               xy-plane, then the dimension of the measurement would be two. This is used by RANSAC when creating model hypotheses. */

    virtual void Init(const SourceParameters& params, unsigned int source_id)=0; /** Initializes the measurement source */
    virtual Eigen::MatrixXd GetLinObsMatState()=0;               /** Returns the jacobian of the observation function w.r.t. the states */
    virtual Eigen::MatrixXd GetLinObsMatSensorNoise()=0;         /** Returns the jacobian of the observation function w.r.t. the sensor noise */
    virtual Eigen::MatrixXd GetEstMeas(const lie_groups::State<G,U>& state)=0; /** Returns an estimated measurement according to the state. */
    

};
} // namespace rransac

#endif // RRANSAC_COMMON_SOURCES_SOURCE_BASE_H_