#ifndef RRANSAC_COMMON_SOURCE_H_
#define RRANSAC_COMMON_SOURCE_H_

#include <Eigen/Core>

namespace rransac
{
/** \class Source 
 * A source is an algorithm that takes data from a sensor and extracts measurements from it. Due to the nature of
 * the source, there is measurement noise and false measurements. The measurement noise is assumed to be sampled from
 * a white-noise, zero-mean, Gaussian distribution. If the measurement covariance is fixed, then it should be added to 
 * the source; otherwise, it is supplied with each measurement. We assume that the false measurements are uniformly
 * distributed in the local surveillance region, and that the number
 * of expected measurements per sensor scan from a source can be modeled using a Poisson distribution. 
 * 
 * Since there can be many different sources, we need to be able to distinguish the different sources using a unique ID.
 * When a measurement is received, the measurement indicates which source it came from using the source ID tag @see Meas.
 * We also need to know the dimension of the measurement.
 * 
 */ 

class Source
{

public:

    bool meas_cov_fixed_;       /** < Flag used to indicate if the measurement covariance is the same for every measurement. 
                                      If it is, then the measurement covariance needs to be given to the source. If it isn't, 
                                      then the measurement covariance needs to be given with the measurement. @see Meas */
    Eigen::MatrixXd meas_cov_;  /** < The fixed measurement covariance.*/
    float expected_num_false_meas_;  /** < The expected number of false measurements. We assume that the expected number of false 
                                        measurements from a source per sensor scan can be modeled using a Poisson distribution. */
    unsigned int source_id_;    /** < An identifier that is unique for every source. When a measurement is given, the measurement
                                      indicates which source it came from using the source ID. @see Meas */

    unsigned int meas_dim_; /** < The dimension of the measurement. You can also consider this as the number of generalized coordinates measured.
                                  For example, if the measurement containted the position of a 2D object in the
                                  xy-plane, then the dimension of the measurement would be two. This is used by RANSAC when creating model hypotheses. */

    
};
} // namespace rransac

#endif // RRANSAC_COMMON_SOURCE_H_