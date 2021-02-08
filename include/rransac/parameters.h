#ifndef RRANSAC_COMMON_PARAMETERS_H_
#define RRANSAC_COMMON_PARAMETERS_H_
#pragma once

#include <vector>
#include <Eigen/Core>


namespace rransac
{
/** \class Parameters
 * This class manages and contains all of the system parameters in RRANSAC except for source specific parameters.
*/
class Parameters
{
public:

/**
 * Sets the values of the member variables to acceptable values.
 * The parameter values can be changed by the member function SetParameters()
 */
Parameters();

/**
 * The copy constructor. It Calls SetParameters().
 * \param[in] new_params The object that contains the new parameters.
 */
Parameters(const Parameters &new_params);

/**
 * The assignment operator. It Calls SetParameters().
 * \param[in] new_params The object that contains the new parameters.
 */
void operator=(const Parameters &new_params);

// Deconstructor
~Parameters();

/**
 * \detail Verifies that all of the new parameters are withing acceptable ranges. If they are, then 
 * the current parameters will be changed to the new parameters.
 * \param[in] new_params The new parameters.
 * \return If the parameters were successfully set (i.e. met all of the prescibed requirements), it returns true. Otherwise false.
 */
bool SetParameters(const Parameters &new_params);

float meas_time_window_; /**< The duration of time in seconds from the current time stamp extending
                              into the past during which measurements should be considered. All
                              measurements outside of the time window will be discarded. For example: a value of 
                              5 means that all the measurement window is from the current time stamp to 
                              curent time stamp -5 seconds. */


bool transform_consensus_set_; /**< A flag that indicates if the measurements in the consensus set should
                                    be transformed when a transformation is provided. For faster performance, it is
                                    recommended that this flag is set to false. */

// Cluster Parameters
double cluster_time_threshold_;             /**< In order for a measurements to be a neighbor to another measurement of a different time stamp, the difference in time
                                              must be less than or equal to  this parameter. */
double cluster_velocity_threshold_;         /**< In order for a measurement to be a neighbor to another measurement of a different time stamp, the distance in pose 
                                              normalized by the time difference between the two measurements must be less than or equal to this parameter. */
double cluster_position_threshold_;         /**< In order for a measurement to be a neighbor to another measurement of the same time stamp, the distance in pose
                                              must be less than or equal to this parameter. */
unsigned int cluster_min_size_requirement_; /**< A cluster must have at least this many measurements before the cluster will be used to initialize a track. */


// RANSAC Parameters
unsigned int RANSAC_max_iters_;  /**< The maximum number of RANSAC iterations per run. */

unsigned int RANSAC_score_stopping_criteria_; /**< During any iteration, if the score of a hypothetical state is greater than this, then RANSAC terminates early 
                                                    and creates a new track using the hypothetical state. */
unsigned int RANSAC_score_minimum_requirement_; /**< A hypothetical state must have a score greater than or equal to this in order to be made into a track. */


unsigned int RANSAC_minimum_subset_;    /**< The minimum number of measurements from different time stamps used to generate the hypothetical state estimate. This value must be
                                            greater than or equal to the number of measurements needed to observer the target. */



// Track Parameters
Eigen::MatrixXd process_noise_covariance_; /**< The process noise covariance of the track. It must be positive definite. */




// Model Manager
double track_good_model_threshold_;      /**< A track whose ModelBase::model_likelihood_ is equal to or greater than this parameter, is considered a 
                                               good model; otherwise, it is considered a poor model. */
double track_max_missed_detection_time_; /**< A track that has not received an associated measurement in the past time in seconds indicated by this parameter
                                               is assumed to be outside the surveillance region of any sensor and will be discarded. */
double track_similar_tracks_threshold_;  /**< The threshold used to determine if two tracks are similar enough to be considered the same. If they are, then they will 
                                               be merged. */
int track_max_num_tracks_;               /**< The maximum number of tracks that RRANSAC will store. */

// Nonlinear LMLE 
bool nonlinear_innov_cov_id_;                      /**< When set to true, the innovation covariance will not be calculated in the nonlinear track initialization process.
                                                    This can drastically speed up performance at the cost of accuracy. */
unsigned int nonlinear_LMLE_Ceres_threads_;        /**< The number of threads the nonlinear track initializer will use. */
unsigned int nonlinear_LMLE_Ceres_max_num_iters_;  /**< When solving the LMLE, this parameter indicates the maximum number of iterations Ceres will use
                                                         to estimate the current state estimate of a hypothetical track. */

};
} // namespace rransac

#endif // RRANSAC_COMMON_PARAMETERS_H_
