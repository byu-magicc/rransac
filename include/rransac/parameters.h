#ifndef RRANSAC_COMMON_PARAMETERS_H_
#define RRANSAC_COMMON_PARAMETERS_H_

#include <vector>
#include <Eigen/Core>


namespace rransac
{
/** \class Parameters
 * This class manages all of the parameters used in R-RANSAC
 * that are set by the user.
*/
class Parameters
{
public:

/**
 * Sets the values of the member variables to acceptable values no prevent
 * numerical errors. These values should be changed by the user via
 * SetParameters(const Parameters &new_params).
 */
Parameters();

/**
 * The copy constructor. It simply calls SetParameters(const Parameters &new_params).
 * \param[in] new_params The object that contains the new parameters.
 */
Parameters(const Parameters &new_params);

void operator=(const Parameters &new_params);

// Deconstructor
~Parameters();

/**
 * \detail Sets all of the parameters
 * \param[in] new_params The new parameters.
 * \return If the parameters were successfully set (i.e. met all of the prescibed requirements), it returns true. Otherwise false.
 */
bool SetParameters(const Parameters &new_params);

float meas_time_window_; /**< The duration of time in seconds from the current moment and extending
                              into the past during which measurements should be considered. All
                              measurements outside of the time window will be discarded.*/

bool fixed_time_interval_; /**< A flag that indicates if measurements are given to R-RANSAC at fixed time intervals. */

double time_interval_; /**< The fixed time interval at which measurements are received.*/

bool transform_consensus_set_; /**< A flag that indicates if the measurements in the consensus set should
                                    be transformed when a transformation is provided. For faster performance, it is
                                    recommended that this flag is set to false.*/

// Cluster Parameters
double cluster_time_threshold_;       /**< In order for a measurements to be a neighbor to another measurement of a different time stamp, the difference in time
                                           must be less than or equal to  cluster_time_threshold_*/
double cluster_velocity_threshold_;   /**< In order for a measurement to be a neighbor to another measurement of a different time stamp, the distance in pose 
                                           normalized by the time difference between the two measurements must be less than or equal to the cluster_velocity_threshold_ */
double cluster_position_threshold_;   /**< In order fo a measurement to be a neighbor to another measurement of the same time stamp, the distance in pose
                                            must be less than or equal to the cluster_position_threshold_ */


// RANSAC Parameters
unsigned int RANSAC_max_iters_;  /**< The maximum number of RANSAC iterations per run. */

unsigned int RANSAC_score_stopping_criteria_; /**< During any iteration, if the score of a hypothetical state is greater than this, then RANSAC terminates early 
                                                    and creates a new track using the hypothetical state*/
unsigned int RANSAC_score_minimum_requirement_; /**< A hypothetical state must have a score greater than or equal to this in order to be made into a track. */


unsigned int RANSAC_minimum_subset_;    /** The minimum number of measurements from different time steps required to observe the system */



// Model Parameters
Eigen::MatrixXd process_noise_covariance_; /**< The process noise covariance of the model */




// Model Manager
double good_model_threshold_;
double max_missed_detection_time_; /** < A model that has  not received an associated measurement in this time (s) will be discarded */
double similar_tracks_threshold_;
int max_num_models_;               /** The maximum number of models that RRANSAC will store */



};
} // namespace rransac

#endif // RRANSAC_COMMON_PARAMETERS_H_
