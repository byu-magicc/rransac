#ifndef RRANSAC_COMMON_PARAMETERS_H_
#define RRANSAC_COMMON_PARAMETERS_H_

#include <vector>

#include "rransac/common/Source.h"

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

  // Deconstructor
  ~Parameters();

  /**
   * \detail Sets the values of the member variables of this object to the
   * values of the member variables of the passed in object. Some values are checked to ensure that
   * they meet the prescribed requirements.
   * \param[in] new_params The object that contains the new parameters.
   * \return If the parameters were successfully set (i.e. met all of the prescibed requirements), it returns true. Otherwise false.
   */
  bool SetParameters(const Parameters &new_params);


  bool fixed_time_interval_; /**< A flag that indicates if measurements are given to R-RANSAC at fixed time intervals. */

  double time_interval_; /**< The fixed time interval at which measurements are received.*/

  bool transform_consensus_set_; /**< A flag that indicates if the measurements in the consensus set should
                                      be transformed when a transformation is provided. For faster performance, it is
                                      recommended that this flag is set to false.*/

  float meas_time_window_; /**< The duration of time in seconds from the current moment and extending
                                into the past during which measurements should be considered. All
                                measurements outside of the time window will be discarded.*/

  // JPDAF parameters
  float probability_of_detection_; /**< The probability that the phenomenon of interest is detected by a sensor during
                                        a single scan. This value must be between 0 and 1.*/



  // TODO:: Remove this
  float expected_num_false_meas_; /**< The expected number of false measurements inside the entire surveillance region
                     per measurement scan.*/

  // RANSAC Parameters
  unsigned int max_RANSAC_iters_;  /**< The maximum number of RANSAC iterations per run. */

  float RANSAC_stopping_criteria_; /**< During any iteration, if the probability of a model hypothesis
                                        a valid model is above this threshold, RANSAC stops early and
                                        uses the current model hypothesis to generate a new model.
                                        This value must be between 0 and 1.*/


  std::vector<Source> sources_; /** < Contains all of the instantiated sources. */

  /**
   * This method adds a source to the vector Parameters::sources_. Before adding a source,
   * it verifies that the source has a unique ID and that the member variables of the source
   * are set properly.
   */
  void AddSource(Source& src);

  // TODO:: Write the definition of the function AddSource().

  };
} // namespace rransac

#endif // RRANSAC_COMMON_PARAMETERS_H_
