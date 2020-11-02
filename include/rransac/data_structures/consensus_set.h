#ifndef RRANSAC_DATA_STRUCTURES_CONSENSUS_SET_H_
#define RRANSAC_DATA_STRUCTURES_CONSENSUS_SET_H_

#include "rransac/common/measurement/measurement_base.h"
#include <list>

namespace rransac
{
/**
 * \class ConsensusSet
 * The consensus set of a model is the set of measurements associated with the model.
 * In order to keep the consensus set becoming arbitrarily large. Expired measurements are 
 * pruned from the consensus set. The data member consensus_set is a list of vectors that contain
 * the measurements. Each vector of measurements contains measurements with the same time stamp. 
 */ 

template <class M>
class ConsensusSet
{
public:

/** < 
 * Add a measurement to the consensus set. 
 * @param[in] meas The measurement to be added.
*/
void AddMeasToConsensusSet(const M& meas);

/** < 
 * Add a measurements to the consensus set. 
 * @param[in] meas The measurement to be added.
*/
void AddMeasurementsToConsensusSet(const std::vector<M>& meas);

/** <
 * Removes all of the measurements from the consensus_set with a time stamp that occurred before expiration_time. 
 * The value of expiration_time is the current time in seconds minus the time window. I.e. old measurements are removed.
 * @param[in] expiration_time The current time in seconds minus the time window.
 */ 

void PruneConsensusSet(double expiration_time);


std::list<std::vector<M>> consensus_set; /** < Contains the measurements associated with the model that have not expired. Each vector of measurements 
                                                contains measurements with the same time stamp. This allows us to efficiently remove expired measurements
                                                by simply removing an entire vector. */
};
} // namespace rransac

#endif // RRANSAC_DATA_STRUCTURES_CONSENSUS_SET_H_