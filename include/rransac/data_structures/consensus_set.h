#include "rransac/common/measurement/measurement.h

/**
 * \class ConsensusSet
 * The consensus set of a model is the set of measurements associated with the model.
 * In order to keep the consensus set becoming arbitrarily large. Expired measurements are 
 * pruned from the consensus set. The data member consensus_set is a list of vectors that contain
 * the measurements. Each vector of measurements contains measurements with the same time stamp. 
 */ 

class ConsensusSet
{
public:

/** < 
 * Add a measurement to the consensus set. 
 * @param[in] meas The measurement to be added.
*/
void AddMeasToConsensusSet(Meas meas);

/** <
 * Removes all of the measurements from the consensus_set with a time stamp that occurred before expiration_time. 
 * The value of expiration_time is the current time in seconds minus the time window. I.e. old measurements are removed.
 * @param[in] expiration_time The current time in seconds minus the time window.
 */ 

void PruneConsensusSet(double expiration_time);


std::list<std::vector<Meas>> consensus_set; /** < Contains the measurements associated with the model that have not expired. Each vector of measurements 
                                                contains measurements with the same time stamp. This allows us to efficiently remove expired measurements
                                                by simply removing an entire vector. */
};
