#ifndef RRANSAC_DATA_STRUCTURES_CONSENSUS_SET_H_
#define RRANSAC_DATA_STRUCTURES_CONSENSUS_SET_H_
#pragma once

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

template <class tMeasurement>
class ConsensusSet
{
public:

/** 
 * Add a measurement to the consensus set. 
 * @param[in] meas The measurement to be added.
*/
void AddMeasToConsensusSet(const tMeasurement& meas);

/** 
 * Add a measurements to the consensus set. 
 * @param[in] meas The measurements to be added.
*/
void AddMeasurementsToConsensusSet(const std::vector<tMeasurement>& meas);

/** 
 * Add a measurements to the consensus set. The measurements have the same time stamp.
 * @param[in] meas The measurements to be added.
*/
void AddMeasurementsToConsensusSetSameTimeStamp(const std::vector<tMeasurement>& meas);

/** 
 * Removes all of the measurements from the consensus_set with a time stamp that occurred before the expiration_time. 
 * The value of expiration_time is the current time in seconds minus the time window. I.e. old measurements are removed.
 * @param[in] expiration_time The current time in seconds minus the time window.
 */ 

void PruneConsensusSet(double expiration_time);

/**
 * Transforms all of the measurements in the consensus set using the provided transformation.
 */ 
template<typename tTransformation>
void TransformConsensusSet(const tTransformation& T);

/**
 * Returns the size of the consensus set.
 */ 
unsigned int Size() { return consensus_set_.size();}

/*
* Merges two consensus sets together. This is used when two models are merged together.
*/
static ConsensusSet MergeConsensusSets(const ConsensusSet& cs1, const ConsensusSet& cs2);


std::list<std::vector<tMeasurement>> consensus_set_; /** < Contains the measurements associated with the model that have not expired. Each vector of measurements 
                                                contains measurements with the same time stamp. This allows us to efficiently remove expired measurements
                                                by simply removing an entire vector. */

};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                            Definitions
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <class tMeasurement>
void ConsensusSet<tMeasurement>::AddMeasToConsensusSet(const tMeasurement& meas) {

    // There are no measurements. Just add it.
    if (consensus_set_.size() == 0) {
        consensus_set_.emplace_back(std::vector<tMeasurement>{meas});
    } 
    // All measurements occurred before the new one. So add it to the back.
    else if (consensus_set_.back().front().time_stamp < meas.time_stamp) {
        consensus_set_.emplace_back(std::vector<tMeasurement>{meas});
    } 
    // The new measurement occurred before all the other measurements
    else if (consensus_set_.front().front().time_stamp > meas.time_stamp) {
        consensus_set_.emplace_front(std::vector<tMeasurement>{meas});
    }
    else {
        // Search from the end to the beginning to find where to place it
        for (auto iter = consensus_set_.rbegin(); iter != consensus_set_.rend(); ++iter) {
            
            // std::cerr << (*iter).front().time_stamp << std::endl;

            if ((*iter).front().time_stamp == meas.time_stamp) {
                (*iter).push_back(meas);
                break;
            } else if ((*iter).front().time_stamp < meas.time_stamp) {
                // std::cout << "add to middle" << std::endl;
                // std::vector<M> tmp{meas};
                consensus_set_.insert(iter.base(),std::vector<tMeasurement>{meas});
                break;
            } 
        }
    }
}

//-----------------------------------------------------------------

template <class tMeasurement>
void ConsensusSet<tMeasurement>::AddMeasurementsToConsensusSet(const std::vector<tMeasurement>& meas) {

    for (tMeasurement m : meas) {
        AddMeasToConsensusSet(m);
    }

}

//-----------------------------------------------------------------

template <class tMeasurement>
void ConsensusSet<tMeasurement>::AddMeasurementsToConsensusSetSameTimeStamp(const std::vector<tMeasurement>& meas) {


 // There are no measurements. Just add it.
    if (consensus_set_.size() == 0) {
        consensus_set_.emplace_back(meas);
    } 
    // All measurements occurred before the new one. So add it to the back.
    else if (consensus_set_.back().front().time_stamp < meas[0].time_stamp) {
        consensus_set_.emplace_back(meas);
    } 
    // The new measurement occurred before all the other measurements
    else if (consensus_set_.front().front().time_stamp > meas[0].time_stamp) {
        consensus_set_.emplace_front(meas);
    }
    else {
        // Search from the end to the beginning to find where to place it
        for (auto iter = consensus_set_.rbegin(); iter != consensus_set_.rend(); ++iter) {
            
            // std::cerr << (*iter).front().time_stamp << std::endl;

            if ((*iter).front().time_stamp == meas[0].time_stamp) {
                iter->insert(iter->end(), meas.begin(), meas.end());
                break;
            } else if ((*iter).front().time_stamp < meas[0].time_stamp) {
                // std::cout << "add to middle" << std::endl;
                // std::vector<M> tmp{meas};
                consensus_set_.insert(iter.base(),meas);
                break;
            } 
        }
    }


}

//-----------------------------------------------------------------

template <class tMeasurement>
void ConsensusSet<tMeasurement>::PruneConsensusSet(double expiration_time) {

// If all of the measurements are expired, clear the list
if(consensus_set_.back().front().time_stamp < expiration_time) {
    consensus_set_.clear();
    return;
} 
// If none of the measurements are expired, don't remove any
if ( consensus_set_.front().front().time_stamp > expiration_time) {
    return;
}

// Some measurements are not expired, so remove only the ones that are
for (auto iter = consensus_set_.begin(); iter != consensus_set_.end(); ++iter) {
    if (iter->front().time_stamp > expiration_time) {
        consensus_set_.erase(consensus_set_.begin(), iter);
        break;
    }

}
}

//-----------------------------------------------------------------

template< class tMeasurement>
template<typename tTransformation>
void ConsensusSet<tMeasurement>::TransformConsensusSet(const tTransformation& T) {

    for (auto iter = consensus_set_.begin(); iter != consensus_set_.end(); ++iter) {
        for (auto& m : (*iter) ) {
            T.TransformMeasurement(m);
        }
    }
}

//-----------------------------------------------------------------

template<class tMeasurement>
ConsensusSet<tMeasurement> ConsensusSet<tMeasurement>::MergeConsensusSets(const ConsensusSet& cs1, const ConsensusSet& cs2) {

    ConsensusSet<tMeasurement> merged_consensus_set = cs1;
    
    for (auto iter = cs2.consensus_set_.begin(); iter != cs2.consensus_set_.end(); ++iter) {
        
        merged_consensus_set.AddMeasurementsToConsensusSetSameTimeStamp((*iter));

    }
    
    return merged_consensus_set;

}

} // namespace rransac
#endif // RRANSAC_DATA_STRUCTURES_CONSENSUS_SET_H_