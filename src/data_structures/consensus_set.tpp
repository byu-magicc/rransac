#include "rransac/data_structures/consensus_set.h"


namespace rransac
{

template <class M>
void ConsensusSet<M>::AddMeasToConsensusSet(const M& meas) {

    // There are no measurements. Just add it.
    if (consensus_set_.size() == 0) {
        consensus_set_.emplace_back(std::vector<M>{meas});
    } 
    // All measurements occurred before the new one. So add it to the back.
    else if (consensus_set_.back().front().time_stamp < meas.time_stamp) {
        consensus_set_.emplace_back(std::vector<M>{meas});
    } 
    // The new measurement occurred before all the other measurements
    else if (consensus_set_.front().front().time_stamp > meas.time_stamp) {
        consensus_set_.emplace_front(std::vector<M>{meas});
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
                consensus_set_.insert(iter.base(),std::vector<M>{meas});
                break;
            } 
        }
    }
}

//-----------------------------------------------------------------

template <class M>
void ConsensusSet<M>::AddMeasurementsToConsensusSet(const std::vector<M>& meas) {

    for (M m : meas) {
        AddMeasToConsensusSet(m);
    }

}

//-----------------------------------------------------------------
template <class M>
void ConsensusSet<M>::AddMeasurementsToConsensusSetSameTimeStamp(const std::vector<M>& meas) {


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
template <class M>
void ConsensusSet<M>::PruneConsensusSet(double expiration_time) {

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



} // namespace rransac


