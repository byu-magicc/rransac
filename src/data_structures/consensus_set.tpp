#include "rransac/data_structures/consensus_set.h"


namespace rransac
{

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






} // namespace rransac


