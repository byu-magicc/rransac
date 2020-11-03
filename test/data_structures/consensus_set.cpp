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
        consensus_set_.emplace_back(std::vector<M>{meas})
    } 
    else {
        // Search from the end to the beginning to find where to place it
        for (auto iter = consensus_set_.end(); iter != consensus_set_.begin(); --iter) {
            

            if (*iter.front().time_stamp == meas.time_stamp) {
                *iter.push_back(meas);
                break;
            } else if (*iter.front().time_stamp < meas.time_stamp) {
                consensus_set_.insert(iter+1,std::vector<M>{Meas});
                break;
            }

            }

        }
    }

    }

    // See if there are existing measurements with the same time stamp
  
}

} // namespace rransac


