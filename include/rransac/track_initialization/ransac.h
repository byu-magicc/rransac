#ifndef RRANSAC_TRACK_INITIALIZATION_RANSAC_H_
#define RRANSAC_TRACK_INITIALIZATION_RANSAC_H_

#include "system.h"
#include "data_containers/cluster.h"
#include <stdlib.h>
#include <time.h>
#include <algorithm> 
#include <numeric>

namespace rransac {

template< typename tLMLEPolicy>
class Ransac : public tLMLEPolicy {

public: 

Ransac();

~Ransac() = default;

/**
 * Randomly selects a measurement from the latest time step and randomly selects 
 * measurements from other time steps until there are num_meas randomly selected.
 * @param num_meas The total number of random numbers to select. 
 * @param cluster The cluster from which to the measurements will be sampled. 
 * @return Returns the indicess of the randomply sampled measurements 
 */ 
static std::vector<Cluster::IteratorPair> GenerateMinimumSubset(const unsigned int num_meas, const Cluster& cluster);

private:


};

template< typename tLMLEPolicy>
Ransac<tLMLEPolicy>::Ransac() {
    srand(time(NULL));
}

//----------------------------------------------------------------------------------------------------------

template< typename tLMLEPolicy>
std::vector<Cluster::IteratorPair> Ransac<tLMLEPolicy>::GenerateMinimumSubset(const unsigned int num_meas, const Cluster& cluster) {
   


    std::vector<Cluster::IteratorPair> meas_index(num_meas);

    // Get a random measurement from current time step
    meas_index.back().outer_it = std::prev(cluster.data_.end());
    int length = cluster.data_.back().size();
    int index = rand() % length;
    meas_index.back().inner_it = std::next(cluster.data_.back().begin(),index);

    if (num_meas == 1)
        return meas_index;

    // Get the time indices for the other randomly sampled measurements
    std::vector<int> time_indices(num_meas-1);
    std::iota(time_indices.begin(), time_indices.end(), 0);
    std::shuffle(time_indices.begin(), time_indices.end(), rand() );

   

    std::list<Meas>::iterator outer_iter;

    // Randomly samples other measurements
    for (unsigned int ii = 0; ii < num_meas-1; ++ii) {

        meas_index[ii].outer_it = std::next(cluster.data_.begin(), time_indices[ii]);
        meas_index[ii].inner_it = rand() % meas_index[ii].outer_it->size();


    }

    return meas_index;
    
}


} //namespace rransac

#endif // RRANSAC_TRACK_INITIALIZATION_RANSAC_H_