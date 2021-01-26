#ifndef RRANSAC_TRACK_INITIALIZATION_SEED_POLOCIES_NULL_POLICY_H_
#define RRANSAC_TRACK_INITIALIZATION_SEED_POLOCIES_NULL_POLICY_H_

#include "system.h"
#include "data_containers/cluster.h"

namespace rransac
{

/**
 * \class NULLSeedPolicy
 * The seed policies are used to seed the LMLE algorithm that finds a current hypothetical state estimate using a minimum subset.
 * This policy is a null policy and sets the inital guess to identity.. 
 */ 

template<typename tModel>
class NULLSeedPolicy {

public:

typedef typename tModel::State State;
typedef typename State::DataType DataType;

/**
 * This function sets the initial guess of the state estimate to zero.
 * @param meas_subset The subset used to generate the seed for the LMLE algorithm
 * @param sys The object that contains the R-RANASAC data
 * @param x The initial state for the LMLE algorithm
 * @param size The size of the pointer x
 * 
 */ 
    static void GenerateSeed(const std::vector<typename Cluster<DataType>::IteratorPair>& meas_subset, const System<tModel>& sys, double x[tModel::cov_dim_], const int size) {
        for (int ii = 0; ii < size; ++ii) {
            x[ii]=0;
        }
    }

};



} // namespace rransac


#endif //RRANSAC_TRACK_INITIALIZATION_SEED_POLOCIES_NULL_POLICY_H_