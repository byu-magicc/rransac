#ifndef RRANSAC_TRACK_INITIALIZATION_SEED_POLOCIES_NULL_POLICY_H_
#define RRANSAC_TRACK_INITIALIZATION_SEED_POLOCIES_NULL_POLICY_H_
#pragma once

#include "rransac/system.h"
#include "rransac/data_containers/cluster.h"

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
typedef typename tModel::Base::TransformDataType TransformDataType;
typedef Cluster<DataType,TransformDataType> ClusterT;



/**
 * This function sets the initial guess of the state estimate to zero.
 * @param meas_subset The subset used to generate the seed for the LMLE algorithm
 * @param sys The object that contains the R-RANASAC data
 * @param x The initial state for the LMLE algorithm
 * @param size The size of the pointer x
 * 
 */ 
    static void GenerateSeedPolicy(const std::vector<typename ClusterT::IteratorPair>& meas_subset, const System<tModel>& sys, double x[tModel::cov_dim_], const int size) {
        std::default_random_engine generator;
        std::uniform_real_distribution<DataType> distribution(0.0,1.0);
        for (int ii = 0; ii < size; ++ii) {
            x[ii]=distribution(generator);
        }
    }

};



} // namespace rransac


#endif //RRANSAC_TRACK_INITIALIZATION_SEED_POLOCIES_NULL_POLICY_H_