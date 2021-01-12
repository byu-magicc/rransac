#ifndef RRANSAC_TRACK_INITIALIZATION_LINEAR_LMLE_POLICY_H_
#define RRANSAC_TRACK_INITIALIZATION_LINEAR_LMLE_POLICY_H_

#include "data_containers/cluster.h"
#include "system.h"
#include "common/measurement/measurement_base.h"
#include <Eigen/Dense>

namespace rransac
{

template<typename tModel>    
class LinearLMLEPolicy {

public:

typedef typename tModel::State State;
typedef typename tModel::Source Source;
typedef tModel Model;

/**
 * Generates a hypothetical state at the current time step using the provided measurements in meas_subset.
 * @param meas_subset The container of iterators to measurements that will be used to estimate the hypothetical state
 * @param curr_time The current time
 * @param sources The vector of sources used. 
 */ 
static State GenerateStateEstimatePolicy(const std::vector<Cluster::ConstIteratorPair>& meas_subset, const double curr_time, const std::vector<Source>& sources);

private:


Eigen::MatrixXd ConstructInnovationCovariance(const State, const double dt);

};



/////////////////////////////////////////////////////////////////////////////////////
//                Definitions
/////////////////////////////////////////////////////////////////////////////////////
template<typename tModel>    
typename tModel::State LinearLMLEPolicy<tModel>::GenerateStateEstimatePolicy(const std::vector<Cluster::ConstIteratorPair>& meas_subset, const double curr_time, const std::vector<typename tModel::Source>& sources){
    
    typename tModel::State x;   // hypothetical state

    Meas& meas_ref;

    for (auto iter = meas_subset.begin(); iter != meas_subset.end(); ++iter) {

    }

}

//-------------------------------------------------------------------------------------------
template<typename tModel>
Eigen::MatrixXd LinearLMLEPolicy<tModel>::ConstructInnovationCovariance(const typename tModel::State, const double dt) {

}

} // namespace rransac

#endif // RRANSAC_TRACK_INITIALIZATION_LINEAR_LMLE_POLICY_H_