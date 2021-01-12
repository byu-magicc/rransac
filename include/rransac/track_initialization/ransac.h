#ifndef RRANSAC_TRACK_INITIALIZATION_RANSAC_H_
#define RRANSAC_TRACK_INITIALIZATION_RANSAC_H_

#include "system.h"
#include "data_containers/cluster.h"
#include <stdlib.h>
#include <time.h>
#include <algorithm> 
#include <numeric>
#include <random>

namespace rransac {

template< typename tModel, template<typename ttModel> typename tLMLEPolicy>
class Ransac : public tLMLEPolicy<tModel> {

public: 

typedef typename tModel::State State;
typedef typename tModel::Source Source;

Ransac();

~Ransac() = default;

/**
 * Randomly selects a measurement from the latest time step and randomly selects 
 * measurements from other time steps until there are num_meas randomly selected.
 * @param num_meas The total number of random numbers to select. 
 * @param cluster The cluster from which to the measurements will be sampled. 
 * @return Returns the indicess of the randomply sampled measurements 
 */ 
static std::vector<Cluster::ConstIteratorPair> GenerateMinimumSubset(const unsigned int num_meas, const Cluster& cluster);

/**
 * Generates a hypothetical state at the current time step using the provided measurements in meas_subset.The method is determined by the policy
 * @param meas_subset The container of iterators to measurements that will be used to estimate the hypothetical state
 * @param curr_time The current time
 * @param sources The vector of sources used. 
 */ 
static State GenerateHypotheticalStateEstimate(const std::vector<Cluster::ConstIteratorPair>& meas_subset, const System<tModel>& sys){
    return GenerateHypotheticalStateEstimatePolicy(meas_subset, sys);
}


/**
 * Finds all of the inliers to the hypothetical state estimate xh from cluster and scores the estimate using the inliers.
 * @param xh The hypothetical state estimate
 * @param cluster The cluster from which xh was generated
 * @param sys Contains the system information
 * @param inliers A vector of iterators from the measurements in cluster that are inliers to the hypothetical state estimate. The inliers are in chronological order
 */ 
static int ScoreHypotheticalStateEstimate(const State& xh, const Cluster& cluster, const System<tModel>& sys, std::vector<Cluster::ConstIteratorPair>& inliers);

private:


};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                            Definitions
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template< typename tModel, template<typename ttModel> typename tLMLEPolicy>
Ransac<tModel, tLMLEPolicy>::Ransac() {
    srand(time(NULL));
}

//----------------------------------------------------------------------------------------------------------

template< typename tModel, template<typename ttModel> typename tLMLEPolicy>
std::vector<Cluster::ConstIteratorPair > Ransac<tModel, tLMLEPolicy>::GenerateMinimumSubset(const unsigned int num_meas, const  Cluster& cluster) {
   


    std::vector<Cluster::ConstIteratorPair> meas_index(num_meas);

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
    std::shuffle(time_indices.begin(), time_indices.end(), std::default_random_engine(time(NULL)) );

    // Randomly samples other measurements
    for (unsigned int ii = 0; ii < num_meas-1; ++ii) {

        meas_index[ii].outer_it = std::next(cluster.data_.begin(), time_indices[ii]);
        meas_index[ii].inner_it = std::next(meas_index[ii].outer_it->begin(), rand() % meas_index[ii].outer_it->size());


    }

    return meas_index;
    
}


//----------------------------------------------------------------------------------------------------------
template< typename tModel, template<typename ttModel> typename tLMLEPolicy>
int Ransac<tModel, tLMLEPolicy>::ScoreHypotheticalStateEstimate(const State& xh, const Cluster& cluster, const System<tModel>& sys, std::vector<Cluster::ConstIteratorPair>& inliers) {

    inliers.clear(); // Make sure it is empty
    int score = 0;
    int current_time = sys.current_time_;
    int dt = 0;
    int src_index = 0;
    double d = 0;          // The distance

    std::vector<Eigen::MatrixXd> innov_cov_inv(sys.sources_.size());
    std::vector<Meas> estimated_meas(sys.sources_.size());
    std::vector<bool> innov_cov_set(sys.sources_.size(),false);

    // Matrices to be used
    Eigen::MatrixXd G;
    Eigen::MatrixXd Q_bar;
    Eigen::MatrixXd V;
    Eigen::MatrixXd H;
    Eigen::MatrixXd e;

    // Find all of the inliers. The outer iterator passes through different time steps and the inner iterator passes through different measurements
    for (auto outer_iter = cluster.data_.begin(); outer_iter != cluster.data_.end(); ++outer_iter) {

        dt = outer_iter->begin()->time_stamp - current_time;                              // Get time difference
        G = tModel::GetLinTransFuncMatNoise(xh,dt);
        Q_bar = G*sys.params_.process_noise_covariance_*G.transpose();
        std::fill(innov_cov_set.begin(), innov_cov_set.end(), false);                     // Reset vector since we are moving to a new time step

        for(auto inner_iter = outer_iter->begin(); inner_iter != outer_iter.end(); ++inner_iter) {

            src_index = inner_iter->source_index;

            // Get the estimate measurement and the innovation covariance
            if(!innov_cov_set[src_index]) {
                H = tModel::GetLinObsMatState(sys.sources_,xh,src_index);
                V = tModel::GetLinObsMatSensorNoise(sys.sources_,xh,src_index);
                innov_cov_inv[src_index] = (V*sys.sources_[src_index].params_.meas_cov_*V.transpose() + Q_bar).inverse();
                estimated_meas[src_index] = sys.sources_[src_index].GetEstMeas(xh);
                innov_cov_set[src_index] = true;
            }

            e = sys.sources_[src_index].OMinus(*inner_iter, estimated_meas[src_index]);
            d = e.transpose()*innov_cov_inv[src_index]*e;


            // If the measurement is an inlier, add it. 
            if (d < sys.sources_[src_index].RANSAC_inlier_threshold_) {

            }

            
        }
    }


}

} //namespace rransac

#endif // RRANSAC_TRACK_INITIALIZATION_RANSAC_H_