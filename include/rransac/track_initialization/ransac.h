#ifndef RRANSAC_TRACK_INITIALIZATION_RANSAC_H_
#define RRANSAC_TRACK_INITIALIZATION_RANSAC_H_
#pragma once

#include "system.h"
#include "data_containers/data_tree/data_tree_cluster.h"
#include "common/models/model_base.h"
#include "data_containers/cluster.h"
#include <stdlib.h>
#include <time.h>
#include <algorithm> 
#include <numeric>
#include <random>
#include <mutex>
#include <thread>

namespace rransac {

template< typename tModel, template <typename > typename tSeed, template<typename , template <typename > typename > typename tLMLEPolicy, template<typename > typename tAssociationPolicy>
class Ransac : public tLMLEPolicy<tModel,tSeed> , tAssociationPolicy<tModel>{

public: 

typedef typename tModel::State State;
typedef typename tModel::Source Source;
typedef typename tModel::DataType DataType;

Ransac();

~Ransac() = default;

/**
 * This is the only function that needs to be called. It performs RANSAC on every cluster with sufficient number of measurements from 
 * different time steps. If a strong enough hypothetical state estimate is generated, then it will create a new track by filtering the 
 * hypothetical state estimate using its inliers.
 * @param system Contains all of the data used by R-RANSAC
 */ 
static void Run(System<tModel>& sys);

/**
 * Randomly selects a measurement from the latest time step and randomly selects 
 * measurements from other time steps until there are num_meas randomly selected.
 * @param num_meas The total number of random numbers to select. 
 * @param cluster The cluster from which to the measurements will be sampled. 
 * @return Returns the indicess of the randomply sampled measurements. They are not gauranteed to be in chronological order
 */ 
static std::vector<typename Cluster<DataType>::IteratorPair> GenerateMinimumSubset(const unsigned int num_meas, Cluster<DataType>& cluster);

/**
 * Generates a hypothetical state at the current time step using the provided measurements in meas_subset.The method is determined by the policy
 * @param meas_subset The container of iterators to measurements that will be used to estimate the hypothetical state
 * @param curr_time The current time
 * @param sources The vector of sources used. 
 */ 
static State GenerateHypotheticalStateEstimate(const std::vector<typename Cluster<DataType>::IteratorPair>& meas_subset, const System<tModel>& sys, bool& success){
    return Ransac::GenerateHypotheticalStateEstimatePolicy(meas_subset, sys,success);
}


/**
 * Finds all of the inliers to the hypothetical state estimate xh from cluster and scores the estimate using the inliers.
 * @param xh The hypothetical state estimate
 * @param cluster The cluster from which xh was generated
 * @param sys Contains the system information
 * @param inliers A vector of iterators from the measurements in cluster that are inliers to the hypothetical state estimate. The inliers are in chronological order
 */ 
static int ScoreHypotheticalStateEstimate(const State& xh, Cluster<DataType>& cluster, const System<tModel>& sys, std::vector<typename Cluster<DataType>::IteratorPair>& inliers);


/**
 * Generates a track using the hypothetical state estimate with the best score and the inliers to the hypothetical state estimate.
 * @param xh The best hypothetical state estimate generated from the cluster
 * @param inliers The measurements that support the hypothetical state estimate
 * 
 */ 
static tModel GenerateTrack(const State&xh, const System<tModel>& sys, const std::vector<typename Cluster<DataType>::IteratorPair>& inliers);

private:

static void CalculateMeasurmentAndLikelihoodData(const System<tModel>& sys, tModel& model) {
    Ransac::CalculateMeasurmentAndLikelihoodDataPolicy(sys, model);
}

/**
 * Runs RANSAC on a single cluster. If a track is created, it will add the track to system 
 * and remove the inliers from the cluster. It is thread safe.
 * @param cluster_iter An iterator to a cluster on which RANSAC will be performed
 * @param sys An object containing all the R-RANSAC data
 */ 
static  void RunSingle(const typename std::list<Cluster<DataType>>::iterator& cluster_iter, System<tModel>& sys);


};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                            Definitions
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template< typename tModel, template <typename > typename tSeed, template<typename , template <typename > typename > typename tLMLEPolicy, template<typename > typename tAssociationPolicy>
Ransac<tModel, tSeed, tLMLEPolicy, tAssociationPolicy>::Ransac() {
    srand(time(NULL));
}

//----------------------------------------------------------------------------------------------------------

template< typename tModel, template <typename > typename tSeed, template<typename , template <typename > typename > typename tLMLEPolicy, template<typename > typename tAssociationPolicy>
std::vector<typename Cluster<typename tModel::DataType>::IteratorPair > Ransac<tModel, tSeed, tLMLEPolicy, tAssociationPolicy>::GenerateMinimumSubset(const unsigned int num_meas,  Cluster<DataType>& cluster) {
   
    if (num_meas > cluster.data_.size() || num_meas < 0)
        throw std::runtime_error("RANSAC: GenerateMinimumSubset: num_meas must be less than the number of different time steps of measurements in cluster and greater than zero.");

    std::vector<typename Cluster<DataType>::IteratorPair> meas_index(num_meas);

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
template< typename tModel, template <typename > typename tSeed, template<typename , template <typename > typename > typename tLMLEPolicy, template<typename > typename tAssociationPolicy>
int Ransac<tModel, tSeed, tLMLEPolicy, tAssociationPolicy>::ScoreHypotheticalStateEstimate(const State& xh, Cluster<DataType>& cluster, const System<tModel>& sys, std::vector<typename Cluster<DataType>::IteratorPair>& inliers) {

    inliers.clear(); // Make sure it is empty
    int score = 0;
    double current_time = sys.current_time_;
    double dt = 0;
    int src_index = 0;
    double d = 0;          // The distance

    typename Cluster<DataType>::IteratorPair pair;
    typename tModel::State xh_p;         // propagated state

    std::vector<Eigen::MatrixXd> innov_cov_inv(sys.sources_.size());
    std::vector<Meas<DataType>> estimated_meas(sys.sources_.size());
    std::vector<bool> innov_cov_set(sys.sources_.size(),false);
    std::vector<bool> src_contributed(sys.sources_.size(),false);

    // Matrices to be used
    Eigen::MatrixXd G;
    Eigen::MatrixXd Q_bar;
    Eigen::MatrixXd V;
    Eigen::MatrixXd H;
    Eigen::MatrixXd e;

    // Find all of the inliers. The outer iterator passes through different time steps and the inner iterator passes through different measurements
    for (auto outer_iter = cluster.data_.begin(); outer_iter != cluster.data_.end(); ++outer_iter) {

        // std::cerr << "outer " << std::endl;

        dt = outer_iter->begin()->time_stamp - current_time;                              // Get time difference
        xh_p = tModel::PropagateState(xh,dt);
        G = tModel::GetLinTransFuncMatNoise(xh_p,dt);
        Q_bar = G*sys.params_.process_noise_covariance_*G.transpose();
        std::fill(innov_cov_set.begin(), innov_cov_set.end(), false);                         // Reset vector since we are moving to a new time step
        std::fill(src_contributed.begin(), src_contributed.end(), false);                     // Reset vector since we are moving to a new time step

        for(auto inner_iter = outer_iter->begin(); inner_iter != outer_iter->end(); ++inner_iter) {

            // std::cerr << "inner " << std::endl;


            src_index = inner_iter->source_index;

            // Get the estimate measurement and the innovation covariance
            if(!innov_cov_set[src_index]) {
                H = tModel::GetLinObsMatState(sys.sources_,xh_p,src_index);
                V = tModel::GetLinObsMatSensorNoise(sys.sources_,xh_p,src_index);
                innov_cov_inv[src_index] = (V*sys.sources_[src_index].params_.meas_cov_*V.transpose() + H*Q_bar*H.transpose()).inverse();
                estimated_meas[src_index] = sys.sources_[src_index].GetEstMeas(xh_p);
                estimated_meas[src_index].type = sys.sources_[src_index].params_.type_;
                innov_cov_set[src_index] = true;
            }

            e = sys.sources_[src_index].OMinus(*inner_iter, estimated_meas[src_index]);
            d = (e.transpose()*innov_cov_inv[src_index]*e)(0,0);


            // If the measurement is an inlier, add it. 
            if (d < sys.sources_[src_index].params_.RANSAC_inlier_threshold_) {
                pair.outer_it = outer_iter;
                pair.inner_it = inner_iter;
                inliers.push_back(pair);

                // If the measurement is the first inlier from the source src_index at this time, increment the score once. 
                // This is becuase multiple measurements from the same measurements source observed at the same time can only count as one total.
                if (!src_contributed[src_index]) {
                    score++;
                    src_contributed[src_index] = true;
                }
            }            
        }
    }

    return score;
}

//-------------------------------------------------------------------------------------------------------------------------------------------
template< typename tModel, template <typename > typename tSeed, template<typename , template <typename > typename > typename tLMLEPolicy, template<typename > typename tAssociationPolicy>
tModel Ransac<tModel, tSeed, tLMLEPolicy, tAssociationPolicy>::GenerateTrack(const State&xh, const System<tModel>& sys, const std::vector<typename Cluster<DataType>::IteratorPair>& inliers) {


    // std::cerr << "State g: " << std::endl << xh.g_.data_ << std::endl;
    // std::cerr << "State u: " << std::endl << xh.u_.data_ << std::endl << std::endl;

    // Create new track with state estimate at the same time step as the oldest inlier measurement
    tModel new_track;
    new_track.Init(sys.params_);
    double dt = inliers.begin()->inner_it->time_stamp - sys.current_time_;
    new_track.state_ = tModel::PropagateState(xh,dt); 
    double curr_time = 0;
    double propagate_time = inliers.begin()->inner_it->time_stamp;

    // std::cerr << "before new_track g: " << std::endl << new_track.state_.g_.data_ << std::endl;
    // std::cerr << "before new_track u: " << std::endl << new_track.state_.u_.data_ << std::endl << std::endl;

    // Initialize objects
    std::vector<ModelLikelihoodUpdateInfo> model_likelihood_update_info(sys.sources_.size());                      


    for (auto iter = inliers.begin(); iter != inliers.end(); ) {
        curr_time = iter->inner_it->time_stamp;
        dt = curr_time - propagate_time;
        propagate_time = curr_time;
        new_track.PropagateModel(dt);


        // if (isnan(new_track.state_.g_.data_(0.0))) {

        //     std::cerr << "current time: " << curr_time << std::endl;
        //     std::cerr << "new_track g: " << std::endl << new_track.state_.g_.data_ << std::endl;
        //     std::cerr << "new_track u: " << std::endl << new_track.state_.u_.data_ << std::endl << std::endl;
        // }

        while(iter != inliers.end() && iter->inner_it->time_stamp == curr_time) {

            new_track.AddNewMeasurement(*iter->inner_it);

            iter++;
        }
        CalculateMeasurmentAndLikelihoodData(sys,new_track);
        new_track.UpdateModel(sys.sources_,sys.params_);
        new_track.UpdateModelLikelihood(sys.sources_);
    }


    // std::cerr << "new_track g: " << std::endl << new_track.state_.g_.data_ << std::endl;
    // std::cerr << "new_track u: " << std::endl << new_track.state_.u_.data_ << std::endl << std::endl;

    return new_track;

}

//-------------------------------------------------------------------------------------------------------------------------------------------


template< typename tModel, template <typename > typename tSeed, template<typename , template <typename > typename > typename tLMLEPolicy, template<typename > typename tAssociationPolicy>
void Ransac<tModel, tSeed, tLMLEPolicy, tAssociationPolicy>::RunSingle(const typename std::list<Cluster<DataType>>::iterator& cluster_iter, System<tModel>& sys) {

    int best_score = 0;
    int score = 0;
    bool success = false;
    unsigned int iterations = 0;
    unsigned int iteration_stopping_criteria = 10; // If after this many iterations the best score is still zero, then it will terminate. 
    std::vector<typename Cluster<DataType>::IteratorPair> meas_subset;
    typename tModel::State hypothetical_state;
    typename tModel::State best_hypothetical_state;
    std::vector<typename Cluster<DataType>::IteratorPair> inliers;
    std::vector<typename Cluster<DataType>::IteratorPair> best_inliers;
    while (best_score < sys.params_.RANSAC_score_stopping_criteria_ && iterations < sys.params_.RANSAC_max_iters_) {
        meas_subset = GenerateMinimumSubset(sys.params_.RANSAC_minimum_subset_, *cluster_iter);
        hypothetical_state = GenerateHypotheticalStateEstimate(meas_subset, sys,success);
        // std::cerr << "hypothetical_state g: " << std::endl << hypothetical_state.g_.data_ << std::endl;
        // std::cerr << "hypothetical_state u: " << std::endl << hypothetical_state.u_.data_ << std::endl << std::endl;
        // std::cerr << "success: " << success << std::endl;
        if (success) {
            score = ScoreHypotheticalStateEstimate(hypothetical_state, *cluster_iter, sys, inliers);
            // std::cerr << "score: " << score << std::endl << std::endl;
        } else
            score = -1;
        

        if (score > best_score) {
            best_hypothetical_state = hypothetical_state;
            best_inliers = inliers;
            best_score = score;
        }

        ++iterations;

        // The first few estimates are so bad that the measurements used to create them aren't even inliers.
        // There is a good chance that the cluster just has a clutter of noisy measurements so terminate early.
        if (iterations >= iteration_stopping_criteria && best_score < sys.params_.RANSAC_minimum_subset_) {
            break;
        }

    }

    if (best_score > sys.params_.RANSAC_score_minimum_requirement_) {
        tModel new_track =  GenerateTrack(best_hypothetical_state, sys, best_inliers);
        typename DataTreeClusters<DataType>::MeasurementLocationInfo info;
        info.cluster_iter = cluster_iter;
        std::mutex m;
        
        for (auto iter = best_inliers.begin(); iter != best_inliers.end(); ++iter) {
            info.iter_pair = *iter;
            {
            std::lock_guard<std::mutex> lockGuard(m);
            sys.data_tree_.RemoveMeasurement(info);
            }
        }
        {
        std::lock_guard<std::mutex> lockGuard(m);
        sys.models_.push_back(new_track);
        }
    }


    
}


//-------------------------------------------------------------------------------------------------------------------------------------------

template< typename tModel, template <typename > typename tSeed, template<typename , template <typename > typename > typename tLMLEPolicy, template<typename > typename tAssociationPolicy>
void Ransac<tModel, tSeed, tLMLEPolicy, tAssociationPolicy>::Run(System<tModel>& sys) {


std::vector<std::thread> threads;

for (auto& cluster_iter : sys.clusters_) {
    
    if (cluster_iter->Size() > sys.params_.RANSAC_score_minimum_requirement_ && cluster_iter->data_.back().front().time_stamp == sys.current_time_) {

        threads.push_back(std::thread(RunSingle,cluster_iter,std::ref(sys)));
    }
}

for (auto& th: threads) {
    th.join();
}

// Delete the pointers and wait for new ones to be given.
sys.clusters_.empty();

}

} //namespace rransac

#endif // RRANSAC_TRACK_INITIALIZATION_RANSAC_H_