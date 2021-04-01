#ifndef RRANSAC_TRACK_INITIALIZATION_RANSAC_H_
#define RRANSAC_TRACK_INITIALIZATION_RANSAC_H_
#pragma once


#include <stdlib.h>
#include <time.h>
#include <algorithm> 
#include <numeric>
#include <random>
#include <mutex>
#include <thread>

#include "rransac/system.h"
#include "rransac/data_containers/data_tree/data_tree_cluster.h"
#include "rransac/common/models/model_base.h"
#include "rransac/data_containers/cluster.h"
#include "rransac/common/models/model_manager.h"


namespace rransac {


/**
 * \class Ransac
 * A target that is not initialized (not being tracked by R-RANSAC) is
 * producing measurements captured by the measurement sources. Provided
 * that the measurements of the un-initialized target are not in the
 * validation region of any initialized target, its measurements are
 * not associated with any existing initialized target; thus, they are
 * associated with a cluster. After enough sensor scans,
 * there will be enough measurements in the cluster to observe the target
 * and initialize a new track. 
 * 
 * This process is done using RANSAC. RANSAC is a regression algorithm
 * that estimates the parameters of a model while mitigating the effect
 * of gross outliers. RANSAC is used to initialize new targets. 
 * 
 * RANSAC uses the measurements in a cluster to try to initialize a
 * new target by first generating many track hypotheses. A track hypothesis
 * is a current hypothetical state \f$x_{k}^{h}\f$, the subscript \f$k\f$ denotes
 * the current time, estimated from a random subset of measurements from
 * the cluster. The number of measurements in the subset is determined by Parameters::RANSAC_minimum_subset_.
 * 
 * Each track hypothesis is scored by finding the number of measurement inliers to 
 * the track hypothesis. A measurement is an inlier if the normalized geodesic distance between
 * the measurement and estimated measurement produced from the track hypothesis is within 
 * the distance perscribed by SourceParameters::RANSAC_inlier_probability_. This distance is normalized
 * by the innovation covariance. 
 * 
 * RANSAC generates up to Parameters::RANSAC_max_iters_ track hypothesis per cluster. If a track hypothesis
 * has a score greater than or equal to Paramters::RANSAC_score_stopping_criteria_, then RANSAC will stop generating
 * other track hypothesis and use most recent track hypothesis to generate a new track. If none of the track hypothesis
 * generated from a cluster has a score that meets the stopping criteria, then RANSAC takes the track hypothesis with the 
 * best score, provided that it's score is at least Parameters::RANSAC_score_minimum_requirement_, and generates a new track.
 * 
 * A new track is generated by propagating the track hypothesis back in time to the time stamp of the oldest measurement inlier. 
 * Measurements are then used to update the track according to the data association policy used by DataAssociationHost. The data 
 * association policy used is one of the template parameters of Ransac. The update 
 * consists of updating the state estimate, the error covariance, the model likelihood, and the consensus set. After the track
 * has been initialized, it is added to System via the member method ModelManager::AddModel. 
 * 
 * RANSAC can work with linear and nonlinear models. Depending on the type of model, the log maximum likelihood estimation (LMLE) 
 * policy needs to be selected. Currently, LinearLMLEPolicy should be used for linear systems, and NonLinearLMLEPolicy should be
 * used for nonlinear systems. Nonlinear estimation solvers can be seeded with initial conditions to speed up the solver and influence
 * what value it converges to. If you don't want to use any seeding methods, use NULLSeedPolicy. These policies are given to Ransac as 
 * template parameters. 
 */

template< typename tModel, template <typename > typename tSeed, template<typename , template <typename > typename > typename tLMLEPolicy, template<typename > typename tAssociationPolicy>
class Ransac : public tLMLEPolicy<tModel,tSeed> , tAssociationPolicy<tModel>{

public: 

typedef typename tModel::State State;                                       /**< The state of the target. @see State. */
typedef typename tModel::DataType DataType;                                 /**< The scalar object for the data. Ex. float, double, etc. */
typedef typename tModel::Source Source;                                     /**< The object type of the source. @see SourceBase. */
typedef tModel Model;                                                       /**< The object type of the model. @see ModelBase. */

Ransac();

~Ransac() = default;

/**
 * This is the only function that needs to be called. It performs RANSAC on every cluster in order
 * to initialize new tracks.
 * @param[in,out] system Contains all of the data used by R-RANSAC
 */ 
static void Run(System<tModel>& sys);

/**
 * Randomly selects a measurement from the latest time step and randomly selects 
 * measurements from other time steps until the number of measurements randomly selected is as perscribed by the parameter num_meas.
 * @param[in] num_meas The total number of random numbers to select. 
 * @param[in] cluster The cluster from which to the measurements will be sampled. 
 * @return Returns the indicess of the randomply sampled measurements. They are not gauranteed to be in chronological order
 */ 
static std::vector<typename Cluster<DataType>::IteratorPair> GenerateMinimumSubset(const unsigned int num_meas, Cluster<DataType>& cluster);

/**
 * Generates a hypothetical state at the current time step using the provided measurements in meas_subset.The method is determined by the log maximum likelihood estimation policy.
 * @param[in] meas_subset The container of iterators to measurements that will be used to estimate the hypothetical state.
 * @param[in] curr_time The current time.
 * @param[in] sources The vector of sources used. 
 */ 
static State GenerateHypotheticalStateEstimate(const std::vector<typename Cluster<DataType>::IteratorPair>& meas_subset, const System<tModel>& sys, bool& success){
    return Ransac::GenerateHypotheticalStateEstimatePolicy(meas_subset, sys,success);
}


/**
 * Finds all of the inliers to the hypothetical state estimate from the cluster and scores the estimate using the inliers.
 * The score is calculated as follows. Starting at the earliest time stamp of an inlier measurement, count the number of measurements
 * sources that produced inlier measurements at that time and add that number to the score. Then repeat this process for each following time stamp. 
 * This method is used since we assume that a measurement source can produce at most one true measurement per target every sensor scan. If a source 
 * produces two inliers at one time stamp, we count it as only one instead of two since both of them cannot be true measurements by our assumption. This 
 * method also weighs each measurement source equally in this way. 
 * @param[in] xh The hypothetical state estimate.
 * @param[in] cluster The cluster from which xh was generated.
 * @param[in] sys Contains the system information.
 * @param[in] inliers A vector of iterators from the measurements in cluster that are inliers to the hypothetical state estimate. The inliers are in chronological order.
 */ 
static int ScoreHypotheticalStateEstimate(const State& xh, Cluster<DataType>& cluster, const System<tModel>& sys, std::vector<typename Cluster<DataType>::IteratorPair>& inliers);


/**
 * Generates a track using the hypothetical state estimate with the best score and the inliers to the hypothetical state estimate.
 * @param[in] xh The best hypothetical state estimate generated from the cluster
 * @param[in] sys Contains the system information.
 * @param[in] inliers The measurements that support the hypothetical state estimate
 * 
 */ 
static tModel GenerateTrack(const State&xh, const System<tModel>& sys, const std::vector<typename Cluster<DataType>::IteratorPair>& inliers);

private:

/**
 * When filtering a new track's state estimate, the measurement weights and model likelihood need to be calculated
 * in order to properly update the track. We use the same policy in the data association method to get the measurement 
 * weights, validation volume, and model likelihood. 
 * 
 */ 
static void CalculateMeasurmentAndLikelihoodData(const System<tModel>& sys, tModel& model) {
    Ransac::CalculateMeasurmentAndLikelihoodDataPolicy(sys, model);
}

/**
 * Runs RANSAC on a single cluster. If a track is created, it will add the track to system 
 * and remove the inliers from the cluster. It is thread safe.
 * @param cluster_iter An iterator to a cluster on which RANSAC will be performed
 * @param sys An object containing all the R-RANSAC data
 */ 
static  void RunSingle(const typename std::list<Cluster<DataType>>::iterator& cluster_iter, System<tModel>& sys, std::mutex& mtx);


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
   
    if (num_meas > cluster.data_.size() || num_meas < 0) {
        std::cerr << "num_meas: " << num_meas << std::endl;
        std::cerr << "cluster size: " << cluster.data_.size() << std::endl;
        throw std::runtime_error("RANSAC: GenerateMinimumSubset: num_meas must be less than the number of different time steps of measurements in cluster and greater than zero.");
    }

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

        dt = outer_iter->begin()->time_stamp - current_time;                              // Get time difference
        xh_p = tModel::PropagateState(xh,dt);
        G = tModel::GetLinTransFuncMatNoise(xh_p,dt);
        Q_bar = G*sys.params_.process_noise_covariance_*G.transpose();
        std::fill(innov_cov_set.begin(), innov_cov_set.end(), false);                         // Reset vector since we are moving to a new time step
        std::fill(src_contributed.begin(), src_contributed.end(), false);                     // Reset vector since we are moving to a new time step

        for(auto inner_iter = outer_iter->begin(); inner_iter != outer_iter->end(); ++inner_iter) {


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



    // Create new track with state estimate at the same time step as the oldest inlier measurement
    tModel new_track;
    new_track.Init(sys.params_);
    double dt = inliers.begin()->inner_it->time_stamp - sys.current_time_;
    new_track.state_ = tModel::PropagateState(xh,dt); 
    double curr_time = 0;
    double propagate_time = inliers.begin()->inner_it->time_stamp;

    // Initialize objects
    std::vector<ModelLikelihoodUpdateInfo> model_likelihood_update_info(sys.sources_.size());                      


    for (auto iter = inliers.begin(); iter != inliers.end(); ) {
        curr_time = iter->inner_it->time_stamp;
        dt = curr_time - propagate_time;
        propagate_time = curr_time;
        new_track.PropagateModel(dt);

#ifdef DEBUG_BUILD
        if (isnan(new_track.state_.g_.data_(0.0))) {
            throw std::runtime_error("RANSAC::GenerateTrack New track state estimate has nan value. ");
        }
#endif


        while(iter != inliers.end() && iter->inner_it->time_stamp == curr_time) {

            new_track.AddNewMeasurement(*iter->inner_it);

            iter++;
        }
        CalculateMeasurmentAndLikelihoodData(sys,new_track);
        new_track.UpdateModel(sys.sources_,sys.params_);
        new_track.UpdateModelLikelihood(sys.sources_);
    }

    return new_track;

}

//-------------------------------------------------------------------------------------------------------------------------------------------


template< typename tModel, template <typename > typename tSeed, template<typename , template <typename > typename > typename tLMLEPolicy, template<typename > typename tAssociationPolicy>
void Ransac<tModel, tSeed, tLMLEPolicy, tAssociationPolicy>::RunSingle(const typename std::list<Cluster<DataType>>::iterator& cluster_iter, System<tModel>& sys, std::mutex& mtx) {

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
        if (success) {
            score = ScoreHypotheticalStateEstimate(hypothetical_state, *cluster_iter, sys, inliers);
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
        
        
        for (auto iter = best_inliers.begin(); iter != best_inliers.end(); ++iter) {
            info.iter_pair = *iter;
            mtx.lock();
            sys.data_tree_.RemoveMeasurement(info);
            mtx.unlock();
        }
        {
        mtx.lock();
        ModelManager<Model>::AddModel(sys,new_track);
        mtx.unlock();
        }
    }


    
}


//-------------------------------------------------------------------------------------------------------------------------------------------

template< typename tModel, template <typename > typename tSeed, template<typename , template <typename > typename > typename tLMLEPolicy, template<typename > typename tAssociationPolicy>
void Ransac<tModel, tSeed, tLMLEPolicy, tAssociationPolicy>::Run(System<tModel>& sys) {


std::vector<std::thread> threads;
std::mutex mtx;

for (auto& cluster_iter : sys.clusters_) {
    
    if (cluster_iter->data_.size() > sys.params_.RANSAC_minimum_subset_ && cluster_iter->Size() > sys.params_.RANSAC_score_minimum_requirement_ && cluster_iter->data_.back().front().time_stamp == sys.current_time_) {

        threads.push_back(std::thread(RunSingle,cluster_iter,std::ref(sys),std::ref(mtx)));
    }
}

for (int ii = 0; ii < threads.size(); ++ii) {
    threads.at(ii).join();
}

// Delete the pointers and wait for new ones to be given.
sys.clusters_.clear();

}

} //namespace rransac

#endif // RRANSAC_TRACK_INITIALIZATION_RANSAC_H_