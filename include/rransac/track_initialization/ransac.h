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
#include "rransac/common/data_association/data_association_host.h"


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

template< typename _Model, template <typename > typename _Seed, template<typename , template <typename > typename > typename _LMLEPolicy, template<class> typename _ValidationRegionPolicy, template<class> typename _UpdateTrackLikelihoodPolicy, template<class> typename _MeasurementWeightPolicy>
class Ransac : public _LMLEPolicy<_Model,_Seed> , _ValidationRegionPolicy<_Model>, _UpdateTrackLikelihoodPolicy<_Model>, _MeasurementWeightPolicy<_Model>{

public: 

typedef _Model Model;                                                       /**< The object type of the model. @see ModelBase. */
typedef typename _Model::State State;                                       /**< The state of the target. @see State. */
typedef typename _Model::DataType DataType;                                 /**< The scalar object for the data. Ex. float, double, etc. */
typedef typename _Model::SourceContainer SourceContainer;                   /**< The object type of the source. @see SourceBase. */
typedef System<Model> Sys;    
typedef typename Sys::ClusterT ClusterT;     
typedef typename Sys::Measurement Measurement; 
typedef std::vector<typename ClusterT::IteratorPair> VecClusterIterPair;  
typedef typename Sys::DataTreeClustersT DataTreeClustersT;   
typedef typename Model::Base::TransformDataType TransformDataType;
typedef DataAssociationInfo<TransformDataType> DataAssociationInfoT;                                    

Ransac();

~Ransac() = default;

/**
 * This is the only function that needs to be called. It performs RANSAC on every cluster in order
 * to initialize new tracks.
 * @param[in,out] system Contains all of the data used by R-RANSAC
 */ 
static void Run(Sys& sys);

/**
 * Randomly selects a measurement from the latest time step and randomly selects 
 * measurements from other time steps until the number of measurements randomly selected is as perscribed by the parameter num_meas.
 * @param[in] num_meas The total number of random numbers to select. 
 * @param[in] cluster The cluster from which to the measurements will be sampled. 
 * @return Returns the indicess of the randomply sampled measurements. They are not gauranteed to be in chronological order
 */ 
static VecClusterIterPair GenerateMinimumSubset(const unsigned int num_meas, ClusterT& cluster);

/**
 * Generates a hypothetical state at the current time step using the provided measurements in meas_subset.The method is determined by the log maximum likelihood estimation policy.
 * @param[in] meas_subset The container of iterators to measurements that will be used to estimate the hypothetical state.
 * @param[in] curr_time The current time.
 * @param[in] sources The vector of sources used. 
 */ 
static State GenerateHypotheticalStateEstimate(const VecClusterIterPair& meas_subset, const Sys& sys, bool& success){
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
static int ScoreHypotheticalStateEstimate(const State& xh, ClusterT& cluster, const Sys& sys, VecClusterIterPair& inliers);


/**
 * Generates a track using the hypothetical state estimate with the best score and the inliers to the hypothetical state estimate.
 * @param[in] xh The best hypothetical state estimate generated from the cluster
 * @param[in] sys Contains the system information.
 * @param[in] inliers The measurements that support the hypothetical state estimate
 * 
 */ 
static Model GenerateTrack(const State&xh, const Sys& sys, const VecClusterIterPair& inliers);

private:

/**
 *  Determines if the measurement falls inside the validation region of the track according to the policy. 
 * @param[in] sys The object that contains all of the data of RRANSAC. This includes the new measurements, tracks, clusters and data tree. 
 * @param[in] meas The measurement
 * @param[in] track The track. It is not a constant reference so be careful!! 
 */ 
static bool InValidationRegion(const Sys& sys, const Measurement& meas, Model& track) {
    return Ransac::PolicyInValidationRegion(sys, meas,track);
}

/** 
 * Uses the new associated measurements to update the track's likelihood. 
 * @param[in,out] sys The object that contains all of the data of RRANSAC. This includes the new measurements, tracks, clusters and data tree. 
 */
static void UpdateTrackLikelihoodSingle(const Sys& sys, Model& track, DataAssociationInfoT& info, const double dt  ) {
    Ransac::PolicyUpdateTrackLikelihoodSingle(sys,track, info, dt);
}

/** 
 * Calculates the weights for each track associated measurement. 
 * @param[in,out] sys The object that contains all of the data of RRANSAC. This includes the new measurements, tracks, clusters and data tree. 
 */
static void CalculateMeasurementWeightSingle(const Sys& sys, Model& track, DataAssociationInfoT& info) {
    Ransac::PolicyCalculateMeasurementWeightSingle(sys,track,info);
}

/**
 * Runs RANSAC on a single cluster. If a track is created, it will add the track to system 
 * and remove the inliers from the cluster. It is thread safe.
 * @param cluster_iter An iterator to a cluster on which RANSAC will be performed
 * @param sys An object containing all the R-RANSAC data
 */ 
static  void RunSingle(const typename std::list<ClusterT>::iterator& cluster_iter, Sys& sys, std::mutex& mtx);


};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                            Definitions
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template< typename _Model, template <typename > typename _Seed, template<typename , template <typename > typename > typename _LMLEPolicy, template<class> typename _ValidationRegionPolicy, template<class> typename _UpdateTrackLikelihoodPolicy, template<class> typename _MeasurementWeightPolicy>
Ransac<_Model, _Seed, _LMLEPolicy, _ValidationRegionPolicy, _UpdateTrackLikelihoodPolicy, _MeasurementWeightPolicy>::Ransac() {
    srand(time(NULL));
}

//----------------------------------------------------------------------------------------------------------

template< typename _Model, template <typename > typename _Seed, template<typename , template <typename > typename > typename _LMLEPolicy, template<class> typename _ValidationRegionPolicy, template<class> typename _UpdateTrackLikelihoodPolicy, template<class> typename _MeasurementWeightPolicy>
typename Ransac<_Model, _Seed, _LMLEPolicy, _ValidationRegionPolicy, _UpdateTrackLikelihoodPolicy, _MeasurementWeightPolicy>::VecClusterIterPair Ransac<_Model, _Seed, _LMLEPolicy, _ValidationRegionPolicy, _UpdateTrackLikelihoodPolicy, _MeasurementWeightPolicy>::GenerateMinimumSubset(const unsigned int num_meas,  ClusterT& cluster) {
   
    if (num_meas > cluster.data_.size() || num_meas < 0) {
        std::cerr << "num_meas: " << num_meas << std::endl;
        std::cerr << "cluster size: " << cluster.data_.size() << std::endl;
        throw std::runtime_error("RANSAC: GenerateMinimumSubset: num_meas must be less than the number of different time steps of measurements in cluster and greater than zero.");
    }

    VecClusterIterPair meas_index(num_meas);

    // Get a random measurement from current time step
    meas_index.back().outer_it = std::prev(cluster.data_.end());
    int length = cluster.data_.back().size();
    int index = rand() % length;
    meas_index.back().inner_it = std::next(cluster.data_.back().begin(),index);

    if (num_meas == 1)
        return meas_index;

    // Get the time indices for the other randomly sampled measurements
    std::vector<int> time_indices(cluster.data_.size()-2);
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
template< typename _Model, template <typename > typename _Seed, template<typename , template <typename > typename > typename _LMLEPolicy, template<class> typename _ValidationRegionPolicy, template<class> typename _UpdateTrackLikelihoodPolicy, template<class> typename _MeasurementWeightPolicy>
int Ransac<_Model, _Seed, _LMLEPolicy, _ValidationRegionPolicy, _UpdateTrackLikelihoodPolicy, _MeasurementWeightPolicy>::ScoreHypotheticalStateEstimate(const State& xh, ClusterT& cluster, const Sys& sys, VecClusterIterPair& inliers) {

    inliers.clear(); // Make sure it is empty
    int score = 0;
    double current_time = sys.current_time_;
    double dt = 0;
    int src_index = 0;
    double d = 0;          // The distance

    typename ClusterT::IteratorPair pair;
    Model track;
    track.Init(sys.params_);

    std::vector<bool> innov_cov_set(sys.source_container_.num_sources_,false);
    std::vector<bool> src_contributed(sys.source_container_.num_sources_,false);


    // Find all of the inliers. The outer iterator passes through different time steps and the inner iterator passes through different measurements
    for (auto outer_iter = cluster.data_.begin(); outer_iter != cluster.data_.end(); ++outer_iter) {

        dt = outer_iter->begin()->time_stamp - current_time;                              // Get time difference
        track.state_ = xh;
        track.PropagateModel(dt);
        std::fill(src_contributed.begin(), src_contributed.end(), false);                     // Reset vector since we are moving to a new time step

        for(auto inner_iter = outer_iter->begin(); inner_iter != outer_iter->end(); ++inner_iter) {

            // If the measurement is an inlier, add it. 
            if (InValidationRegion(sys,*inner_iter,track)) {
                pair.outer_it = outer_iter;
                pair.inner_it = inner_iter;
                inliers.push_back(pair);

                // If the measurement is the first inlier from the source src_index at this time, increment the score once. 
                // This is becuase multiple measurements from the same measurements source observed at the same time can only count as one total.
                if (!src_contributed[inner_iter->source_index]) {
                    score++;
                    src_contributed[inner_iter->source_index] = true;
                }
            }            
        }
    }

    return score;
}

//-------------------------------------------------------------------------------------------------------------------------------------------
template< typename _Model, template <typename > typename _Seed, template<typename , template <typename > typename > typename _LMLEPolicy, template<class> typename _ValidationRegionPolicy, template<class> typename _UpdateTrackLikelihoodPolicy, template<class> typename _MeasurementWeightPolicy>
_Model Ransac<_Model, _Seed, _LMLEPolicy, _ValidationRegionPolicy, _UpdateTrackLikelihoodPolicy, _MeasurementWeightPolicy>::GenerateTrack(const State&xh, const Sys& sys, const VecClusterIterPair& inliers) {


    // std::cout << "ransac state: " << std::endl << xh.g_.data_ << std::endl << xh.u_.data_ << std::endl;

    TransformDataType empty_transform;
    DataAssociationInfoT data_association_info;
    data_association_info.source_produced_measurements_.resize(sys.source_container_.num_sources_,false);
    data_association_info.transform_state_.resize(sys.source_container_.num_sources_,false);
    data_association_info.transform_data_t_m_.resize(sys.source_container_.num_sources_,empty_transform);
    // Create new track with state estimate at the same time step as the oldest inlier measurement
    Model new_track;
    new_track.Init(sys.params_);
    double dt = inliers.begin()->inner_it->time_stamp - sys.current_time_;
    new_track.state_ = xh;
    Model::PropagateState(new_track.state_,dt); 
    double curr_time = 0;
    double propagate_time = inliers.begin()->inner_it->time_stamp;

    // Initialize objects
    std::vector<ModelLikelihoodUpdateInfo> model_likelihood_update_info(sys.source_container_.num_sources_);                      


    for (auto iter = inliers.begin(); iter != inliers.end(); ) {
        curr_time = iter->inner_it->time_stamp;
        dt = curr_time - propagate_time;
        propagate_time = curr_time;
        new_track.PropagateModel(dt);

        // Reset the data associaiton for each timestep
        std::fill(data_association_info.source_produced_measurements_.begin(),data_association_info.source_produced_measurements_.end(),false);
        std::fill(data_association_info.transform_state_.begin(),data_association_info.transform_state_.end(),false);
        std::fill(data_association_info.transform_data_t_m_.begin(),data_association_info.transform_data_t_m_.end(),empty_transform);

#ifdef DEBUG_BUILD
        if (isnan(new_track.state_.g_.data_(0.0))) {
            throw std::runtime_error("RANSAC::GenerateTrack New track state estimate has nan value. ");
        }
#endif


        while(iter != inliers.end() && iter->inner_it->time_stamp == curr_time) {

            new_track.AddNewMeasurement(*iter->inner_it);

            if (data_association_info.source_produced_measurements_[iter->inner_it->source_index]==false) {
                data_association_info.source_produced_measurements_[iter->inner_it->source_index]= true;
                data_association_info.transform_state_[iter->inner_it->source_index] = iter->inner_it->transform_state;
                data_association_info.transform_data_t_m_[iter->inner_it->source_index] = iter->inner_it->transform_data_t_m;
            }

            iter++;
        }
        
        // std::cout << "det: " << new_track.err_cov_.determinant() << std::endl;
        // std::cout << "norm: " << new_track.err_cov_.norm() << std::endl;
        UpdateTrackLikelihoodSingle(sys,new_track,data_association_info,dt);
        CalculateMeasurementWeightSingle(sys,new_track,data_association_info);
        new_track.UpdateModel(sys.source_container_,sys.params_);
    }

    return new_track;

}

//-------------------------------------------------------------------------------------------------------------------------------------------


template< typename _Model, template <typename > typename _Seed, template<typename , template <typename > typename > typename _LMLEPolicy, template<class> typename _ValidationRegionPolicy, template<class> typename _UpdateTrackLikelihoodPolicy, template<class> typename _MeasurementWeightPolicy>
void Ransac<_Model, _Seed, _LMLEPolicy, _ValidationRegionPolicy, _UpdateTrackLikelihoodPolicy, _MeasurementWeightPolicy>::RunSingle(const typename std::list<ClusterT>::iterator& cluster_iter, Sys& sys, std::mutex& mtx) {

    int best_score = 0;
    int score = 0;
    bool success = false;
    unsigned int iterations = 0;
    unsigned int iteration_stopping_criteria = 10; // If after this many iterations the best score is still zero, then it will terminate. 
    VecClusterIterPair meas_subset;
    State hypothetical_state;
    State best_hypothetical_state;
    VecClusterIterPair inliers;
    VecClusterIterPair best_inliers;
    while (best_score < sys.params_.RANSAC_score_stopping_criteria_ && iterations < sys.params_.RANSAC_max_iters_) {
        meas_subset = GenerateMinimumSubset(sys.params_.RANSAC_minimum_subset_, *cluster_iter);
        hypothetical_state = GenerateHypotheticalStateEstimate(meas_subset, sys,success);

        // std::cout << "hypothetical state pose: " << std::endl << hypothetical_state.g_.data_ << std::endl;
        // std::cout << "hypothetical state twist: " << std::endl << hypothetical_state.u_.data_ << std::endl;

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
        Model new_track =  GenerateTrack(best_hypothetical_state, sys, best_inliers);
        typename DataTreeClustersT::MeasurementLocationInfo info;
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

template< typename _Model, template <typename > typename _Seed, template<typename , template <typename > typename > typename _LMLEPolicy, template<class> typename _ValidationRegionPolicy, template<class> typename _UpdateTrackLikelihoodPolicy, template<class> typename _MeasurementWeightPolicy>
void Ransac<_Model, _Seed, _LMLEPolicy, _ValidationRegionPolicy, _UpdateTrackLikelihoodPolicy, _MeasurementWeightPolicy>::Run(Sys& sys) {


std::vector<std::thread> threads;
std::mutex mtx;

for (auto& cluster_iter : sys.clusters_) {
    
    if (cluster_iter->data_.size() > sys.params_.RANSAC_minimum_subset_ && cluster_iter->Size() > sys.params_.RANSAC_score_minimum_requirement_ && cluster_iter->data_.back().front().time_stamp == sys.current_time_) {

        threads.push_back(std::thread(RunSingle,cluster_iter,std::ref(sys),std::ref(mtx)));
        // RunSingle(cluster_iter,sys,mtx);
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