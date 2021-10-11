#ifndef RRANSAC_COMMON_MODEL_MANAGER_H_
#define RRANSAC_COMMON_MODEL_MANAGER_H_
#pragma once


#include <random>
#include <unsupported/Eigen/MatrixFunctions>

#include "rransac/system.h"
#include "rransac/common/utilities.h"
#include "rransac/parameters.h"

namespace rransac
{


/**
 * \class ModelManager
 * This class is responsible for managing all of the tracks. It is designed to 
 * facilitate adding new tracks, propagating and updating all the tracks, transforming the tracks
 * and managing the tracks. Track management consits of pruning the track's consensus sets,
 * merging similar tracks, pruning tracks, and ranking the tracks. 
 * 
 */ 
template <typename _Model>
class ModelManager {



public:

typedef _Model Model;
typedef System<Model> Sys;
typedef typename Sys::State State;

/**
* Add a new model. If the number of models is greater than the max number of models, then
* the function PruneModels is called.
* @param[in,out] sys The object that contains all of the data of RRANSAC. Thus it contains all of the tracks. 
* @param[in] model The model to be added.
*/
static void AddModel(Sys& sys, const Model& model);

/**
* Propagates every track by calling ModelBase::PropagateModel method on each 
* track.
* @param[in,out] sys The object that contains all of the data of RRANSAC. Thus it contains all of the tracks. 
* @param[in] dt The amount of time the track needs to be propagated. 
*/
static void PropagateModels(Sys& sys, const double dt);

/**
* Updates every track by calling ModelBase::UpdateModel method on each track.
* @param[in,out] sys The object that contains all of the data of RRANSAC. Thus it contains all of the tracks. 
*/
static void UpdateModels(Sys& sys);


/**
* Updates the parameters of every track. 
* @param[in,out] sys The object that contains all of the data of RRANSAC. Thus it contains all of the tracks and the system parameters 
*/
static void SetModelParameters(Sys& sys){ 
    for (auto iter = sys.models_.begin(); iter!=sys.models_.end(); ++iter) {
        iter->SetParameters(sys.params_);
    }
}

/**
 * Transforms all of the tracks and their consensus sets from the previous tracking frame to the current tracking frame. 
 * The consensus sets are only transformed if the flag Parameters::transform_consensus_set_ is set to true.
 * @param[in,out] sys The object that contains all of the data of RRANSAC. Thus it contains all of the tracks. 
 */ 
static void TransformModels(Sys& sys){
    for (auto iter = sys.models_.begin(); iter!=sys.models_.end(); ++iter) {
        iter->TransformModel(sys.transformaion_);
        if (sys.params_.transform_consensus_set_)
            iter->TransformConsensusSet(sys.transformaion_);
    }
}


/**
 * This function manages the tracks by pruning their consensus sets, merging similar tracks, pruning tracks, and ranking tracks.
 * @param[in,out] sys The object that contains all of the data of RRANSAC. Thus it contains all of the tracks. 
 * @param[in] expiration_time The expiration time of the measurements in the consensus sets. All measurements with a time stamp before the expiration 
 * time are removed.
 */
static void ManageModels(Sys& sys, const double expiration_time) {
    PruneConsensusSets(sys, expiration_time);
    MergeModels(sys);
    PruneModels(sys);
    RankModels(sys);

}



private:

/**
* Prunes the consensus set for each track by calling the method ConsensusSet::PruneConsensusSet on each consensus set.
* @param[in,out] sys The object that contains all of the data of RRANSAC. Thus it contains all of the tracks. 
* @param[in] expiration_time The expiration time of the measurements in the consensus sets. All measurements with a time stamp before the expiration 
* time are removed.
*/
static void PruneConsensusSets(Sys& sys, const double expiration_time);


/**
* Looks for similar models and merges them. This method is similar to the track 2 track fusion method discussed in 
* Tracking and Data Fusion by Bar-Shalom 2011.
* @param[in,out] sys The object that contains all of the data of RRANSAC. Thus it contains all of the tracks. 
*/
static void MergeModels(Sys& sys);

/**
 * Ranks the models according to their model likelihood. If the model likelihood is above the threshold Parameters::track_good_model_threshold_, then
 * it is considered a good model; otherwise, a poor model. Good models are added to System::good_models_. If there are more tracks than Parameters::track_max_num_tracks_, the tracks with
 * the lowest model likelihood are removed until there are only Parameters::track_max_num_tracks_ number of tracks. 
 * @param[in,out] sys The object that contains all of the data of RRANSAC. Thus it contains all of the tracks. 
 */
static void RankModels(Sys& sys); 

/**
* If there are more tracks than Parameters::track_max_num_tracks_, the tracks with
* the lowest model likelihood are removed until there are only Parameters::track_max_num_tracks_ number of tracks.
* @param[in,out] sys The object that contains all of the data of RRANSAC. Thus it contains all of the tracks.  
*/
static void PruneModels(Sys& sys);

/**
* Tests to see if two tracks are similar by weighing the geodesic distance between the track's states by their error covariances.
* @param[in] sys The object that contains all of the data of RRANSAC. Thus it contains all of the tracks. 
* @param[in] model1 One of the two tracks that will be compared.
* @param[in] model2 One of the two tracks that will be compared.
* @return Returns true if the tracks are similar
*/
static bool SimilarModels(const Sys& sys, const Model& model1, const Model& model2);

/**
 * Fuse two tracks together using the sampled covariance intersection method. The method is described in "A no-loss covariance intersection algorithm
 * for track-to-track fusion" by Xin Tian 2010. The method also merges the the two consensus sets together. If one or both of the tracks are good tracks, the lowest 
 * label is kept because it corresponds to the first track that was promoted to a good track. We do not properly fuse the model likelihoods together. To properly fuse them 
 * would require keeping a history of their ModelBase::model_likelihood_update_info_, which could require a lot of data storage depending on how long the tracks have existed.
 * So we simplified the prosecces by assigning the largetest model likelihood to the fused track.
 * @param[in,out] model1 One of the two tracks that will be fused together. model1 will become the merged track.
 * @param[in] model2 One of the two tracks that will be fused together.
 */ 
static void FuseModels(Model& model1, const Model& model2);



};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                            Definitions
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename _Model>
void ModelManager<_Model>::AddModel(Sys& sys, const Model& model) {

    sys.models_.push_back(model);
    sys.accumulative_number_of_tracks_++;

}

//-------------------------------------------------------------------------------------------------------------------

template <typename _Model>
void ModelManager<_Model>::PruneModels(Sys& sys) {

    //
    // IDEA:: WE MIGHT WANT TO PRUNE MODELS WHOSE LIKELIHOOD FALL BELOW A THRESHOLD
    // THE REASON BEHIND THIS IS THAT A MODEL MIGHT NOT GET A MEASUREMENT FOR A WHILE
    // BECUASE IT IS HIDDEN OR SOMETHING, BUT BY THE TIME IT DOES RECEIVE A NEW MEASUREMENT
    // ITS LIKELIHOOD COULD BE SO LOW THAT IT TAKES A WHILE TO HAVE A POSITIVE LIKELIHOOD AGAIN

    // Remove the models that have not received a measurement for a while
    for (auto it_expired_model = sys.models_.begin(); it_expired_model !=sys.models_.end();) {
        // std::cerr << "current time: " << sys.current_time_ << std::endl;
        // std::cerr << "it_expired_model->newest_measurement_time_stamp): " << it_expired_model->newest_measurement_time_stamp << std::endl;
        if((sys.current_time_ - it_expired_model->newest_measurement_time_stamp) > sys.params_.track_max_missed_detection_time_) {
            it_expired_model = sys.models_.erase(it_expired_model); // erases the current iterator and returns the next iterator      
        } else {
            ++it_expired_model;
        }
    }
 

    // Remove the worst ones until there are only the desired number of models or less
    while (sys.models_.size() > sys.params_.track_max_num_tracks_) {

        auto iter_to_remove = sys.models_.begin(); // The iterator will point to the model to be removed

        for(auto it = sys.models_.begin(); it != sys.models_.end(); ++it) {
            if(it->model_likelihood_ < iter_to_remove->model_likelihood_) {
                iter_to_remove = it;
            }
        }
#ifdef DEBUG_BUILD  // Had an issue with the begining iterator being the same as the end iterator when size wasn't zero. We believe the cause was due to not locking resources when threading in track initialization.
        if (iter_to_remove != sys.models_.end())
            sys.models_.erase(iter_to_remove);
        else {
            throw std::runtime_error("ModelManager::PruneModels: Iter of model to be removed is pointing to the end of the list. There is no model there!!");
        }
#else
    sys.models_.erase(iter_to_remove);
#endif

        
    }
}

//-------------------------------------------------------------------------------------------------------------------

template <typename _Model>
void ModelManager<_Model>::PruneConsensusSets(Sys& sys, const double expiration_time) {

    for(auto it = sys.models_.begin(); it != sys.models_.end(); ++it) {
        it->PruneConsensusSet(expiration_time);
    }
}

//-------------------------------------------------------------------------------------------------------------------

template <typename _Model>
void ModelManager<_Model>::PropagateModels(Sys& sys, const double dt) {
    for(auto it = sys.models_.begin(); it != sys.models_.end(); ++it) {
        it->PropagateModel(dt);
    }
}

//-------------------------------------------------------------------------------------------------------------------

template <typename _Model>
void ModelManager<_Model>::UpdateModels(Sys& sys) {
    for(auto it = sys.models_.begin(); it != sys.models_.end(); ++it) {
        it->UpdateModel(sys.source_container_, sys.params_);
    }
}

//-------------------------------------------------------------------------------------------------------------------

template <typename _Model>
void ModelManager<_Model>::MergeModels(Sys& sys) {

for (auto iter1 = sys.models_.begin(); iter1 != sys.models_.end(); ++iter1) {
    
    for(auto iter2 = std::next(iter1,1); iter2 != sys.models_.end();) {
        if (SimilarModels(sys, *iter1, *iter2)) {
            FuseModels(*iter1, *iter2);   // Fuse the models and store them as iter1
            iter2 = sys.models_.erase(iter2);      // Erase the element at iter2 and point iter2 to the next one
                                                   
        } else {
            ++iter2;
        }
    }
   
}

}

//-------------------------------------------------------------------------------------------------------------------

template <typename _Model>
void ModelManager<_Model>::RankModels(Sys& sys) {

    sys.good_models_.clear();

    for (auto iter = sys.models_.begin(); iter != sys.models_.end(); ++iter) {

        // See if the model is a good model
        if (iter->model_likelihood_ >= sys.params_.track_good_model_threshold_) { 

            // See if it needs a label
            if (iter->label_ == -1) {
                iter->label_ = sys.model_label_;
                sys.model_label_++;
            }

            sys.good_models_.push_back(&(*iter));

        }


    }

} 

//-------------------------------------------------------------------------------------------------------------------

template <typename _Model>
bool ModelManager<_Model>::SimilarModels(const Sys& sys, const Model& model1, const Model& model2) {

    bool similar = false;

    Eigen::Matrix<double,Model::cov_dim_,1> err_21 = Model::OMinus(model1, model2);
    typename Model::MatModelCov Jr_inv = Model::JrInv(err_21);
    typename Model::MatModelCov Jl_inv = Model::JlInv(err_21);
    typename Model::MatModelCov P_21;
    P_21 = Jl_inv*model2.err_cov_*Jl_inv.transpose() + Jr_inv*model1.err_cov_*Jr_inv.transpose();
    
    
    double d = err_21.transpose()*P_21.inverse()*err_21;

    if(sys.params_.track_similar_tracks_threshold_ > d) {
        similar = true;
    } 
    return similar;

}

//-------------------------------------------------------------------------------------------------------------------

template <typename _Model>
void ModelManager<_Model>::FuseModels(Model& model1, const Model& model2) {



//////////////////////
// Error covariance
//////////////////////


typename Model::MatModelCov P1_inv = model1.err_cov_.inverse();
typename Model::MatModelCov P2_inv = model2.err_cov_.inverse();
Eigen::Matrix<double,Model::cov_dim_,1> err_21 = Model::OMinus(model1, model2);

typename Model::MatModelCov Jr_inv = Model::JrInv(err_21);


typename Model::MatModelCov P_inv = P1_inv + Jr_inv.transpose()*P2_inv*Jr_inv;
typename Model::MatModelCov P = P_inv.inverse();


Eigen::Matrix<double,Model::cov_dim_,1> mu = -P*Jr_inv.transpose()*P2_inv*err_21;

typename Model::MatModelCov Jr = Model::Jr(mu);


model1.OPlusEQ(mu);
model1.err_cov_ = Jr*P*Jr.transpose();

// set the missed_detection_time to the lowest between the states
if (model1.newest_measurement_time_stamp < model2.newest_measurement_time_stamp)
    model1.newest_measurement_time_stamp = model2.newest_measurement_time_stamp;

// Set the label to the most recent label
if (model1.label_ == -1)
    model1.label_ = model2.label_;
else if (model1.label_ > model2.label_ && model2.label_ >=0) 
    model1.label_ = model2.label_;

// An accurate model_likelihood would require keeping a history of the number of associated measurements, probability of detection, etc which we dont do.
// thus we set the model likelihood to the best one. 
if (model1.model_likelihood_ < model2.model_likelihood_)
    model1.model_likelihood_ = model2.model_likelihood_;


for (auto iter = model2.cs_.consensus_set_.begin(); iter != model2.cs_.consensus_set_.end(); ++iter) {
    
    model1.cs_.AddMeasurementsToConsensusSetSameTimeStamp((*iter));

}
}



} // namespace rransac


#endif // RRANSAC_COMMON_MODEL_MANAGER_H_
