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
template <typename tModel>
class ModelManager {

public:

typedef tModel Model; /**< The object type of the track. */

/**
* Add a new model. If the number of models is greater than the max number of models, then
* the function PruneModels is called.
* @param[in,out] sys The object that contains all of the data of RRANSAC. Thus it contains all of the tracks. 
* @param[in] model The model to be added.
*/
static void AddModel(System<tModel>& sys, const tModel& model);

/**
* Propagates every track by calling ModelBase::PropagateModel method on each 
* track.
* @param[in,out] sys The object that contains all of the data of RRANSAC. Thus it contains all of the tracks. 
* @param[in] dt The amount of time the track needs to be propagated. 
*/
static void PropagateModels(System<tModel>& sys, const double dt);

/**
* Updates every track by calling ModelBase::UpdateModel method on each track.
* @param[in,out] sys The object that contains all of the data of RRANSAC. Thus it contains all of the tracks. 
*/
static void UpdateModels(System<tModel>& sys);


/**
* Updates the parameters of every track. 
* @param[in,out] sys The object that contains all of the data of RRANSAC. Thus it contains all of the tracks and the system parameters 
*/
static void SetModelParameters(System<tModel>& sys){ 
    for (auto iter = sys.models_.begin(); iter!=sys.models_.end(); ++iter) {
        iter->SetParameters(sys.params_);
    }
}

/**
 * Transforms all of the tracks and their consensus sets from the previous tracking frame to the current tracking frame. 
 * The consensus sets are only transformed if the flag Parameters::transform_consensus_set_ is set to true.
 * @param[in,out] sys The object that contains all of the data of RRANSAC. Thus it contains all of the tracks. 
 */ 
static void TransformModels(System<tModel>& sys){
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
static void ManageModels(System<tModel>& sys, const double expiration_time) {
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
static void PruneConsensusSets(System<tModel>& sys, const double expiration_time);


/**
* Looks for similar models and merges them. This method is similar to the track 2 track fusion method discussed in 
* Tracking and Data Fusion by Bar-Shalom 2011.
* @param[in,out] sys The object that contains all of the data of RRANSAC. Thus it contains all of the tracks. 
*/
static void MergeModels(System<tModel>& sys);

/**
 * Ranks the models according to their model likelihood. If the model likelihood is above the threshold Parameters::track_good_model_threshold_, then
 * it is considered a good model; otherwise, a poor model. Good models are added to System::good_models_. If there are more tracks than Parameters::track_max_num_tracks_, the tracks with
 * the lowest model likelihood are removed until there are only Parameters::track_max_num_tracks_ number of tracks. 
 * @param[in,out] sys The object that contains all of the data of RRANSAC. Thus it contains all of the tracks. 
 */
static void RankModels(System<tModel>& sys); 

/**
* If there are more tracks than Parameters::track_max_num_tracks_, the tracks with
* the lowest model likelihood are removed until there are only Parameters::track_max_num_tracks_ number of tracks.
* @param[in,out] sys The object that contains all of the data of RRANSAC. Thus it contains all of the tracks.  
*/
static void PruneModels(System<tModel>& sys);

/**
* Tests to see if two tracks are similar by weighing the geodesic distance between the track's states by their error covariances.
* @param[in] sys The object that contains all of the data of RRANSAC. Thus it contains all of the tracks. 
* @param[in] model1 One of the two tracks that will be compared.
* @param[in] model2 One of the two tracks that will be compared.
* @return Returns true if the tracks are similar
*/
static bool SimilarModels(const System<tModel>& sys, const tModel& model1, const tModel& model2);

/**
 * Fuse two tracks together using the sampled covariance intersection method. The method is described in "A no-loss covariance intersection algorithm
 * for track-to-track fusion" by Xin Tian 2010. The method also merges the the two consensus sets together. If one or both of the tracks are good tracks, the lowest 
 * label is kept because it corresponds to the first track that was promoted to a good track. We do not properly fuse the model likelihoods together. To properly fuse them 
 * would require keeping a history of their ModelBase::model_likelihood_update_info_, which could require a lot of data storage depending on how long the tracks have existed.
 * So we simplified the prosecces by assigning the largetest model likelihood to the fused track.
 * @param[in] model1 One of the two tracks that will be fused together.
 * @param[in] model2 One of the two tracks that will be fused together.
 * @return The fused track.
 */ 
static tModel FuseModels(const tModel& model1, const tModel& model2);



};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                            Definitions
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename tModel>
void ModelManager<tModel>::AddModel(System<tModel>& sys, const tModel& model) {

    sys.models_.push_back(model);
    sys.accumulative_number_of_tracks_++;

}

//-------------------------------------------------------------------------------------------------------------------

template <typename tModel>
void ModelManager<tModel>::PruneModels(System<tModel>& sys) {

    //
    // IDEA:: WE MIGHT WANT TO PRUNE MODELS WHOSE LIKELIHOOD FALL BELOW A THRESHOLD
    // THE REASON BEHIND THIS IS THAT A MODEL MIGHT NOT GET A MEASUREMENT FOR A WHILE
    // BECUASE IT IS HIDDEN OR SOMETHING, BUT BY THE TIME IT DOES RECEIVE A NEW MEASUREMENT
    // ITS LIKELIHOOD COULD BE SO LOW THAT IT TAKES A WHILE TO HAVE A POSITIVE LIKELIHOOD AGAIN

    // Remove the models that have not received a measurement for a while
    for(auto it = sys.models_.begin(); it != sys.models_.end(); ++it) {
        if(it->missed_detection_time_ > sys.params_.track_max_missed_detection_time_) {
            auto it_copy = it;
            --it;                       // We are going to remove it so we need to back up to the previous
            sys.models_.erase(it_copy);
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
        sys.models_.erase(iter_to_remove);
    }
}

//-------------------------------------------------------------------------------------------------------------------

template <typename tModel>
void ModelManager<tModel>::PruneConsensusSets(System<tModel>& sys, const double expiration_time) {

    for(auto it = sys.models_.begin(); it != sys.models_.end(); ++it) {
        it->PruneConsensusSet(expiration_time);
    }
}

//-------------------------------------------------------------------------------------------------------------------

template <typename tModel>
void ModelManager<tModel>::PropagateModels(System<tModel>& sys, const double dt) {
    for(auto it = sys.models_.begin(); it != sys.models_.end(); ++it) {
        it->PropagateModel(dt);
    }
}

//-------------------------------------------------------------------------------------------------------------------

template <typename tModel>
void ModelManager<tModel>::UpdateModels(System<tModel>& sys) {
    for(auto it = sys.models_.begin(); it != sys.models_.end(); ++it) {
        it->UpdateModel(sys.sources_, sys.params_);
    }
}

//-------------------------------------------------------------------------------------------------------------------

template <typename tModel>
void ModelManager<tModel>::MergeModels(System<tModel>& sys) {

for (auto iter1 = sys.models_.begin(); iter1 != sys.models_.end(); ++iter1) {

    for (auto iter2 = std::next(iter1,1); iter2 != sys.models_.end(); ++iter2) {

        if (SimilarModels(sys, *iter1, *iter2)) {
            *iter1 = FuseModels(*iter1, *iter2);   // Fuse the models and store them as iter1
            iter2 = sys.models_.erase(iter2);      // Erase the element at iter2 and point iter2 to the next one
            --iter2;                               // Bring iter2 back one b/c the for loop will increment it.
        }


    }
}

}

//-------------------------------------------------------------------------------------------------------------------

template <typename tModel>
void ModelManager<tModel>::RankModels(System<tModel>& sys) {

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

template <typename tModel>
bool ModelManager<tModel>::SimilarModels(const System<tModel>& sys, const tModel& model1, const tModel& model2) {

    bool similar = false;

    Eigen::Matrix<double,tModel::cov_dim_,1> err = tModel::OMinus(model1, model2);

    typename tModel::Mat T = model1.err_cov_ + model2.err_cov_;
    double d = err.transpose()*T.inverse()*err;

    if(sys.params_.track_similar_tracks_threshold_ > d) {
        similar = true;
    } 
    return similar;

}

//-------------------------------------------------------------------------------------------------------------------

template <typename tModel>
tModel ModelManager<tModel>::FuseModels(const tModel& model1, const tModel& model2) {

tModel fused_model = model1;

//////////////////////
// Error covariance
//////////////////////
typename tModel::Mat P1_inv = model1.err_cov_.inverse();
typename tModel::Mat P2_inv = model2.err_cov_.inverse();
typename tModel::Mat P_inv = P1_inv + P2_inv;
typename tModel::Mat P = P_inv.inverse();
const typename tModel::Mat P_sqrt = P.sqrt();


// The sample covariance intersection method needs 100 samples
std::vector<Eigen::Matrix<double,tModel::cov_dim_,1>> samples(100);
for (auto& sample : samples) {
    sample = P_sqrt*utilities::GaussianRandomGenerator(tModel::cov_dim_);
}

double r_max = -1;
double r_min = 1e10;
double r_candidate;
for (auto& sample : samples) {

    if (sample.norm() != 0 ) {
        double tmp1 = sample.transpose()*P_inv*sample;
        double tmp2 = sample.transpose()*P1_inv*sample;
        double tmp3 = sample.transpose()*P2_inv*sample;

        if (tmp2 > tmp3) {
            r_candidate = tmp1/tmp2;
        }  else {
            r_candidate = tmp1/tmp3;
        }   

        if (r_candidate > r_max)
            r_max = r_candidate;
        if (r_candidate < r_min) 
            r_min = r_candidate;
    }


}

// Scale the error covariance
P = P/(0.5*(r_min + r_max));

fused_model.OPlusEQ(P*P2_inv*tModel::OMinus(model2, model1));


fused_model.err_cov_ = P;

// set the missed_detection_time to the lowest between the states
if (model1.missed_detection_time_ > model2.missed_detection_time_)
    fused_model.missed_detection_time_ = model2.missed_detection_time_;

// Set the label to the most recent label
if (model1.label_ == -1)
    fused_model.label_ = model2.label_;
else if (model1.label_ > model2.label_ && model2.label_ >=0) 
    fused_model.label_ = model2.label_;

// An accurate model_likelihood would require keeping a history of the number of associated measurements, probability of detection, etc which we dont do.
// thus we set the model likelihood to the best one. 
if (model1.model_likelihood_ < model2.model_likelihood_)
    fused_model.model_likelihood_ = model2.model_likelihood_;



fused_model.cs_ = fused_model.cs_.MergeConsensusSets(model1.cs_, model2.cs_);
return fused_model;
}




} // namespace rransac


#endif // RRANSAC_COMMON_MODEL_MANAGER_H_
