#include <random>
#include <unsupported/Eigen/MatrixFunctions>

#include "system.h"
#include "common/utilities.h"

namespace rransac
{
    
template <typename tModel>
class ModelManager {

public:

typedef tModel Model;

/**
* Add a new model. If the number of models is greater than the max number of models, then
* the function PruneModels is called.
*/
static void AddModel(System<tModel>& sys, const tModel& model);

/**
* Propagates every model
*/
static void PropagateModels(System<tModel>& sys, const double dt);

/**
 * Updates the models, prune the consensus set, merge models, ranks models, and prunes models
 */
static void ManageModels(System<tModel>& sys, const double expiration_time) {
    UpdateModels(sys);
    PruneConsensusSets(sys, expiration_time);
    MergeModels(sys);
    PruneModels(sys);
    RankAndPruneModels(sys);
}



private:

/**
* Prunes the consensus set for each model.
*/
static void PruneConsensusSets(System<tModel>& sys, const double expiration_time);

/**
* Updates every model
*/
static void UpdateModels(System<tModel>& sys);

/**
* Looks for similar models and merges them.
*/
static void MergeModels(System<tModel>& sys);

/**
 * Promotes models to good models, demotes good models to poor models, and prunes models
 */
static void RankModels(System<tModel>& sys); 

/**
* Removes the models with the lowest model likelihood.
*/
static void PruneModels(System<tModel>& sys);

/**
* Tests to see if two models are similar
* @return Returns true if the models are similar
*/
static bool SimilarModels(const System<tModel>& sys, const tModel& model1, const tModel& model2);

/**
 *  Fuse two models together using the sampled covariance intersection method
 */ 
static tModel FuseModels(const tModel& model1, const tModel& model2);



};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                            Definitions
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename tModel>
void ModelManager<tModel>::AddModel(System<tModel>& sys, const tModel& model) {

    sys.models_.push_back(model);

}

//-------------------------------------------------------------------------------------------------------------------

template <typename tModel>
void ModelManager<tModel>::PruneModels(System<tModel>& sys) {

    // Remove the models that have not received a measurement for a while
    for(auto it = sys.models_.begin(); it != sys.models_.end(); ++it) {
        if(it->missed_detection_time_ > sys.params_.max_missed_detection_time_) {
            auto it_copy = it;
            --it;                       // We are going to remove it so we need to back up to the previous
            sys.models_.erase(it_copy);
        }
    }

    // Remove the worst ones until there are only the desired number of models or less
    while (sys.models_.size() > sys.params_.max_num_models_) {

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
        it->UpdateModel(sys.params_);
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
        if (iter->model_likelihood_ >= sys.params_.good_model_threshold_) { 

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

    if(sys.params_.similar_tracks_threshold_ > d) {
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

// Update the state
fused_model.state_.OPlusEQ(P*P2_inv*tModel::OMinus(model2, model1));


fused_model.err_cov_ = P;

// set the missed_detection_time to the lowest between the states
if (model1.missed_detection_time_ > model2.missed_detection_time_)
    fused_model.missed_detection_time_ = model2.missed_detection_time_;

// Set the label to the most recent label
if (model1.label_ > model2.label_ && model2.label_ >=0) {
    fused_model.label_ = model2.label_;

// An accurate model_likelihood would require keeping a history of the number of associated measurements, probability of detection, etc which we dont do.
// thus we set the model likelihood to the best one. 
if (model1.model_likelihood_ < model2.model_likelihood_)
    fused_model.model_likelihood_ = model2.model_likelihood_;

}

fused_model.cs_ = fused_model.cs_.MergeConsensusSets(model1.cs_, model2.cs_);
return fused_model;
}




} // namespace rransac



