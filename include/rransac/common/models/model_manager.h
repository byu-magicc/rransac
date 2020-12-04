#include <random>

#include "system.h"

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
* Prunes the consensus set for each model.
*/
static void PruneConsensusSets(System<tModel>& sys, const double expiration_time);

/**
* Propagates every model
*/
static void PropagateModels(System<tModel>& sys, const double dt);

/**
* Updates every model
*/
static void UpdateModels(System<tModel>& sys);



private:

/**
* Removes the models with the lowest model likelihood.
*/
static void PruneModels(System<tModel>& sys);

/**
* Looks for similar models and merges them.
*/
static void MergeModels(System<tModel>& sys);

/**
* Tests to see if two models are similar
* @return Returns true if the models are similar
*/
static bool SimilarModels(const System<tModel>& sys, const tModel& model1, const tModel& model2);

/**
 *  Fuse two models together using the sampled covariance intersection method
 */ 
tModel FuseTracks(const tModel& model1, const tModel& model2);



};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                            Definitions
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename tModel>
void ModelManager<tModel>::AddModel(System<tModel>& sys, const tModel& model) {

    sys.models_.push_back(model);

    if (sys.models_.size() > sys.params_.max_num_models_) {
    
        PruneModels(sys);
    }

}

//-------------------------------------------------------------------------------------------------------------------

template <typename tModel>
void ModelManager<tModel>::PruneModels(System<tModel>& sys) {

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
tModel ModelManager<tModel>::FuseTracks(const tModel& model1, const tModel& model2) {

typename tModel::State fused_state;

tModel::Mat P1_inv = model1.err_cov_.inverse();
tModel::Mat P2_inv = model2.err_cov_.inverse()
tModel::Mat P_inv = P1_inv + P2_inv;
tModel::Mat P = P_inv.inverse();
tModel::Mat P_sqrt = P.sqrt();

// Generate 100 random samples
std::vector<Eigen::Matrix<double,tModel::cov_dim_,1>> samples(100);

for (auto& sample : samples) {
    sample = P_sqrt*utilities::GaussianRandomGenerator(tModel::cov_dim_);
}

double r_max = -1;
double r_min = -1;
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

P = P/(0.5*(r_min + r_max));

fused_state = x1 oplus P*P2_inv*tModel::OMinus(model1, model2)

}




} // namespace rransac



