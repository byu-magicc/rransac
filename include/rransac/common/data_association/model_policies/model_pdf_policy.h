#ifndef RRANSAC__COMMON__DATA_ASSOCIATION__MODEL_POLICIES__MODEL_PDF_POLICY_H_
#define RRANSAC__COMMON__DATA_ASSOCIATION__MODEL_POLICIES__MODEL_PDF_POLICY_H_

#include <math.h>
#include <numeric>
#include <vector>
#include "system.h"




namespace rransac
{
    
/**
 * \class
 * This policy associates measurements to models using the probabilistic data association (PDA) method. 
 * If a measurement is associated to a model, it is removed from System<Model>::new_meas_ and added to Model<State>::new_assoc_meas_. 
 * It must also update the model member variable Model<State>::model_likelihood_update_info_ with the proper information.
 * The policy must expose the function static void PolicyDataAssociationModel(System<Model>& sys) which is called by the host class.
 */ 

template<typename tModel>
class ModelPDFPolicy {

public:

static void  PolicyDataAssociationModel(System<tModel>& sys);


/**
 * This policy is used after RANSAC has found a hypothetical state estimate with sufficient inliers.
 * In chronological order, the inliers are added as new measurements to the model. Once added, the weights
 * and model update info is calculated for all the inliers of the same time step. Then the model is updated.
 * This process is repeated for every time step until all inliers are added. NOTE: This method is not compatible for measurements with non fixed measurement covariance.
 * @param model The model that is being filtered.
 */ 
static void PolicyModelFiltering(const System<tModel>& sys, tModel& model);

private: 

/**
 * Used with the filtering process to calculate the measurement likelihood, validation volume, and model likelihood update info
 * 
 */ 
static void FilteringCalculateMeasStatistics(const System<tModel>& sys, tModel& model);

static void AssociateMeasurements(System<tModel>& sys, std::vector<bool>& source_produced_meas);

static void CalculateWeightsForModel(const System<tModel>& sys, tModel& model);

/**
 * Calls the function CalculateWeightsForModel for every model;
 */ 
static void CalculateWeights(System<tModel>& sys);


static double GetDistance(const std::vector<typename tModel::Source>& sources, const Meas& meas, const tModel& model, const Eigen::MatrixXd& innovation_covariance);

static void CalculateModelUpdateInfo(System<tModel>& sys, std::vector<bool>& source_produced_meas);

static double GetVolume(const System<tModel>& sys, const double det_inn_cov_sqrt, const double source_index);

/**
 * Calculates a measurement likelihood
 * @param distance The normalized distance the measurement is from the model
 * @param dimensions The dimensions of the measurement
 * @param det_inn_cov_sqrt The square root of the determinant of the innovation
 */ 
static double GetLikelihood(const double distance, const double dimensions, const double det_inn_cov_sqrt);

};




//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                            Definitions
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template<typename tModel>
void ModelPDFPolicy<tModel>::PolicyDataAssociationModel(System<tModel>& sys) {

    std::vector<bool> source_produced_meas; /**< Flag to indicate if a source produced a measurement*/

    AssociateMeasurements(sys, source_produced_meas);
    CalculateModelUpdateInfo(sys, source_produced_meas);
    CalculateWeights(sys);
}

//-------------------------------------------------------------------------------------------------------------------------------------------------

template<typename tModel>
void ModelPDFPolicy<tModel>::PolicyModelFiltering(const System<tModel>& sys, tModel& model) {

    FilteringCalculateMeasStatistics(sys,model);
    CalculateWeightsForModel(sys, model);
}

//-------------------------------------------------------------------------------------------------------------------------------------------------

template<typename tModel>
void ModelPDFPolicy<tModel>::FilteringCalculateMeasStatistics(const System<tModel>& sys, tModel& model) {

    Eigen::MatrixXd S; // Innovation covariance
    double det_S_sqrt = 0; 
    double distance = 0;
    int src_index = 0;

    ModelLikelihoodUpdateInfo info;

    model.model_likelihood_update_info_.clear();

    // Calculate measurement statistics. The outer iterator iterates through different sources. The inner iterator
    // iterates through different measurements of the same source
    for (auto outer_iter = model.new_assoc_meas_.begin(); outer_iter != model.new_assoc_meas_.end(); ++outer_iter) {

        src_index = outer_iter->front().source_index;
        info.num_assoc_meas = 0;
        info.source_index = src_index;
        info.in_local_surveillance_region = true;

        // The innovation covariance will be the same for all measurements of the same measurement source
        S = model.GetInnovationCovariance(sys.sources_, outer_iter->front());
        det_S_sqrt = sqrt(S.determinant());
        info.volume = GetVolume(sys, det_S_sqrt, src_index);

        for (auto inner_iter = outer_iter->begin(); inner_iter != outer_iter->end(); ++ inner_iter) {

            distance = GetDistance(sys.sources_, *inner_iter, model, S);
            inner_iter->likelihood = GetLikelihood(distance, S.rows(), det_S_sqrt); 

            info.num_assoc_meas++; // increment the number of associated measurements for the source

        }

        model.model_likelihood_update_info_.push_back(info);
    }

    

}

//-------------------------------------------------------------------------------------------------------------------------------------------------

template<typename tModel>
void ModelPDFPolicy<tModel>::AssociateMeasurements(System<tModel>& sys, std::vector<bool>& source_produced_meas) {

    // Initialize objects
    Eigen::MatrixXd innovation_covariance;
    double distance = 0;
    double det_inn_cov_sqrt = 0;

    // Reset the source produced meas flag
    source_produced_meas.clear();
    source_produced_meas.resize(sys.sources_.size(), false);

    // Clear the measurements
    for(auto model_iter = sys.models_.begin(); model_iter != sys.models_.end(); ++ model_iter) {
        model_iter->new_assoc_meas_.clear();
    }

    // Add the new measurements
    for(auto meas_iter = sys.new_meas_.begin(); meas_iter != sys.new_meas_.end(); ++meas_iter) {
        bool meas_associated = false;

        // Indicates that a source produced a measurement
        source_produced_meas[meas_iter->source_index] = true;

        // Iterate through all of the models
        for(auto model_iter = sys.models_.begin(); model_iter != sys.models_.end(); ++ model_iter) {
            
            innovation_covariance = model_iter->GetInnovationCovariance(sys.sources_, *meas_iter);
            distance = GetDistance(sys.sources_, *meas_iter, *model_iter, innovation_covariance);
            det_inn_cov_sqrt = sqrt(innovation_covariance.determinant());

            // In validation region
            if (distance <= sys.sources_[meas_iter->source_index].params_.gate_threshold_) {
                meas_associated = true;
                meas_iter->likelihood = GetLikelihood(distance, innovation_covariance.rows(), det_inn_cov_sqrt); 
                meas_iter->vol = GetVolume(sys, det_inn_cov_sqrt, meas_iter->source_index);
                model_iter->AddNewMeasurement(*meas_iter);
            }
        }

        // Remove the measurement if it was associated
        if(meas_associated) {
            meas_iter = sys.new_meas_.erase(meas_iter);
            --meas_iter;
        }
    
    }

}

//-------------------------------------------------------------------------------------------------------------------------------------------------

template<typename tModel>
void ModelPDFPolicy<tModel>::CalculateWeightsForModel(const System<tModel>& sys, tModel& model) {
    for(auto outer_meas_iter = model.new_assoc_meas_.begin(); outer_meas_iter != model.new_assoc_meas_.end(); ++outer_meas_iter) {
        double total_likelihood_ratio = 0;

        const typename tModel::Source& source = sys.sources_[outer_meas_iter->begin()->source_index];

        // Calculate the partial likelihood ratios
        for(auto inner_meas_iter = outer_meas_iter->begin(); inner_meas_iter != outer_meas_iter->end(); ++inner_meas_iter) {                

            total_likelihood_ratio += inner_meas_iter->likelihood;
        }

        // Get the total likelihood ratio
        total_likelihood_ratio *= (source.params_.probability_of_detection_/ source.params_.expected_num_false_meas_);

        // calculate the weights            
        double && denominator = 1.0 - source.params_.probability_of_detection_*source.params_.gate_probability_ + total_likelihood_ratio;
        for(auto inner_meas_iter = outer_meas_iter->begin(); inner_meas_iter != outer_meas_iter->end(); ++inner_meas_iter) {                

            inner_meas_iter->weight = inner_meas_iter->likelihood*source.params_.probability_of_detection_/ (source.params_.expected_num_false_meas_*denominator);

        }
    }
}

//-------------------------------------------------------------------------------------------------------------------------------------------------

template<typename tModel>
void ModelPDFPolicy<tModel>::CalculateWeights(System<tModel>& sys) {

    for(auto model_iter = sys.models_.begin(); model_iter != sys.models_.end(); ++model_iter) {
        CalculateWeightsForModel(sys, *model_iter);
    }
        

}

//-------------------------------------------------------------------------------------------------------------------------------------------------

template<typename tModel>
void ModelPDFPolicy<tModel>::CalculateModelUpdateInfo(System<tModel>& sys, std::vector<bool>& source_produced_meas) {

    ModelLikelihoodUpdateInfo info;
    

    for(auto model_iter = sys.models_.begin(); model_iter != sys.models_.end(); ++model_iter) {

        model_iter->model_likelihood_update_info_.clear();

        for(auto source_iter = sys.sources_.begin(); source_iter != sys.sources_.end(); ++source_iter) {

            if(source_produced_meas[source_iter->params_.source_index_]) {
                info.source_index = source_iter->params_.source_index_;
                if( source_iter->StateInsideSurveillanceRegion(model_iter->state_) ) {

                    info.in_local_surveillance_region = true;
                    info.num_assoc_meas = 0; // set default value
                    info.volume = 1;

                    // see if there are measurements with the source
                    for(auto outer_meas_iter = model_iter->new_assoc_meas_.begin(); outer_meas_iter != model_iter->new_assoc_meas_.end(); ++ outer_meas_iter) {
                        if (outer_meas_iter->begin()->source_index == source_iter->params_.source_index_) {
                            info.num_assoc_meas = outer_meas_iter->size();
                            info.volume = outer_meas_iter->begin()->vol;
                            break;
                        }
                    }


                } else {
                    info.in_local_surveillance_region = false;
                    info.num_assoc_meas = 0; // set default value
                    info.volume = 1;
                }

                model_iter->model_likelihood_update_info_.push_back(info);
            }
        }
    }
}

//-------------------------------------------------------------------------------------------------------------------------------------------------

template<typename tModel>
double ModelPDFPolicy<tModel>::GetDistance(const std::vector<typename tModel::Source>& sources, const Meas& meas, const tModel& model, const Eigen::MatrixXd& innovation_covariance) {
    Eigen::MatrixXd err = sources[meas.source_index].OMinus(meas, sources[meas.source_index].GetEstMeas(model.state_));
    return (err.transpose()*innovation_covariance.inverse()*err)(0,0);
}



//-------------------------------------------------------------------------------------------------------------------------------------------------

template<typename tModel>
double ModelPDFPolicy<tModel>::GetVolume(const System<tModel>& sys, const double det_inn_cov_sqrt, const double source_index) {
    return sys.sources_[source_index].params_.vol_unit_hypershpere_ * sys.sources_[source_index].params_.gate_threshold_sqrt_*det_inn_cov_sqrt;
}


//-------------------------------------------------------------------------------------------------------------------------------------------------

template<typename tModel>
double ModelPDFPolicy<tModel>::GetLikelihood(const double distance, const double dimensions, const double det_inn_cov_sqrt) {
    return exp(-distance/2.0)/(pow(2.0*M_PI,dimensions/2.0)*det_inn_cov_sqrt);
}

 } // namespace rransac

#endif //RRANSAC__COMMON__DATA_ASSOCIATION__MODEL_POLICIES__MODEL_PDF_POLICY_H_
