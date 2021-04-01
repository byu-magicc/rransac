#ifndef RRANSAC__COMMON__DATA_ASSOCIATION__MODEL_POLICIES__MODEL_PDF_POLICY_H_
#define RRANSAC__COMMON__DATA_ASSOCIATION__MODEL_POLICIES__MODEL_PDF_POLICY_H_
#pragma once


#include <math.h>
#include <numeric>
#include <vector>
#include "rransac/system.h"




namespace rransac
{
    
/**
 * \class ModelPDFPolicy
 * This policy associates measurements to tracks using the probabilistic data association (PDA) method. 
 * If a measurement is associated to a track, it is removed from System::new_meas_ and added to Model::new_assoc_meas_ using the method ModelBase::AddNewMeasurement. 
 * It must also update the model member variable Model::model_likelihood_update_info_ with the proper information.
 * The policy must expose the function static void PolicyDataAssociationModel which is called by the host class. TODO::clean up this class. The current version was for optimization purposes, but it's ugly. So clean it up.
 */ 

template<typename tModel>
class ModelPDFPolicy {

public:

typedef typename tModel::DataType DataType;   /**< The scalar object for the data. Ex. float, double, etc. */

/**
 * Implements the probabilistic data association method. * 
 * @param[in,out] sys The object that contains all of the data of RRANSAC. This includes the models. 
 */ 
static void  PolicyDataAssociationModel(System<tModel>& sys);


/**
 * This policy is used after RANSAC has found a hypothetical state estimate with sufficient inliers.
 * In chronological order, the inliers are added as new measurements to the model. Once added, the weights
 * and model update info is calculated for all the inliers of the same time step. Then the model is updated.
 * This process is repeated for every time step until all inliers are added. NOTE: This method is not compatible for measurements with non fixed measurement covariance.
 * @param[in,out] sys The object that contains all of the data of RRANSAC. This includes the models. 
 * @param model The model that is being filtered.
 */ 
static void CalculateMeasurmentAndLikelihoodDataPolicy(const System<tModel>& sys, tModel& model);

private: 

/**
 * Used with the filtering process to calculate the measurement likelihood, validation volume, and model likelihood update info
 * @param[in,out] sys The object that contains all of the data of RRANSAC. This includes the models. 
 * @param model The model that is being filtered.
 */ 
static void FilteringCalculateMeasStatistics(const System<tModel>& sys, tModel& model);

/**
 * Tries to associate the new measurements to tracks. A measurement is associated to a track is it is 
 * inside the validation region of the track. All of the associated measurements are given weights according
 * to the PDA algorithm to be later used in updating the track's state estimate. 
 * @param[in,out] sys The object that contains all of the data of RRANSAC. This includes the models. 
 * @param[in,out] source_produced_meas Not every source will produce measurements at every sensor scan. This parameter 
 *                              indicates which sources produced measurements during the current sensor scan.
 */ 
static void AssociateMeasurements(System<tModel>& sys, std::vector<bool>& source_produced_meas);

/**
 * During the track initialization process, measurement inliers are added to the track's ModelBase::new_assoc_meas_ member variable.
 * This method calculates the weights for those measurements. 
 * @param[in,out] sys The object that contains all of the data of RRANSAC. This includes the models. 
 * @param[in] model The track whose new associated measurements are being weighed.
 */ 
static void CalculateWeightsForModel(const System<tModel>& sys, tModel& model);

/**
 * Calls the function CalculateWeightsForModel for every track;
 * @param[in,out] sys The object that contains all of the data of RRANSAC. This includes the models. 
 */ 
static void CalculateWeights(System<tModel>& sys);

/**
 * Computes the geodesic distance between the measurement and the state of a track normalized by the innovation covariance.
 * @param[in] sources All of the measurement sources.
 * @param[in] meas The measurement whose distance from a track is being calculated.
 * @param[in] model The track whose distance from a measurement is being calculated.
 * @param[in] innovation_covariance The innovation covariance associated with the probability the measurement originated from the track. 
 * @return The geodesic distance between the measurement and the state of the track normalized by the innovation covariance. 
 */ 
static double GetDistance(const std::vector<typename tModel::Source>& sources, const Meas<DataType>& meas, const tModel& model, const Eigen::MatrixXd& innovation_covariance);

/**
 * Calculates the information needed to update the track's model likelihood. 
 * @param[in,out] sys The object that contains all of the data of RRANSAC. This includes the models. 
 * @param[in,out] source_produced_meas Not every source will produce measurements at every sensor scan. This parameter 
 *                              indicates which sources produced measurements during the current sensor scan.
 */ 
static void CalculateModelUpdateInfo(System<tModel>& sys, std::vector<bool>& source_produced_meas);

/**
 * Calculates the volume of the validation region. The volume of the validation region depends on the track and source.
 * @param[in,out] sys The object that contains all of the data of RRANSAC. This includes the models. 
 * @param[in] det_inn_cov_sqrt The square root of the innovation covariance's determinant. 
 * @param[in] source_index The index of the source. 
 */ 
static double GetVolume(const System<tModel>& sys, const double det_inn_cov_sqrt, const double source_index);

/**
 * Calculates the likelihood the measurement originated from a track.
 * @param[in] distance The normalized distance the measurement is from the model.
 * @param[in] dimensions The dimensions of the measurement.
 * @param[in] det_inn_cov_sqrt The square root of the determinant of the innovation covariance corresponding to a measurement source and track.
 * @return The likelihood the measurement originated from a track
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
void ModelPDFPolicy<tModel>::CalculateMeasurmentAndLikelihoodDataPolicy(const System<tModel>& sys, tModel& model) {

    FilteringCalculateMeasStatistics(sys,model);
    CalculateWeightsForModel(sys, model);
}

//-------------------------------------------------------------------------------------------------------------------------------------------------

template<typename tModel>
void ModelPDFPolicy<tModel>::FilteringCalculateMeasStatistics(const System<tModel>& sys, tModel& model) {

    Eigen::MatrixXd S; // Innovation covariance
    double det_S_sqrt = 0; 
    double det_S = 0;
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
        S = model.GetInnovationCovariance(sys.sources_, outer_iter->front().source_index);
        det_S = S.determinant();
        det_S_sqrt = sqrt(det_S);
        
        info.volume = GetVolume(sys, det_S_sqrt, src_index);

        for (auto inner_iter = outer_iter->begin(); inner_iter != outer_iter->end(); ++ inner_iter) {

            distance = GetDistance(sys.sources_, *inner_iter, model, S);
            inner_iter->likelihood = GetLikelihood(distance, S.rows(), det_S_sqrt); 
            info.num_assoc_meas++; // increment the number of associated measurements for the source

#ifdef DEBUG_BUILD
            if (isnan(inner_iter->likelihood)) {
                std::cerr << "Model PDF policy: measurement likelihood is NAN. Showing some values" << std::endl;
                std::cerr << "meas pose: "  << std::endl << inner_iter->pose << std::endl;
                std::cerr << "meas twist: " << std::endl << inner_iter->twist << std::endl;
                std::cerr << "distance: " << distance << std::endl;
                std::cerr << "S: " << std::endl << S << std::endl;
                std::cerr << "det S sqrt: " << det_S_sqrt << std::endl;
                std::cerr << "S.determinant(): " << S.determinant() << std::endl;
                
                std::cerr << "P.determinant(): " << model.err_cov_.determinant() << std::endl;
                std::cerr << "here: " << std::endl;
            }
#endif


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
    auto meas_iter = sys.new_meas_.begin();
    while(meas_iter != sys.new_meas_.end()) {
        bool meas_associated = false;

        // Indicates that a source produced a measurement
        source_produced_meas[meas_iter->source_index] = true;

        // Iterate through all of the models
        for(auto model_iter = sys.models_.begin(); model_iter != sys.models_.end(); ++ model_iter) {
            
            innovation_covariance = model_iter->GetInnovationCovariance(sys.sources_, meas_iter->source_index);
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
        } else {
            ++meas_iter;
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
        total_likelihood_ratio *= (source.params_.probability_of_detection_/ source.params_.spacial_density_of_false_meas_);

        // calculate the weights            
        double && denominator = 1.0 - source.params_.probability_of_detection_*source.params_.gate_probability_ + total_likelihood_ratio;
        for(auto inner_meas_iter = outer_meas_iter->begin(); inner_meas_iter != outer_meas_iter->end(); ++inner_meas_iter) {                

            inner_meas_iter->weight = inner_meas_iter->likelihood*source.params_.probability_of_detection_/ (source.params_.spacial_density_of_false_meas_*denominator);

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
    Eigen::MatrixXd innovation_covariance;
    double distance = 0;
    double det_inn_cov_sqrt = 0;
    

    for(auto model_iter = sys.models_.begin(); model_iter != sys.models_.end(); ++model_iter) {

        model_iter->model_likelihood_update_info_.clear();

        for(auto source_iter = sys.sources_.begin(); source_iter != sys.sources_.end(); ++source_iter) {

            if(source_produced_meas[source_iter->params_.source_index_]) {
                info.source_index = source_iter->params_.source_index_;
                if( source_iter->StateInsideSurveillanceRegion(model_iter->state_) ) {

                    info.in_local_surveillance_region = true;
                    info.num_assoc_meas = 0; // set default value
                    info.volume = 0;

                    // see if there are measurements with the source
                    for(auto outer_meas_iter = model_iter->new_assoc_meas_.begin(); outer_meas_iter != model_iter->new_assoc_meas_.end(); ++ outer_meas_iter) {
                        if (outer_meas_iter->begin()->source_index == source_iter->params_.source_index_) {
                            info.num_assoc_meas = outer_meas_iter->size();
                            info.volume = outer_meas_iter->begin()->vol;
                            break;
                        }
                    }

                    // Indicates that there was not a measurement pertaining to this source
                    if (info.volume == 0) {
                        innovation_covariance = model_iter->GetInnovationCovariance(sys.sources_, source_iter->params_.source_index_);
                        det_inn_cov_sqrt = sqrt(innovation_covariance.determinant());
                        info.volume = GetVolume(sys, det_inn_cov_sqrt, source_iter->params_.source_index_);

                    }


                } else {
                    info.in_local_surveillance_region = false;
                    info.num_assoc_meas = 0; // set default value
                    info.volume = 1;
                }

#ifdef DEBUG_BUILD
                if (info.volume <= 0) {
                    throw std::runtime_error("ModelPDFPolicy::CalculateModelUpdateInfo Volume of validation region is less than or equal to 0");
                }
#endif

                model_iter->model_likelihood_update_info_.push_back(info);
            }
        }
    }
}

//-------------------------------------------------------------------------------------------------------------------------------------------------

template<typename tModel>
double ModelPDFPolicy<tModel>::GetDistance(const std::vector<typename tModel::Source>& sources, const Meas<DataType>& meas, const tModel& model, const Eigen::MatrixXd& innovation_covariance) {
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
