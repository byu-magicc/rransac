#ifndef RRANSAC__COMMON__DATA_ASSOCIATION__MODEL_POLICIES__MODEL_PDF_POLICY_H_
#define RRANSAC__COMMON__DATA_ASSOCIATION__MODEL_POLICIES__MODEL_PDF_POLICY_H_

#include <math.h>
#include <numeric>
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

private: 

static void AssociateMeasurements(System<tModel>& sys);

static void CalculateWeights(System<tModel>& sys);

static void CalculateModelUpdateInfo(System<tModel>& sys);

static bool InValidationRegion(const Meas& meas, const tModel& model, const Eigen::MatrixXd& innovation_covariance, double& distance);

static double GetVolume(const System<tModel>& sys, const double det_inn_cov_sqrt, const double source_index);

};




//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                            Definitions
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template<typename tModel>
void DataTreeClusterAssociationPolicy<tModel>::PolicyDataAssociationModel(System<tModel>& sys) {
    AssociateMeasurements(sys);
    CalculateModelUpdateInfo(sys);
    CalculateWeights(sys);
}

//-------------------------------------------------------------------------------------------------------------------------------------------------

template<typename tModel>
void DataTreeClusterAssociationPolicy<tModel>::AssociateMeasurements(System<tModel>& sys) {

    // Clear the measurements
    for(auto model_iter = sys.models_.begin(); model_iter != sys.models_.end(); ++ model_iter) {
        model_iter->new_assoc_meas_.clear();
    }

    // Add the new measurements
    for(auto outer_meas_iter = sys.new_meas_.begin(); outer_meas_iter != sys.new_meas_.end(); ++outer_meas_iter) {
        for(auto inner_meas_iter = outer_meas_iter.begin(); inner_meas_iter != outer_meas_iter->end(); ++inner_meas_iter) {

            bool meas_associated = false;

            // Iterate through all of the models
            for(auto model_iter = sys.models_.begin(); model_iter != sys.models_.end(); ++ model_iter) {
                
                const Eigen::MatrixXd& innovation_covariance = model_iter->GetInnovationCovariance(*inner_iter);
                double distance = 0;
                double det_inn_cov_sqrt = sqrt(innovation_covariance.determinant());

                if(InValidationRegion(*inner_meas_iter, *model_iter, innovation_covariance,distance)) {
                    meas_associated = true;
                    inner_meas_iter->likelihood = exp(-distance/2.0)/(pow(2*M_PI,innovation_covariance.rows/2.0)*det_inn_cov_sqrt);
                    inner_meas_iter->vol = GetVolume(sys, det_inn_cov_sqrt, inner_meas_iter->source_index);
                    model_iter->AddNewMeasurement(inner_meas_iter);
                }
            }

            // Remove the measurement if it was associated
            if(meas_associated) {
                inner_meas_iter = outer_meas_iter->erase(inner_meas_iter);
                --inner_meas_iter;
            }
        }
    }

}

//-------------------------------------------------------------------------------------------------------------------------------------------------

template<typename tModel>
void DataTreeClusterAssociationPolicy<tModel>::CalculateWeights(System<tModel>& sys) {

    for(auto model_iter = sys.models_.begin(); model_iter != sys.models_.end(); ++model_iter) {
        for(auto outer_meas_iter = model_iter->new_assoc_meas_.begin(); outer_meas_iter != model_iter->new_assoc_meas_.end(); ++outer_meas_iter) {
            std::vector<double> likelihood_ratios(outer_meas_iter->size());
            double total_likelihood_ratio = 0;
            auto vec_iter = likelihood_ratios.begin();

            typename tModel::Source& source = sys.sources_[outer_meas_iter->begin()->source_index];

            // Calculate the partial likelihood ratios
            for(auto inner_meas_iter = outer_meas_iter->begin(); inner_meas_iter != outer_meas_iter->end(); ++inner_meas_iter) {                

                *vec_iter = inner_meas_iter->likelihood*source.params_.probability_of_detection_/ source.params_.expected_num_false_meas_;
                ++vec_iter;
            }

            // Get the total likelihood ratio
            double total_likelihood_ratio = std::accumulate(likelihood_ratios.begin(), likelihood_ratios.end(), 0.0);

            // calculate the weights
            vec_iter = likelihood_ratios.begin();
            double && denominator = 1 - source.params_.probability_of_detection_*source.params_.gate_probability_ + total_likelihood_ratio;
            for(auto inner_meas_iter = outer_meas_iter->begin(); inner_meas_iter != outer_meas_iter->end(); ++inner_meas_iter) {                

                inner_meas_iter->weight = *vec_iter/denominator;

                ++vec_iter;
            }


        }
    }

}

//-------------------------------------------------------------------------------------------------------------------------------------------------

template<typename tModel>
void DataTreeClusterAssociationPolicy<tModel>::CalculateModelUpdateInfo(System<tModel>& sys) {

    typename tModel::ModelLikelihoodUpdateInfo info;
    info.in_local_surveillance_region = true;

    for(auto model_iter = sys.models_.begin(); model_iter != sys.models_.end(); ++model_iter) {

        model_iter->model_likelihood_update_info_.clear();

        for(auto source_iter = sys.sources_.begin(); source_iter != sys.sources_.end(); ++ source_iter) {

            if( source_iter->StateInsideSurveillanceRegion(model_iter->state_) ) {
                info.source_index = source_iter->params_.source_index_;

                // Set default values
                info.num_assoc_meas = 0;
                info.source_index = 0;
                info.volume = 1;
                
                // see if there are measurements with the source
                for(auto outer_meas_iter = model_iter->new_assoc_meas_.begin(); outer_meas_iter != model_iter->new_assoc_meas_.end(); ++ outer_meas_iter) {
                    if (outer_meas_iter->begin()->source_index == source_iter->params_.source_index_) {
                        info.num_assoc_meas = outer_meas_iter->size();
                        info.vol = outer_meas_iter->begin()->vol;
                        break;
                    }
                }


            }
        }
    }
}

//-------------------------------------------------------------------------------------------------------------------------------------------------

template<typename tModel>
bool DataTreeClusterAssociationPolicy<tModel>::InValidationRegion(const Meas& meas, const tModel& model, const Eigen::MatrixXd& innovation_covariance, double& distance) {

    Eigen::MatrixXd err = (*(model.sources_))[meas.source_index].OMinus(meas, (*(model.sources_))[meas.source_index].GetEstMeas(model.state_));

    distance = err*innovation_covariance.inverse()*err.transpose();

    if(distance < (*(model.sources_))[meas.source_index].params_.gate_threshold_) {
        return true;
    } else {
        return false;
    }

}

//-------------------------------------------------------------------------------------------------------------------------------------------------

template<typename tModel>
double DataTreeClusterAssociationPolicy<tModel>::GetVolume(const System<tModel>& sys, const double det_inn_cov_sqrt, const double source_index) {
    return sys.sources_[source_index].params_.vol_unit_hypershpere_ * sys.sources_[source_index].params_.gate_threshold_sqrt_*det_inn_cov_sqrt;
}


//-------------------------------------------------------------------------------------------------------------------------------------------------



 } // namespace rransac

#endif //RRANSAC__COMMON__DATA_ASSOCIATION__MODEL_POLICIES__MODEL_PDF_POLICY_H_
