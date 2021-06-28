#ifndef RRANSAC_COMMON_DATA_ASSOCIATION_DATA_ASSOCIATION_HOST_H_
#define RRANSAC_COMMON_DATA_ASSOCIATION_DATA_ASSOCIATION_HOST_H_
#pragma once

#include <vector>

#include "rransac/system.h"

namespace rransac
{

struct DataAssociationInfo {

    std::vector<bool> source_produced_measurements_;                                         /**< When a sensor scan occurs one or more sources can produce measurements. This flag indicates if a source 
                                                                                                  during the latest sensor scan produced measurements. */

};

/**
 * \class DataAssociationHost
 * 
 * Associates the new measurements System::new_measurements_ to the tracks, clusters, and data tree. 
 * This class is host to three policies. The first policy is the validation region policy and it returns true if a measurement falls inside a track's validation region. The 
 * second ploicy is the track's likelihood policy and it is responsible
 * for updating the track's likelihood. The last policy is the measurement weight policy and is responsible for calculating weights for each measurement associated to a model.
 *  
 * The validation region policy takes in two arguments: a measurement and a track. It returns true if the measurement falls inside of the track's validation region; otherwise false. 
 *
 * The track likelihood update policy uses the new measurements associated with a track to update the track's likelihood. The track's likelihood is the probability that the target which the track represents
 * exists.
 *
 * The measurement weight policy assigns a weight to each track associated measurement that is used to update the state estimate and error covariance.
 *
 */ 
template<typename tModel, template<class> typename tValidationRegionPolicy, template<class> typename tUpdateTrackLikelihoodPolicy, template<class> typename tMeasurementWeightPolicy>
class DataAssociationHost : public tValidationRegionPolicy<tModel>, tUpdateTrackLikelihoodPolicy<tModel>, tMeasurementWeightPolicy<tModel>{

// Some policies are only compatible with others. These assertions enforce compatability. 
static_assert(tMeasurementWeightPolicy<tModel>::template CompatiblityCheck<tValidationRegionPolicy,tUpdateTrackLikelihoodPolicy>::value,"MW_IPDAFPolicy::CompatiblityCheck The policy MW_IPDAFPolicy is only compatible with the update track likelihood policy TLI_IPDAFPolicy." );



public:

    /**
     * Associates new measurements to the tracks, clusters and the data tree. 
     * @param[in,out] sys The object that contains all of the data of RRANSAC. This includes the new measurements, tracks, clusters and data tree. 
     */ 
    void AssociateNewMeasurements(System<tModel>& sys) {
        SourceProducedMeasurement(sys, this->data_association_info_);                    // Determine which sources produced measurements
        ResetModelLikelihoodAndNewMeasurements(sys);  
        DataAssociation(sys,this->data_association_info_);
        UpdateTrackLikelihood(sys,this->data_association_info_);
        CalculateMeasurementWeight(sys,this->data_association_info_);
        

        if(sys.new_meas_.size() != 0)   
            throw std::runtime_error("DataAssociationHost::AssociateNewMeasurements: All new measurements should have been copied to a model, cluster, or data tree and removed from System<Model>::new_meas_");
    }

private:

    /**
     * Associates measurements to at least one track or the data tree. If a measurement falls within the validation region of a track, then it is associated with the track; otherwise it is added
     * to the data tree. The validation region is determined by the policy. 
     * @param[in,out] sys The object that contains all of the data of RRANSAC. This includes the new measurements, tracks, clusters and data tree. 
     */ 
    void DataAssociation(System<tModel>& sys, DataAssociationInfo& info);

    /**
     *  Determines if the measurement falls inside the validation region of the track according to the policy. 
     * @param[in] sys The object that contains all of the data of RRANSAC. This includes the new measurements, tracks, clusters and data tree. 
     * @param[in] meas The measurement
     * @param[in] track The track. It is not a constant reference so be careful!! 
     */ 
    bool InValidationRegion(const System<tModel>& sys, const Meas<typename tModel::DataType>& meas, tModel& track) {
        return DataAssociationHost::PolicyInValidationRegion(sys, meas,track);
    }

    /** 
     * Uses the new associated measurements to update the track's likelihood. 
     * @param[in,out] sys The object that contains all of the data of RRANSAC. This includes the new measurements, tracks, clusters and data tree. 
     */
    void UpdateTrackLikelihood(System<tModel>& sys, DataAssociationInfo& info ) {
        DataAssociationHost::PolicyUpdateTrackLikelihood(sys,info);
    }

    /** 
     * Calculates the weights for each track associated measurement. 
     * @param[in,out] sys The object that contains all of the data of RRANSAC. This includes the new measurements, tracks, clusters and data tree. 
     */
    void CalculateMeasurementWeight(System<tModel>& sys, DataAssociationInfo& info) {
        DataAssociationHost::PolicyCalculateMeasurementWeight(sys,info);
    }

    /**
     * Sets the member vairable SourceBase::source_produced_measurements_ to true if the source produced a measurement this sensor scan; otherwise false.
     * @param[in,out] sys The object that contains all of the data of RRANSAC. This includes the new measurements, tracks, clusters and data tree. 
     */
    void SourceProducedMeasurement(System<tModel>& sys, DataAssociationInfo& info);

    /** 
     * Ensures that none of the tracks have any new associated measurements in ModelBase::new_assoc_meas_, and it resets the
     * model likelihood update info.
     */ 
    void ResetModelLikelihoodAndNewMeasurements(System<tModel>& sys);

    DataAssociationInfo data_association_info_;

};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                            Definitions
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//------------------------------------------------------------------------------------------------------------------------------------

template<typename tModel, template<class> typename tValidationRegionPolicy, template<class> typename tTrackLikelihoodUpdatePolicy, template<class> typename tMeasurementWeightPolicy>
void DataAssociationHost<tModel,tValidationRegionPolicy,tTrackLikelihoodUpdatePolicy,tMeasurementWeightPolicy>::DataAssociation(System<tModel>& sys, DataAssociationInfo& info) {

    bool associated_with_track = false;;
    for (auto meas_iter = sys.new_meas_.begin(); meas_iter != sys.new_meas_.end(); ++meas_iter) {
        associated_with_track = false;
        for (auto track_iter = sys.models_.begin(); track_iter != sys.models_.end(); ++track_iter) {
            if(InValidationRegion(sys,*meas_iter,*track_iter)) {
                track_iter->AddNewMeasurement(*meas_iter);
                associated_with_track = true;
            }
        }

        if(!associated_with_track) {
            sys.data_tree_.AddMeasurement(sys,*meas_iter);
        }

    }

    sys.new_meas_.clear();

}

//---------------------------------------------------------------------------------------------------------------
template<typename tModel, template<class> typename tValidationRegionPolicy, template<class> typename tTrackLikelihoodUpdatePolicy, template<class> typename tMeasurementWeightPolicy>
void DataAssociationHost<tModel,tValidationRegionPolicy,tTrackLikelihoodUpdatePolicy,tMeasurementWeightPolicy>::SourceProducedMeasurement(System<tModel>& sys, DataAssociationInfo& info) {

    // Reset the vector
    if(info.source_produced_measurements_.size() != sys.source_container_.num_sources_) {
        info.source_produced_measurements_.clear();
        info.source_produced_measurements_.resize(sys.source_container_.num_sources_,false);
    } else {
        std::fill(info.source_produced_measurements_.begin(), info.source_produced_measurements_.end(), false);
    }



    // If there is only one source and at least one measurement, 
    // then the source produced the measurement
    if (sys.source_container_.num_sources_ ==1 && sys.new_meas_.size() >0) {
         info.source_produced_measurements_[0] = true;
    } else {
        for(auto& meas : sys.new_meas_) {
            info.source_produced_measurements_[0] = true;
        }
    }

}

//---------------------------------------------------------------------------------------------------------------

template<typename tModel, template<class> typename tValidationRegionPolicy, template<class> typename tTrackLikelihoodUpdatePolicy, template<class> typename tMeasurementWeightPolicy>
void DataAssociationHost<tModel,tValidationRegionPolicy,tTrackLikelihoodUpdatePolicy,tMeasurementWeightPolicy>::ResetModelLikelihoodAndNewMeasurements(System<tModel>& sys) {
    
    for (auto track_iter = sys.models_.begin(); track_iter != sys.models_.end(); ++track_iter) {

        for (int ii =0; ii < track_iter->new_assoc_meas_.size(); ++ii) {
            track_iter->model_likelihood_update_info_[ii].Reset();
            track_iter->new_assoc_meas_[ii].clear();
        }
    }
}


} // namespace rransac



#endif // RRANSAC_COMMON_DATA_ASSOCIATION_DATA_ASSOCIATION_HOST_H_