#ifndef RRANSAC_COMMON_DATA_ASSOCIATION_DATA_ASSOCIATION_HOST_H_
#define RRANSAC_COMMON_DATA_ASSOCIATION_DATA_ASSOCIATION_HOST_H_
#pragma once

#include <vector>

#include "rransac/system.h"

namespace rransac
{

template<typename _TransformDataType>
struct DataAssociationInfo {

    std::vector<bool> source_produced_measurements_;                                         /**< When a sensor scan occurs one or more sources can produce measurements. This flag indicates if a source 
                                                                                                  during the latest sensor scan produced measurements. */
    std::vector<bool> transform_state_;                                                      /**< Indicates if a measurement from a source requires the state to be transformed into the measurement frame. */
    std::vector<_TransformDataType> transform_data_t_m_;                                        /**< The data required to transform the state into the measurement frame. */

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
template<typename _Model, template<class> typename _ValidationRegionPolicy, template<class> typename _UpdateTrackLikelihoodPolicy, template<class> typename _MeasurementWeightPolicy>
class DataAssociationHost : public _ValidationRegionPolicy<_Model>, _UpdateTrackLikelihoodPolicy<_Model>, _MeasurementWeightPolicy<_Model>{

// Some policies are only compatible with others. These assertions enforce compatability. 
static_assert(_MeasurementWeightPolicy<_Model>::template CompatiblityCheck<_ValidationRegionPolicy,_UpdateTrackLikelihoodPolicy>::value,"MW_IPDAFPolicy::CompatiblityCheck The policy MW_IPDAFPolicy is only compatible with the update track likelihood policy TLI_IPDAFPolicy." );



public:

    typedef _Model Model;
    typedef typename Model::Base::Measurement Measurement;
    typedef typename Model::Base::TransformDataType TransformDataType;
    typedef DataAssociationInfo<TransformDataType> DataAssociationInfoT;
    typedef System<Model> Sys;

    /**
     * Associates new measurements to the tracks, clusters and the data tree. 
     * @param[in,out] sys The object that contains all of the data of RRANSAC. This includes the new measurements, tracks, clusters and data tree. 
     */ 
    void AssociateNewMeasurements(Sys& sys) {
        CalculateDataAssociationInfo(sys, this->data_association_info_);                    // Determine which sources produced measurements
        ResetModelLikelihoodAndNewMeasurements(sys);  
        DataAssociation(sys,this->data_association_info_);
        UpdateTrackLikelihood(sys,this->data_association_info_);
        CalculateMeasurementWeight(sys,this->data_association_info_);
        

        if(sys.new_meas_.size() != 0)   
            throw std::runtime_error("DataAssociationHost::AssociateNewMeasurements: All new measurements should have been copied to a model, cluster, or data tree and removed from System<Model>::new_meas_");
    }

    const DataAssociationInfoT& GetDataAssociationInfo() const {return data_association_info_;}

private:

    /**
     * Associates measurements to at least one track or the data tree. If a measurement falls within the validation region of a track, then it is associated with the track; otherwise it is added
     * to the data tree. The validation region is determined by the policy. 
     * @param[in,out] sys The object that contains all of the data of RRANSAC. This includes the new measurements, tracks, clusters and data tree. 
     */ 
    void DataAssociation(Sys& sys, DataAssociationInfoT& info);

    /**
     *  Determines if the measurement falls inside the validation region of the track according to the policy. 
     * @param[in] sys The object that contains all of the data of RRANSAC. This includes the new measurements, tracks, clusters and data tree. 
     * @param[in] meas The measurement
     * @param[in] track The track. It is not a constant reference so be careful!! 
     */ 
    bool InValidationRegion(const Sys& sys, const Measurement& meas, Model& track) {
        return DataAssociationHost::PolicyInValidationRegion(sys, meas,track);
    }

    /** 
     * Uses the new associated measurements to update the track's likelihood. 
     * @param[in,out] sys The object that contains all of the data of RRANSAC. This includes the new measurements, tracks, clusters and data tree. 
     */
    void UpdateTrackLikelihood(Sys& sys, DataAssociationInfoT& info ) {
        DataAssociationHost::PolicyUpdateTrackLikelihood(sys,info);
    }

    /** 
     * Calculates the weights for each track associated measurement. 
     * @param[in,out] sys The object that contains all of the data of RRANSAC. This includes the new measurements, tracks, clusters and data tree. 
     */
    void CalculateMeasurementWeight(Sys& sys, DataAssociationInfoT& info) {
        DataAssociationHost::PolicyCalculateMeasurementWeight(sys,info);
    }

    /**
     * Calculates the data association info.
     * Sets the member vairable DataAssociationInfo::source_produced_measurements_ to true if the source produced a measurement this sensor scan; otherwise false.
     * Sets the member vairable DataAssociationInfo::transform_state_ to true if the track needs to be transformed into the measurement frame. This info is taken from the measurement.
     * Sets the member vairable DataAssociationInfo::transform_data_t_m_ to contain the transformation data necessary to transform the track to the measurement frame. This info is taken from the measurement.
     * @param[in,out] sys The object that contains all of the data of RRANSAC. This includes the new measurements, tracks, clusters and data tree. 
     */
    void CalculateDataAssociationInfo(Sys& sys, DataAssociationInfoT& info);

    /** 
     * Ensures that none of the tracks have any new associated measurements in ModelBase::new_assoc_meas_, and it resets the
     * model likelihood update info.
     */ 
    void ResetModelLikelihoodAndNewMeasurements(Sys& sys);

    DataAssociationInfoT data_association_info_;

};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                            Definitions
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//------------------------------------------------------------------------------------------------------------------------------------

template<typename _Model, template<class> typename _ValidationRegionPolicy, template<class> typename _UpdateTrackLikelihoodPolicy, template<class> typename _MeasurementWeightPolicy>
void DataAssociationHost<_Model,_ValidationRegionPolicy,_UpdateTrackLikelihoodPolicy,_MeasurementWeightPolicy>::DataAssociation(Sys& sys, DataAssociationInfoT& info) {

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
template<typename _Model, template<class> typename _ValidationRegionPolicy, template<class> typename _UpdateTrackLikelihoodPolicy, template<class> typename _MeasurementWeightPolicy>
void DataAssociationHost<_Model,_ValidationRegionPolicy,_UpdateTrackLikelihoodPolicy,_MeasurementWeightPolicy>::CalculateDataAssociationInfo(Sys& sys, DataAssociationInfoT& info) {

    TransformDataType transform_data;

    // Reset the vector
    if(info.source_produced_measurements_.size() != sys.source_container_.num_sources_) {
        info.source_produced_measurements_.clear();
        info.transform_state_.clear();
        info.source_produced_measurements_.resize(sys.source_container_.num_sources_,false);
        info.transform_state_.resize(sys.source_container_.num_sources_,false);
        info.transform_data_t_m_.resize(sys.source_container_.num_sources_,transform_data);
    } else {
        std::fill(info.source_produced_measurements_.begin(), info.source_produced_measurements_.end(), false);
        std::fill(info.transform_state_.begin(), info.transform_state_.end(), false);
        std::fill(info.transform_data_t_m_.begin(), info.transform_data_t_m_.end(), transform_data);
    }



    // If there is only one source and at least one measurement, 
    // then the source produced the measurement
    for(auto& meas :sys.new_meas_) {
        if (info.source_produced_measurements_[meas.source_index] == false) {
            info.source_produced_measurements_[meas.source_index] = true;
            info.transform_state_[meas.source_index] = meas.transform_state;
            info.transform_data_t_m_[meas.source_index] = meas.transform_data_t_m;
        }
        // Stop searching if we have all the information we need.
        if(std::all_of(info.source_produced_measurements_.begin(),info.source_produced_measurements_.end(),[](bool v){return v;})) {
            break;
        }
    }


}

//---------------------------------------------------------------------------------------------------------------

template<typename _Model, template<class> typename _ValidationRegionPolicy, template<class> typename _UpdateTrackLikelihoodPolicy, template<class> typename _MeasurementWeightPolicy>
void DataAssociationHost<_Model,_ValidationRegionPolicy,_UpdateTrackLikelihoodPolicy,_MeasurementWeightPolicy>::ResetModelLikelihoodAndNewMeasurements(Sys& sys) {
    
    for (auto track_iter = sys.models_.begin(); track_iter != sys.models_.end(); ++track_iter) {

        for (int ii =0; ii < track_iter->new_assoc_meas_.size(); ++ii) {
            track_iter->model_likelihood_update_info_[ii].Reset();
            track_iter->new_assoc_meas_[ii].clear();
        }
    }
}


} // namespace rransac



#endif // RRANSAC_COMMON_DATA_ASSOCIATION_DATA_ASSOCIATION_HOST_H_