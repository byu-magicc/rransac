#ifndef RRANSAC__COMMON__DATA_ASSOCIATION__MEASUREMENT_WEIGHT_POLICY_MW_IPDAF_POLICY
#define RRANSAC__COMMON__DATA_ASSOCIATION__MEASUREMENT_WEIGHT_POLICY_MW_IPDAF_POLICY


#include "rransac/common/data_association/track_likelihood_info_policies/tli_ipdaf_policy.h"
#include "rransac/common/data_association/data_association_host.h"
#include "rransac/system.h"

namespace rransac
{

/** \class MW_IPDAFPolicy 
 * This policy class is inherited by DataAssociationHost when selected. It is responsible to calculate the weights for each track associated measurement. 
 * This class assumes that the track likelihood info policy used is TLI_IPDAFPolicy; otherwise, the behavior is undefined. The weights are calculated according 
 * to the integrated probabilistic data association filter in 
 * Musicki, D., Evans, R., & Stankovic, S. (1994). Integrated probabilistic data association. IEEE Transactions on Automatic Control, 39(6), 1237–1241.
 * The value of delta in equation 2.14 and the measurement probabilities are calculated by the policy TLI_IPDAFPolicy and stored in ModelBase::model_likelihood_update_info_. 
 * This value of delta and the measurement probabilities are used to calculate the weights according to equation 2.18 in the referenced paper. 
 * 
 * The Measurement Weight Policies must implement two member functions PolicyCalculateMeasurementWeight and PolicyCalculateMeasurementWeightSingle. The former is responsible 
 * for calculating the measurement weights for all of the track associated measurements and is used by the DataAssociationHost. The latter is responsible for calculating
 * the measurement weights for all of the measurements associated to a single track. 
 * 
 * Lastly, this policy must have a struct titled CompatiblityCheck. It is given the validation region policy and the update track likelihood policy as template parameters. This struct 
 * calculates if these other policies are compatible with this policy. 
 */ 

template<typename tModel>
class MW_IPDAFPolicy {

public:

    /** 
     * Calculates the weights for each track associated measurement. 
     * @param[in,out] sys The object that contains all of the data of RRANSAC. This includes the new measurements, tracks, clusters and data tree. 
     */
    static void PolicyCalculateMeasurementWeight(System<tModel>& sys, DataAssociationInfo& info) {

        for (auto& track: sys.models_) {
            PolicyCalculateMeasurementWeightSingle(sys,track,info);
        }

    }

    /**
     * Calculates the weights for each measurement associated to the track.
     * @param[in] info DataAssociationInfo
     * @param[in] track the track whose new measurement weights are being calculated
     * @param[in,out] sys The object that contains all of the data of RRANSAC. This includes the new measurements, tracks, clusters and data tree. 
     */ 
    static void PolicyCalculateMeasurementWeightSingle(const System<tModel>& sys, tModel& track, DataAssociationInfo& info) {

        for (int source_index =0; source_index < sys.source_container_.num_sources_; ++source_index) {
            for(auto& meas : track.new_assoc_meas_[source_index]) {
                meas.weight = sys.source_container_.GetParams(source_index).probability_of_detection_*meas.probability/sys.source_container_.GetParams(source_index).spacial_density_of_false_meas_/(1.0 - track.model_likelihood_update_info_[source_index].delta);
            }
        }

    }


    template<template<class> typename tValidationRegionPolicy, template<class> typename tUpdateTrackLikelihoodPolicy>
    struct CompatiblityCheck {
        static constexpr bool value = std::is_same<tUpdateTrackLikelihoodPolicy<tModel>,TLI_IPDAFPolicy<tModel>>::value;
        
    }; 



};
    
} // namespace rransac


#endif //RRANSAC__COMMON__DATA_ASSOCIATION__MEASUREMENT_WEIGHT_POLICY_MW_IPDAF_POLICY