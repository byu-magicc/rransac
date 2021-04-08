#ifndef RRANSAC__COMMON__DATA_ASSOCIATION_TRACK_LIKELIHOOD_INFO_POLICIES_TRACK_LIKELIHOOD_INFO_INTEGRATED_PROBABILISTIC_DATA_ASSOCIATION_FILTER_POLICY
#define RRANSAC__COMMON__DATA_ASSOCIATION_TRACK_LIKELIHOOD_INFO_POLICIES_TRACK_LIKELIHOOD_INFO_INTEGRATED_PROBABILISTIC_DATA_ASSOCIATION_FILTER_POLICY

#include <math.h>
#include <numeric>
#include <vector>
#include "rransac/system.h"
#include "rransac/common/models/model_base.h"
#include "rransac/common/data_association/data_association_host.h"

namespace rransac
{

/** \class TLI_IPDAFPolicy
 * This class updates the track likelihood according to the integrated probabilistic data association filter in 
 * Musicki, D., Evans, R., & Stankovic, S. (1994). Integrated probabilistic data association. IEEE Transactions on Automatic Control, 39(6), 1237–1241.
 * The information used to calculate the track likelihood is also used to calculate the weights. Thus the value of delta in equation 2.14 of the paper is 
 * saved in the update track likelihood info. 
 * 
 * This policy should only modify the following member variables of a track: ModelBase::model_likelihood_ and ModelBase::model_likelihood_update_info_. Every time measurements
 * are associated using DataAssociationHost, ModelBase::model_likelihood_update_info_ is reset. 
 * 
 * The Update Track Likelihood Policies must implement two member functions PolicyUpdateTrackLikelihood and PolicyUpdateTrackLikelihoodSingle. The former is responsible 
 * for updating the track likelihood of all the tracks in system and is used by DataAssociationHost. The latter is responsible for updating the likelihood of a single track and is used in Ransac.
 * 
 * This policy should also assume that if there is a measurement from a source, then the track is inside the source's surveillance region.
 * 
 */ 
template<typename tModel>
class TLI_IPDAFPolicy {


public:

static constexpr double decay_rate_ = 0.02;

/**
* Update the track likelihood for all of the tracks. 
* @param[in] sys The object that contains all of the data of RRANSAC. This includes the new measurements, tracks, clusters and data tree. 
* @param[in] info The data association info. It contains which source produced measurements at this sensor scan.
*/ 
static void PolicyUpdateTrackLikelihood(System<tModel>& sys, DataAssociationInfo& info );

/**
 * Updates the track likelihood of a track. 
 * @param[in] sys The object that contains all of the data of RRANSAC. This includes the new measurements, tracks, clusters and data tree. 
 * @param[in,out] track The track's whose likelihood is being updated
 * @param[in] source_index 
 * @param[in] info The data association info. It contains which source produced measurements at this sensor scan.
 * @param[in] dt The time interval from the previous update step to the current update step
 */ 
static void PolicyUpdateTrackLikelihoodSingle(const System<tModel>& sys, tModel& track, DataAssociationInfo& info, const double dt );


private: 

/**
 * Calculates the probability of a measurement given the track.
 * @param[in] sys The object that contains all of the data of RRANSAC. This includes the new measurements, tracks, clusters and data tree. 
 * @param[in] meas The measurement
 * @param[in] track The track   
 */     
static double CalculateMeasurementLikelihood(const System<tModel>& sys, tModel& track, const Meas<typename tModel::DataType>& meas);



};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                            Definitions
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename tModel>
void TLI_IPDAFPolicy<tModel>::PolicyUpdateTrackLikelihood(System<tModel>& sys, DataAssociationInfo& info )  {

    for (auto& track : sys.models_) {

        // Decrease the track likelihood according to the time interval between sensor scans.         
        PolicyUpdateTrackLikelihoodSingle(sys,track,info,sys.dt_);    

    }

}

//---------------------------------------------------------------------------------------------------------------------------

template<typename tModel>
void TLI_IPDAFPolicy<tModel>::PolicyUpdateTrackLikelihoodSingle(const System<tModel>& sys, tModel& track, DataAssociationInfo& info, const double dt ) {

    double alpha = decay_rate_ * dt;
    if (alpha > 0.5) {
        alpha = 0.5;
    }
    track.model_likelihood_ *= (1.0 - alpha);

    for (int source_index = 0; source_index < sys.sources_.size(); ++ source_index) {
        const double& PD = sys.sources_[source_index].params_.probability_of_detection_;
        const double& PG = sys.sources_[source_index].params_.gate_probability_;
        const double& lambda = sys.sources_[source_index].params_.spacial_density_of_false_meas_;

        // If a measurement from a source was associated, it is assumed that the track is inside the surveillance region of the source
        // and that the source produced measurements this sensor scan.
        if(track.new_assoc_meas_[source_index].size() > 0) {

            track.model_likelihood_update_info_[source_index].in_lsr_and_produced_meas = true;
            track.model_likelihood_update_info_[source_index].num_assoc_meas = track.new_assoc_meas_[source_index].size();

            double prob_sum = 0;

            for (auto& meas : track.new_assoc_meas_[source_index]) {
                meas.probability = CalculateMeasurementLikelihood(sys,track,meas);
                prob_sum += meas.probability;
            }


            // Store the value of delta to be used in calculating the measurement weights
            track.model_likelihood_update_info_[source_index].delta  = PD*PG - PD*prob_sum/lambda;


#ifdef DEBUG_BUILD
        if (!sys.sources_[source_index].StateInsideSurveillanceRegion(track.state_)) 
            std::cerr << "TLI_IPDAFPolicy::Single A source produced a measurement for a track that is not in its surveillance region. ";
#endif

        // If the source produced measurements, the track was inside the sources surveillance region, and source didn't observe the track
        } else if (sys.sources_[source_index].StateInsideSurveillanceRegion(track.state_) && info.source_produced_measurements_[source_index]) {
            track.model_likelihood_update_info_[source_index].in_lsr_and_produced_meas = true;
            track.model_likelihood_update_info_[source_index].num_assoc_meas = 0;
            track.model_likelihood_update_info_[source_index].delta = PD*PG;
        } else {
            track.model_likelihood_update_info_[source_index].in_lsr_and_produced_meas = false;
            track.model_likelihood_update_info_[source_index].num_assoc_meas = 0;
            track.model_likelihood_update_info_[source_index].delta =0;
        }

        // If the target is in the local surveillance region of the source and the source produced measurements this sensor scan, update the track likelihood. 
        if (track.model_likelihood_update_info_[source_index].in_lsr_and_produced_meas) {

            track.model_likelihood_  = (1.0-track.model_likelihood_update_info_[source_index].delta)*track.model_likelihood_/(1.0 - track.model_likelihood_update_info_[source_index].delta*track.model_likelihood_);

        }
    }
    

}

//--------------------------------------------------------------------------------------------------------------------------------

template<typename tModel>
double TLI_IPDAFPolicy<tModel>::CalculateMeasurementLikelihood(const System<tModel>& sys, tModel& track, const Meas<typename tModel::DataType>& meas) {

    Eigen::MatrixXd err = sys.sources_[meas.source_index].OMinus(meas, sys.sources_[meas.source_index].GetEstMeas(track.state_));
    Eigen::MatrixXd S = track.GetInnovationCovariance(sys.sources_,meas.source_index);
    double det_inn_cov_sqrt = sqrt(S.determinant());

    return exp( - (err.transpose()*S.inverse()*err)(0,0)/2.0)   /(pow(2.0*M_PI,S.cols()/2.0)*det_inn_cov_sqrt) ;

}



} // namespace rransac


#endif // RRANSAC__COMMON__DATA_ASSOCIATION_TRACK_LIKELIHOOD_INFO_POLICIES_TRACK_LIKELIHOOD_INFO_INTEGRATED_PROBABILISTIC_DATA_ASSOCIATION_FILTER_POLICY
