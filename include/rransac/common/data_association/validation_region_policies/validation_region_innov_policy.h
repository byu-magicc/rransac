#ifndef RRANSAC__COMMON__DATA_ASSOCIATION__VALIDATION_REGION_POLICIES_VALIDATION_REGION_INNOV_POLICY_H_
#define RRANSAC__COMMON__DATA_ASSOCIATION__VALIDATION_REGION_POLICIES_VALIDATION_REGION_INNOV_POLICY_H_

#include <math.h>
#include <numeric>
#include <vector>
#include "rransac/system.h"

namespace rransac
{

/** \class ValidationiRegionFixedPolicy
 * This class specifies the validation region of a track as an ellipse in the measurement space centered at
 * the estimated measurement of radius SourceParameters::gate_threshold_. Thus, it is a fixed size validation region.
 * 
 */ 
template<typename tModel>
class ValidationRegionInnovPolicy {


public:

typedef typename tModel::Base::Measurement Measurement;

/**
*  Determines if the measurement falls inside the validation region of the track. The validation region is the ellipse in the measurement space of radius SourceParameters::gate_threshold_.
* @param[in] sys The object that contains all of the data of RRANSAC. This includes the new measurements, tracks, clusters and data tree. 
* @param[in] meas The measurement
* @param[in] track The track  
*/ 
static bool PolicyInValidationRegion(const System<tModel>& sys, const Measurement& meas, tModel& track)  {



    Eigen::MatrixXd err = sys.source_container_.OMinus(meas.source_index, meas, sys.source_container_.GetEstMeas(meas.source_index,track.state_,meas.transform_state,meas.transform_data_t_m));
    Eigen::MatrixXd S = track.GetInnovationCovariance(sys.source_container_,meas.source_index,meas.transform_state,meas.transform_data_t_m);

    if( (err.transpose()*S.inverse()*err)(0,0) <= sys.source_container_.GetParams(meas.source_index).gate_threshold_ ) {
        return true;
    } else {
        return false;
    }

}
    

};



} // namespace rransac


#endif // RRANSAC__COMMON__DATA_ASSOCIATION__VALIDATION_REGION_POLICIES_VALIDATION_REGION_INNOV_POLICY_H_
