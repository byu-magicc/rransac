#ifndef RRANSAC__COMMON__DATA_ASSOCIATION__VALIDATION_REGION_POLICIES_VALIDATION_REGION_FIXED_POLICY_H_
#define RRANSAC__COMMON__DATA_ASSOCIATION__VALIDATION_REGION_POLICIES_VALIDATION_REGION_FIXED_POLICY_H_

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
class ValidationRegionFixedPolicy {


public:

/**
*  Determines if the measurement falls inside the validation region of the track. The validation region is the ellipse in the measurement space of radius SourceParameters::gate_threshold_.
* @param[in] sys The object that contains all of the data of RRANSAC. This includes the new measurements, tracks, clusters and data tree. 
* @param[in] meas The measurement
* @param[in] track The track  
*/ 
static bool PolicyInValidationRegion(const System<tModel>& sys, const Meas<typename tModel::DataType>& meas, tModel& track)  {

    Eigen::MatrixXd err = sys.sources_[meas.source_index].OMinus(meas, sys.sources_[meas.source_index].GetEstMeas(track.state_));

    if(err.norm() <= sys.sources_[meas.source_index].params_.gate_threshold_ ) {
        return true;
    } else {
        return false;
    }

}
    

};



} // namespace rransac


#endif // RRANSAC__COMMON__DATA_ASSOCIATION__VALIDATION_REGION_POLICIES_VALIDATION_REGION_FIXED_POLICY_H_
