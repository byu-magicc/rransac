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
template<typename _Model>
class ValidationRegionFixedPolicy {


public:

typedef _Model Model;
typedef System<Model> Sys;
typedef typename Model::Base::Measurement Measurement;

/**
*  Determines if the measurement falls inside the validation region of the track. The validation region is the ellipse in the measurement space of radius SourceParameters::gate_threshold_.
* @param[in] sys The object that contains all of the data of RRANSAC. This includes the new measurements, tracks, clusters and data tree. 
* @param[in] meas The measurement
* @param[in] track The track  
*/ 
static bool PolicyInValidationRegion(const Sys& sys, const Measurement& meas, Model& track)  {

    Eigen::MatrixXd err = sys.source_container_.OMinus(meas.source_index, meas, sys.source_container_.GetEstMeas(meas.source_index,track.state_,meas.transform_state,meas.transform_data_t_m));

    if(err.norm() <= sys.source_container_.GetParams(meas.source_index).gate_threshold_ ) {
        return true;
    } else {
        return false;
    }

}
    

};



} // namespace rransac


#endif // RRANSAC__COMMON__DATA_ASSOCIATION__VALIDATION_REGION_POLICIES_VALIDATION_REGION_FIXED_POLICY_H_
