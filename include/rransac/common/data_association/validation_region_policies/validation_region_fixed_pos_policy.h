#ifndef RRANSAC__COMMON__DATA_ASSOCIATION__VALIDATION_REGION_POLICIES_VALIDATION_REGION_FIXED_POS_POLICY_H_
#define RRANSAC__COMMON__DATA_ASSOCIATION__VALIDATION_REGION_POLICIES_VALIDATION_REGION_FIXED_POS_POLICY_H_

#include <math.h>
#include <numeric>
#include <vector>
#include "rransac/system.h"
#include "lie_groups/utilities.h"
#include "rransac/common/measurement/measurement_base.h"
#include "rransac/common/models/model_RN.h"
#include "rransac/common/models/model_SEN_pos_vel.h"
#include "rransac/common/models/model_SEN_pose_twist.h"
#include "rransac/common/sources/source_container.h"
#include "lie_groups/state.h"
#include "rransac/common/utilities.h"

namespace rransac
{

/** \class ValidationiRegionFixedPolicy
 * This class specifies the validation region of a track as an ellipse in the position dimension of the measurement space centered at
 * the position dimension of the estimated measurement of radius SourceParameters::gate_threshold_. Thus, it is a fixed size validation region. This policy is compatible with tracks
 * whose configuration manifold is RN, SE2, or SE3 whose measurement space is position and its derivatives with the measured position being Meas::pose
 */ 


template<typename tModel>
class ValidationRegionFixedPosPolicy {




public:


/**
*  Determines if the measurement falls inside the validation region of the track. The validation region is the ellipse in the measurement space of radius SourceParameters::gate_threshold_.
* @param[in] sys The object that contains all of the data of RRANSAC. This includes the new measurements, tracks, clusters and data tree. 
* @param[in] meas The measurement
* @param[in] track The track  
*/ 
static bool PolicyInValidationRegion(const System<tModel>& sys, const Meas<typename tModel::DataType>& meas, tModel& track)  {

    static_assert(utilities::AlwaysFalse<tModel>::value,"ValidationRegionFixedPosPolicy::InValidationRegion This function is not implemented for the specified model type.");

    return false;
}   

};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Class Specializations
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename tDataType, int tN, template <typename > typename tTransformation, MeasurementTypes tMeasurementType, template <typename, MeasurementTypes, template <typename > typename> typename tS0, typename tS1, typename tS2, typename tS3, typename tS4>
class ValidationRegionFixedPosPolicy<ModelRN<SourceContainer<tS0<lie_groups::State<lie_groups::Rn,tDataType,tN>,tMeasurementType,tTransformation>,tS1,tS2,tS3,tS4>>>{

    typedef ModelRN<SourceContainer<tS0<lie_groups::State<lie_groups::Rn,tDataType,tN>,tMeasurementType,tTransformation>,tS1,tS2,tS3,tS4>> Model;
public:
    static bool PolicyInValidationRegion(const System<Model>& sys, const Meas<typename Model::DataType>& meas, Model& track)  {

        Meas<tDataType> m = sys.source_container_.GetEstMeas(meas.source_index,track.state_,meas.transform_state,meas.transform_data_t_m);

        Eigen::Matrix<double,tN,1> err = m.pose - meas.pose;

        if(err.norm() <= sys.source_container_.GetParams(meas.source_index).gate_threshold_) {
            return true;
        } else {
            return false;
        }

}

};

//-------------------------------------------------------------------------------------------------------------------------------------

template<typename tDataType, int tN, template <typename > typename tTransformation, MeasurementTypes tMeasurementType, template <typename, MeasurementTypes, template <typename > typename> typename tS0, typename tS1, typename tS2, typename tS3, typename tS4>
class ValidationRegionFixedPosPolicy<ModelSENPosVel<SourceContainer<tS0<lie_groups::State<lie_groups::SE2,tDataType,tN>,tMeasurementType,tTransformation>,tS1,tS2,tS3,tS4>>>{

    typedef ModelSENPosVel<SourceContainer<tS0<lie_groups::State<lie_groups::SE2,tDataType,tN>,tMeasurementType,tTransformation>,tS1,tS2,tS3,tS4>> Model;
public:
    static bool PolicyInValidationRegion(const System<Model>& sys, const Meas<typename Model::DataType>& meas, Model& track)  {

#ifdef DEBUG_BUILD
       if (meas.type != SEN_POS || meas.type != SEN_POS_VEL)
        throw std::runtime_error("ValidationRegionFixedPosPolicy::PolicyInValidationRegion Measurement type is invalid. It must be SEN_POS or SEN_POS_VEL")
#endif

        Meas<tDataType> m = sys.source_container_.GetEstMeas(meas.source_index,track.state_,meas.transform_state,meas.transform_data_t_m);
        Eigen::MatrixXd err = m.pose - meas.pose;

        if(err.norm() <= sys.source_container_.GetParams(meas.source_index).gate_threshold_) {
            return true;
        } else {
            return false;
        }

}

};

//-------------------------------------------------------------------------------------------------------------------------------------

template<typename tDataType, int tN, template <typename > typename tTransformation, MeasurementTypes tMeasurementType, template <typename, MeasurementTypes, template <typename > typename> typename tS0, typename tS1, typename tS2, typename tS3, typename tS4>
class ValidationRegionFixedPosPolicy<ModelSENPosVel<SourceContainer<tS0<lie_groups::State<lie_groups::SE3,tDataType,tN>,tMeasurementType,tTransformation>,tS1,tS2,tS3,tS4>>>{

     typedef ModelSENPosVel<SourceContainer<tS0<lie_groups::State<lie_groups::SE3,tDataType,tN>,tMeasurementType,tTransformation>,tS1,tS2,tS3,tS4>> Model;
public:
    static bool PolicyInValidationRegion(const System<Model>& sys, const Meas<typename Model::DataType>& meas, Model& track)  {

#ifdef DEBUG_BUILD
       if (meas.type != SEN_POS || meas.type != SEN_POS_VEL)
        throw std::runtime_error("ValidationRegionFixedPosPolicy::PolicyInValidationRegion Measurement type is invalid. It must be SEN_POS or SEN_POS_VEL")
#endif

        Meas<tDataType> m = sys.source_container_.GetEstMeas(meas.source_index,track.state_,meas.transform_state,meas.transform_data_t_m);
        Eigen::MatrixXd err = m.pose - meas.pose;

        if(err.norm() <= sys.source_container_.GetParams(meas.source_index).gate_threshold_) {
            return true;
        } else {
            return false;
        }

}

};



} // namespace rransac


#endif // RRANSAC__COMMON__DATA_ASSOCIATION__VALIDATION_REGION_POLICIES_VALIDATION_REGION_FIXED_POS_POLICY_H_
