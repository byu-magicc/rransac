#include "common/sources/source_base.h"
#include "common/measurement/measurement_base.h"

namespace rransac
{
    
/**
 * Implements the spatial distance for a measurement of type MeasurementTypes::R2_POSE
 */ 
template<>
float SourceBase::GetSpatialDistance<MeasurementTypes::R2_POSE,MeasurementTypes::R2_POSE>(const Meas& meas1, const Meas meas2, const Parameters& params) {
    return (meas1.data-meas2.data).norm();
}

/**
 * Implements the spatial distance for a measurement of type MeasurementTypes::R2_POSE_TWIST
 */ 
template<>
float SourceBase::GetSpatialDistance<MeasurementTypes::R2_POSE_TWIST,MeasurementTypes::R2_POSE_TWIST>(const Meas& meas1, const Meas meas2, const Parameters& params) {
        return (meas1.data.block(0,0,2,1)-meas2.data.block(0,0,2,1)).norm();
}

// /**
//  * Implements the spatial distance for a measurement of type MeasurementTypes::R2_POSE and MeasurementTypes::R2_POSE_TWIST
//  */ 
// template<>
// float SourceBase::GetSpatialDistance(const MeasR2PosVel& meas1, const MeasR2Pos meas2, const Parameters& params) {
//     return (meas1.data.block(0,0,2,1)-meas2.data).norm();
// }

// /**
//  * Implements the spatial distance for a measurement of type MeasR2PosVel and MeasR2Pos
//  */ 
// template<>
// float SourceBase::GetSpatialDistance(const MeasR2Pos& meas1, const MeasR2PosVel meas2, const Parameters& params) {
//     return (meas1.data-meas2.data.block(0,0,2,1)).norm();
// }


} // namespace rransac
