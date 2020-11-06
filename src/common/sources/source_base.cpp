#include "common/sources/source_base.h"

namespace rransac
{
template< class S>
SourceBase<S>::SourceBase() {

    gsd_ptr_ = new GSDFuncPTR *[MeasurementTypes::NUM_TYPES];
    for(int i = 0; i < MeasurementTypes::NUM_TYPES; ++i)
    {
        gsd_ptr_[i] = new GSDFuncPTR[MeasurementTypes::NUM_TYPES];
    }

    gsd_ptr_[MeasurementTypes::R2_POSE][MeasurementTypes::R2_POSE]             = &GSDR2PoseR2Pose;
    gsd_ptr_[MeasurementTypes::R2_POSE][MeasurementTypes::R2_POSE_TWIST]       = &GSDR2PoseR2PoseTwist;
    gsd_ptr_[MeasurementTypes::R2_POSE_TWIST][MeasurementTypes::R2_POSE]       = &GSDR2PoseTwistR2Pose;
    gsd_ptr_[MeasurementTypes::R2_POSE_TWIST][MeasurementTypes::R2_POSE_TWIST] = &GSDR2PoseTwistR2PoseTwist;



}

//---------------------------------------------------

template< class S>
SourceBase<S>::~SourceBase() {

    for (int i = 0; i < MeasurementTypes::NUM_TYPES; i++) {
        delete [] gsd_ptr_[i];
    }
    delete [] gsd_ptr_;

}

///////////////////////////////// Private ////////////////////////////////////////


} // namespace rransac



