#include "common/sources/source_base.h"

namespace rransac {

template<>
inline void SourceBase::Init<MeasurementTypes::R2_POSE_TWIST>(const SourceParameters& params) {
    params_ = params;
    Eigen::Matrix<double,4,4> H;
    Eigen::Matrix<double, 4,4> V;
    H.setIdentity();
    V.setIdentity();

    H_ = H;
    V_ = V;
}

//-------------------------------------------

template <>
inline Eigen::MatrixXd SourceBase::GetLinObsMatState<MeasurementTypes::R2_POSE_TWIST>(const lie_groups::R2_r2& state) {
    return H_;
}                             

//-------------------------------------------

template <>
inline Eigen::MatrixXd SourceBase::GetLinObsMatSensorNoise<MeasurementTypes::R2_POSE_TWIST>(const lie_groups::R2_r2& state) {
    return V_;
}

//-------------------------------------------

template <>
inline Eigen::MatrixXd SourceBase::GetEstMeas<MeasurementTypes::R2_POSE_TWIST>(const lie_groups::R2_r2& state) {
    Eigen::Matrix<double,4,1> tmp;
    tmp.block(0,0,2,1) = state.g_.data_;
    tmp.block(2,0,2,1) = state.u_.data_;
    return tmp;
}
} // namesapce rransac