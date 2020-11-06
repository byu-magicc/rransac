#include "common/sources/source_base.h"

namespace rransac {

template<>
inline void SourceBase::Init<MeasurementTypes::R2_POSE>(const SourceParameters& params) {
    params_ = params;
    Eigen::Matrix<double,2,4> H;
    Eigen::Matrix<double, 2,2> V;
    H.block(0,0,2,2).setIdentity();
    H.block(0,2,2,2).setZero();
    V.setIdentity();

    H_ = H;
    V_ = V;
}

//-------------------------------------------

template <>
inline Eigen::MatrixXd SourceBase::GetLinObsMatState<MeasurementTypes::R2_POSE>(const lie_groups::R2_r2& state) {
    return H_;
}                             

//-------------------------------------------

template <>
inline Eigen::MatrixXd SourceBase::GetLinObsMatSensorNoise<MeasurementTypes::R2_POSE>(const lie_groups::R2_r2& state) {
    return V_;
}

//-------------------------------------------

template <>
inline Eigen::MatrixXd SourceBase::GetEstMeas<MeasurementTypes::R2_POSE>(const lie_groups::R2_r2& state) {
    return state.g_.data_;
}

} // namesapce rransac