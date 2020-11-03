#include "common/sources/source_base.h"

namespace rransac
{

// ///////////////////////////////////////////////////////////////////////
// //                             R2_POS_VEL
// ///////////////////////////////////////////////////////////////////////
// template<>
// inline void SourceBase::Init<SourceTypes::R2_POS>(const SourceParameters& params, unsigned int source_id) {
//     params_ = params;
//     source_id_ = source_id;
//     Eigen::Matrix<double,2,4> H;
//     Eigen::Matrix<double, 3,3> V;
//     H.block(0,0,2,2).setIdentity();
//     H.block(0,2,2,2).setZero();
//     V.setIdentity();

//     H_ = H;
//     V_ = V;
// }

// //-------------------------------------------

// template <>
// inline Eigen::MatrixXd SourceBase::GetLinObsMatState<SourceTypes::R2_POS>(const lie_groups::R2_r2& state) {
//     return H_;
// }                             

// //-------------------------------------------

// template <>
// inline Eigen::MatrixXd SourceBase::GetLinObsMatSensorNoise<SourceTypes::R2_POS,lie_groups::R2_r2>(const lie_groups::R2_r2& state) {
//     return V_;
// }

// //-------------------------------------------

// template <>
// inline Eigen::MatrixXd SourceBase::GetEstMeas<SourceTypes::R2_POS,lie_groups::R2_r2>(const lie_groups::R2_r2& state) {
//     return state.g_.data_;
// }

} // namespace rransac



