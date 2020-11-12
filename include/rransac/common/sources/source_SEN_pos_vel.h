#ifndef RRANSAC_COMMON_SOURCES_SOURCE_SEN_POS_VEL_H_
#define RRANSAC_COMMON_SOURCES_SOURCE_SEN_POS_VEL_H_

#include <typeinfo>

#include "common/sources/source_base.h"

namespace rransac {


/**
 * \class SourceSENPosVel
 * This source is meant to be used for a target that evolves on SEN and whose
 * measurements are on RN or on RN and its tangent space. Compatible with 
 * MeasurementType::SEN_POS and MeasurementType::SEN_POS_VEL 
 */ 
template<class S>
class SourceSENPosVel: public SourceBase<S,SourceSENPosVel<S>> {

public:

typedef S type_;

/** Initializes the measurement source. This function must set the parameters.  */
void Init(const SourceParameters& params);      

/** Returns the jacobian of the observation function w.r.t. the states 
 * 
*/
Eigen::MatrixXd GetLinObsMatState(S const& state);                        

/** Returns the jacobian of the observation function w.r.t. the sensor noise */
Eigen::MatrixXd GetLinObsMatSensorNoise(const S& state){return this->V_;}                        

/** Computes the estimated measurement given a state */
Meas GetEstMeas(const S& state) {
    Meas m;
    m.pose = state.g_.t_;
    m.twist = state.u_.p_;
    return m;
} 

/**
 * Maps the pose to Euclidean space. In this case, it just returns the pose.
 * @param Meas The measurement whose pose needs to be transformed
 */
Eigen::MatrixXd ToEuclidean(const Meas& m)  {
    return m.pose;
}



};

//-----------------------------------------------------------------

template<class S>
void SourceSENPosVel<S>::Init(const SourceParameters& params) {

    // Verify state
    if (typeid(S).name() != typeid(lie_groups::SE2_se2).name() && typeid(S).name() != typeid(lie_groups::SE3_se3).name())
    {
        throw std::runtime_error("SourceSENPosVel::Init State is not supported by this source type");
    }

    // Verify measurement type
    switch (params.type_)
    {
    case MeasurementTypes::SEN_POS:
        this->V_ = Eigen::Matrix<double,S::g_type_::dim_pos_,S::g_type_::dim_pos_>::Identity();
        this->H_ = Eigen::Matrix<double, S::g_type_::dim_pos_, S::dim_>::Zero();
        break;
    case MeasurementTypes::SEN_POS_VEL:
        this->V_ = Eigen::Matrix<double,S::g_type_::dim_pos_ + S::u_type_::dim_t_vel_,S::g_type_::dim_pos_+ S::u_type_::dim_t_vel_>::Identity();
        this->H_ = Eigen::Matrix<double,S::g_type_::dim_pos_+S::u_type_::dim_t_vel_, S::dim_>::Zero();
        this->H_.block(S::g_type_::dim_pos_,S::g_type_::dim_,S::u_type_::dim_t_vel_,S::u_type_::dim_).setIdentity();
        break;
    default:
        throw std::runtime_error("SourceSENPosVel::Init Measurement type not supported.");
        break;
    }

    this->params_ = params;

}

//-----------------------------------------------------------------
template<class S>
Eigen::MatrixXd SourceSENPosVel<S>::GetLinObsMatState(const S& state) {

switch (this->params_.type_)
{
case MeasurementTypes::SEN_POS:
    
    this->H_.block(0,0,S::g_type_::dim_pos_, S::g_type_::dim_pos_) = state.g_.R_;
    return this->H_;
    break;
case MeasurementTypes::SEN_POS_VEL:
    this->H_.block(0,0,S::g_type_::dim_pos_, S::g_type_::dim_pos_) = state.g_.R_;
    return this->H_;
    break;
default:
    throw std::runtime_error("SourceSENPosVel::GetLinObsMatState Measurement type not supported.");
    break;
}


}
} // namespace rransac

#endif // RRANSAC_COMMON_SOURCES_SOURCE_SEN_POS_VEL_H_
