#ifndef RRANSAC_COMMON_SOURCES_SOURCE_SEN_POSE_TWIST_H_
#define RRANSAC_COMMON_SOURCES_SOURCE_SEN_POSE_TWIST_H_

#include <typeinfo>

#include "common/sources/source_base.h"
namespace rransac {

/**
 * \class SourceSENPoseTwist
 * This souce is meant to be used for a target that evolves on SEN and whose
 * measurements are also on SEN or SEN and its tangent space. Compatible with 
 * MeasurementType::SEN_POSE and MeasurementType::SEN_POSE_TWIST
 */ 

template <class S>
class SourceSENPoseTwist : public SourceBase<S,SourceSENPoseTwist<S>> {

public:

typedef S type_;

/** Initializes the measurement source. This function must set the parameters.  */
void Init(const SourceParameters& params);      

/** Returns the jacobian of the observation function w.r.t. the states 
 * 
*/
Eigen::MatrixXd GetLinObsMatState(S const& state){return this->H_;};                        

/** Returns the jacobian of the observation function w.r.t. the sensor noise */
Eigen::MatrixXd GetLinObsMatSensorNoise(const S& state){return this->V_;}                        

/** Computes the estimated measurement given a state */
Meas GetEstMeas(const S& state) {
    Meas m;
    m.pose = state.g_.data_;
    m.twist = state.u_.data_;
    return m;
    } 

/**
 * Maps the pose to Euclidean space. The translation is unchanged; however, the rotation is transformed using Cayley coordinates of the first kind.
 * @param Meas The measurement whose pose needs to be transformed
 */
Eigen::MatrixXd ToEuclidean(const Meas& m);

};



//-----------------------------------------------------------------

template<class S>
void SourceSENPoseTwist<S>::Init(const SourceParameters& params) {

    // Verify state
    if (typeid(S).name() != typeid(lie_groups::SE2_se2).name() && typeid(S).name() != typeid(lie_groups::SE3_se3).name())
    {
        throw std::runtime_error("SourceSENPoseTwist::Init State is not supported by this source type");
    }

    // Verify measurement type
    switch (params.type_)
    {
    case MeasurementTypes::SEN_POSE:
        this->V_ = Eigen::Matrix<double,S::g_type_::dim_,S::g_type_::dim_>::Identity();
        this->H_ = Eigen::Matrix<double, S::g_type_::dim_, S::dim_>::Zero();
        this->H_.block(0,0,S::g_type_::dim_,S::g_type_::dim_).setIdentity();
        break;
    case MeasurementTypes::SEN_POSE_TWIST:
        this->V_ = Eigen::Matrix<double,S::dim_,S::dim_>::Identity();
        this->H_ = Eigen::Matrix<double,S::dim_,S::dim_>::Identity();
        break;
    default:
        throw std::runtime_error("SourceSENPoseTwist::Init Measurement type not supported.");
        break;
    }

    this->params_ = params;

}


//----------------------------------------------------------------------------------------
template<class S>
Eigen::MatrixXd SourceSENPoseTwist<S>::ToEuclidean(const Meas& m)  {
    typename S::g_type_ g(m.pose);
    Eigen::Matrix<double,  S::g_type_::dim_pos_,S::g_type_::dim_pos_> I = Eigen::Matrix<double,  S::g_type_::dim_pos_,S::g_type_::dim_pos_>::Identity();
    Eigen::Matrix<double, S::g_type_::dim_pos_,S::g_type_::dim_pos_>&& u = 2.0*(g.R_-I)*(g.R_+I).inverse();
    Eigen::Matrix<double, S::g_type_::dim_,1> pose_euclidean;
    pose_euclidean.block(0,0,S::g_type_::dim_pos_,1) = g.t_;
    pose_euclidean.block(S::g_type_::dim_pos_,0,S::g_type_::dim_rot_,1) = S::g_type_::rot_algebra::Vee(u);


    return pose_euclidean;
}





} // namespace rransac
#endif // RRANSAC_COMMON_SOURCES_SOURCE_SEN_POSE_TWIST_H_