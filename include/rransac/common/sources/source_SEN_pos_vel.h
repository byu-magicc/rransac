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
static constexpr unsigned int l_dim_ =  S::u_type_::dim_a_vel_ + 1;
static constexpr unsigned int dim = S::g_type_::dim_pos_;
static constexpr unsigned int cov_dim_ = S::g_type_::dim_ + S::u_type_::dim_ - S::u_type_::dim_t_vel_ + 1;



/** Initializes the measurement source. This function must set the parameters.  */
void Init(const SourceParameters& params);      

/** Returns the jacobian of the observation function w.r.t. the states 
 * 
*/
Eigen::MatrixXd GetLinObsMatState(S const& state);                        

/** Returns the jacobian of the observation function w.r.t. the sensor noise */
Eigen::MatrixXd GetLinObsMatSensorNoise(const S& state){return this->V_;}                        

/** Computes the estimated measurement given a state */
Meas GetEstMeas(const S& state);

/**
 * Returns the error between the estimated measurement and the measurement
 */
Eigen::MatrixXd OMinus(const Meas& m1, const Meas& m2);

/**
 * Maps the pose to Euclidean space. In this case, it just returns the pose.
 * @param Meas The measurement whose pose needs to be transformed
 */
Eigen::MatrixXd ToEuclidean(const Meas& m)  { return m.pose;}

/**
 * Generates a random measurement from a Gaussian distribution with mean defined by the state and covariance defined by meas_cov
 * @param state The state that serves as the mean in the Gaussian distribution
 * @param meas_std The measurement standard deviation
 */ 
Meas GenerateRandomMeasurement(const S& state, const Eigen::MatrixXd& meas_std);


};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                            Definitions
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template<class S>
void SourceSENPosVel<S>::Init(const SourceParameters& params) {

    // Verify state
    if (typeid(S).name() != typeid(lie_groups::SE2_se2).name() && typeid(S).name() != typeid(lie_groups::SE3_se3).name())
    {
        throw std::runtime_error("SourceSENPosVel::Init State is not supported by this source type");
    }

    // Verify measurement type
    switch (this->params_.type_)
    {
    case MeasurementTypes::SEN_POS:
        this->V_ = Eigen::Matrix<double,S::g_type_::dim_pos_,S::g_type_::dim_pos_>::Identity();
        this->H_ = Eigen::Matrix<double, S::g_type_::dim_pos_, this->cov_dim_>::Zero();
        break;
    case MeasurementTypes::SEN_POS_VEL:
        this->V_ = Eigen::Matrix<double,S::g_type_::dim_pos_ + S::u_type_::dim_t_vel_,S::g_type_::dim_pos_+ S::u_type_::dim_t_vel_>::Identity();
        this->H_ = Eigen::Matrix<double,S::g_type_::dim_pos_+S::u_type_::dim_t_vel_, this->cov_dim_>::Zero();
        // this->H_.block(S::g_type_::dim_pos_,S::g_type_::dim_,S::u_type_::dim_t_vel_,S::u_type_::dim_).setIdentity();
        break;
    default:
        throw std::runtime_error("SourceSENPosVel::Init Measurement type not supported.");
        break;
    }
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

        this->H_.block(0,0,S::g_type_::dim_pos_, S::g_type_::dim_pos_) = state.g_.R_; //dt/dt
        this->H_.block(S::g_type_::dim_pos_, S::g_type_::dim_, S::u_type_::dim_t_vel_,1) = state.g_.R_.block(0,0,S::u_type_::dim_t_vel_,1); //dtd/rho_x
        if ( S::g_type_::dim_pos_== 2) { 
            this->H_.block(S::g_type_::dim_pos_,S::g_type_::dim_pos_, S::u_type_::dim_t_vel_, S::g_type_::rot_algebra::dim_) = state.g_.R_ * S::g_type_::rot_algebra::Wedge(Eigen::Matrix<double,S::g_type_::rot_algebra::dim_,1>::Ones()) * state.u_.p_;
        } else {
            this->H_.block(S::g_type_::dim_pos_,S::g_type_::dim_pos_, S::u_type_::dim_t_vel_, S::g_type_::rot_algebra::dim_) = -state.g_.R_ * S::g_type_::rot_algebra::Wedge(state.u_.p_.block(0,0,S::g_type_::rot_algebra::dim_,1));
        }
        return this->H_;
        break;
    default:
        throw std::runtime_error("SourceSENPosVel::GetLinObsMatState Measurement type not supported.");
        break;
    }

}

//-----------------------------------------------------------------
template<class S>
Meas SourceSENPosVel<S>::GetEstMeas(const S& state) {
    Meas m;
    m.pose = state.g_.t_;

    if(this->params_.type_ == MeasurementTypes::SEN_POS_VEL)
        m.twist = state.g_.R_*state.u_.p_;

    return m;
} 

//-----------------------------------------------------------------
template<class S>
Eigen::MatrixXd SourceSENPosVel<S>::OMinus(const Meas& m1, const Meas& m2) {

    if (this->params_.type_ == MeasurementTypes::SEN_POS) {
        return m1.pose - m2.pose;
    } else if (this->params_.type_ == MeasurementTypes::SEN_POS_VEL){
        Eigen::Matrix<double, S::g_type_::dim_pos_*2,1> error;
        error.block(0,0,S::g_type_::dim_pos_,1) = m1.pose - m2.pose;
        error.block(S::g_type_::dim_pos_,0,S::g_type_::dim_pos_,1) = m1.twist - m2.twist;
        return error;
    } else {
        throw std::runtime_error("SourceSENPosVel::OMinus Measurement type not supported.");
    }


}

//-----------------------------------------------------------------
template<class S>
Meas SourceSENPosVel<S>::GenerateRandomMeasurement(const S& state, const Eigen::MatrixXd& meas_std){
    Meas m;
    m.source_index = this->params_.source_index_;

    Eigen::MatrixXd deviation = meas_std*this->GaussianRandomGenerator(meas_std.rows());

    switch (this->params_.type_)
    {
    case MeasurementTypes::SEN_POS:        
        m.pose = state.g_.t_ + deviation;
        m.type = MeasurementTypes::SEN_POS;
        break;
    case MeasurementTypes::SEN_POS_VEL:
        m.pose = state.g_.t_ + deviation.block(0,0,S::g_type_::dim_pos_,1);   
        m.twist = state.g_.R_*state.u_.p_ + deviation.block(S::g_type_::dim_pos_,0,S::g_type_::dim_pos_,1);
        m.type = MeasurementTypes::SEN_POS_VEL;
        break;
    default:
        throw std::runtime_error("SourceSENPosVel::GenerateRandomMeasurement Measurement type not supported.");
        break;
    }

    return m;
}

} // namespace rransac

#endif // RRANSAC_COMMON_SOURCES_SOURCE_SEN_POS_VEL_H_
