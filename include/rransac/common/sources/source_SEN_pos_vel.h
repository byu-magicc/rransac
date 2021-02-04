#ifndef RRANSAC_COMMON_SOURCES_SOURCE_SEN_POS_VEL_H_
#define RRANSAC_COMMON_SOURCES_SOURCE_SEN_POS_VEL_H_
#pragma once


#include <typeinfo>

#include "common/sources/source_base.h"
#include "utilities.h"

namespace rransac {


/**
 * \class SourceSENPosVel
 * This source is meant to be used for a target that evolves on SEN and whose
 * measurements space is RN. It is assumed that the Jacobian of the observation
 * function w.r.t. the state is identity. Compatible with 
 * MeasurementType::SEN_POS and MeasurementType::SEN_POS_VEL 
 */ 
template<class tState>
class SourceSENPosVel: public SourceBase<tState,SourceSENPosVel<tState>> {

public:


typedef tState State;                                                            /**< The state of the target. @see State. */
typedef typename tState::DataType DataType;                                      /**< The scalar object for the data. Ex. float, double, etc. */
typedef Eigen::Matrix<DataType,Eigen::Dynamic,Eigen::Dynamic> MatXd;             /**< The object type of the Jacobians. */
static constexpr unsigned int meas_space_dim_ = tState::Group::dim_pos_;         /**< The dimension of the measurement space. */
static constexpr unsigned int meas_pose_rows_ = tState::Group::dim_pos_;         /**< The number of rows in the pose measurement. */
static constexpr unsigned int meas_pose_cols_ = 1;                               /**< The number of columns in the pose measurement. */
static constexpr unsigned int meas_twist_rows_ = tState::Group::dim_pos_;        /**< The number of rows in the twist measurement. */
static constexpr unsigned int meas_twist_cols_ = 1;                              /**< The number of columns in the twist measurement. */
typedef utilities::CompatibleWithModelSENPosVel ModelCompatibility;              /**< Indicates which model this source is compatible with. */


static constexpr unsigned int l_dim_ =  tState::Algebra::dim_a_vel_ + 1;         /**< The dimension of the angular velocity of the target plus one. */
static constexpr unsigned int cov_dim_ = tState::Group::dim_ + tState::Algebra::dim_ - tState::Algebra::dim_t_vel_ + 1; /**< The dimension of the state covariance. */


static_assert(lie_groups::utilities::StateIsSEN_seN<tState>::value, "The state is not compatible with the model");



/** 
 * Initializes the measurement source by initializing the Jacobians. 
 * @param[in] params The source parameters.
 */
void DerivedInit(const SourceParameters& params);      

/** 
 * Returns the jacobian of the observation function w.r.t. the states. This is an optimized version. 
 * @param[in] state The state of a target at which the Jacobian is be evaluated.
 */
MatXd DerivedGetLinObsMatState(tState const& state) const;  

/** 
 * Returns the jacobian of the observation function w.r.t. the states. This method is not optimized. 
 * @param[in] state The state of a target at which the Jacobian is be evaluated.
 */
static MatXd DerivedGetLinObsMatState(const State& state, const MeasurementTypes type);

/** 
 * Returns the jacobian of the observation function w.r.t. the sensor noise. This is an optimized version. 
 * @param[in] state The state of a target at which the Jacobian is be evaluated.
 */
MatXd DerivedGetLinObsMatSensorNoise(const tState& state) const {return this->V_;}   

/** 
 * Returns the jacobian of the observation function w.r.t. the sensor noise. This method is not optimized
 * @param[in] state The state of a target at which the Jacobian is be evaluated.
 */
static MatXd DerivedGetLinObsMatSensorNoise(const State& state, const MeasurementTypes type);

/** 
 * This is an optimized function that implements the observation function and returns an estimated measurement based on the state.
 * @param[in] state A state of the target.
 */
Meas<DataType> DerivedGetEstMeas(const tState& state) const ;

/**
 *  Implements the observation function and returns an estimated measurement based on the state. 
 * @param[in] state A state of the target.
 */
static Meas<DataType> DerivedGetEstMeas(const State& state, const MeasurementTypes type);

/**
 * Performs the OMinus operation between two measurement (m1 ominus m2) of the same type. In other words, this
 * method computes the geodesic distance between two measurements of the same type.
 * @param[in] m1 a measurement
 * @param[in] m2 a measurement
 */
static MatXd DerivedOMinus(const Meas<DataType>& m1, const Meas<DataType>& m2);

/**
 * Generates a random measurement from a Gaussian distribution with mean defined by the state and covariance defined by meas_cov
 * @param[in] state The state that serves as the mean in the Gaussian distribution
 * @param[in] meas_std The measurement standard deviation
 */ 
Meas<DataType> DerivedGenerateRandomMeasurement(const tState& state, const MatXd& meas_std);


};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                            Definitions
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template<class tState>
void SourceSENPosVel<tState>::DerivedInit(const SourceParameters& params) {

    // Verify measurement type
    switch (this->params_.type_)
    {
    case MeasurementTypes::SEN_POS:
        this->V_ = Eigen::Matrix<DataType,meas_space_dim_,meas_space_dim_>::Identity();
        this->H_ = Eigen::Matrix<DataType, meas_space_dim_, this->cov_dim_>::Zero();
        break;
    case MeasurementTypes::SEN_POS_VEL:
        this->V_ = Eigen::Matrix<DataType,meas_space_dim_ + meas_space_dim_,meas_space_dim_+ tState::Algebra::dim_t_vel_>::Identity();
        this->H_ = Eigen::Matrix<DataType,meas_space_dim_+tState::Algebra::dim_t_vel_, this->cov_dim_>::Zero();
        break;
    default:
        throw std::runtime_error("SourceSENPosVel::Init Measurement type not supported.");
        break;
    }
}

//-----------------------------------------------------------------
template<class tState>
Eigen::Matrix<typename tState::DataType,Eigen::Dynamic,Eigen::Dynamic> SourceSENPosVel<tState>::DerivedGetLinObsMatState(const tState& state) const {

    MatXd H = this->H_;

    switch (this->params_.type_)
    {
    case MeasurementTypes::SEN_POS:
        
        
        H.block(0,0,meas_space_dim_, meas_space_dim_) = state.g_.R_;
        return H;
        break;
    case MeasurementTypes::SEN_POS_VEL:

        H.block(0,0,meas_space_dim_, meas_space_dim_) = state.g_.R_; //dt/dt
        H.block(meas_space_dim_, tState::Group::dim_, meas_space_dim_,1) = state.g_.R_.block(0,0,tState::Algebra::dim_t_vel_,1); //dtd/rho_x
        if ( meas_space_dim_== 2) { 
            H.block(meas_space_dim_,meas_space_dim_, meas_space_dim_, tState::Group::RotAlgebra::dim_) = state.g_.R_ * tState::Group::RotAlgebra::Wedge(Eigen::Matrix<DataType,tState::Group::RotAlgebra::dim_,1>::Ones()) * state.u_.p_;
        } else {
            H.block(meas_space_dim_,meas_space_dim_, meas_space_dim_, tState::Group::RotAlgebra::dim_) = -state.g_.R_ * tState::Group::RotAlgebra::Wedge(state.u_.p_.block(0,0,tState::Group::RotAlgebra::dim_,1));
        }
        return H;
        break;
    default:
        throw std::runtime_error("SourceSENPosVel::GetLinObsMatState Measurement type not supported.");
        break;
    }

}

//-----------------------------------------------------------------
template<class tState>
Meas<typename tState::DataType> SourceSENPosVel<tState>::DerivedGetEstMeas(const tState& state) const{
    Meas<DataType> m;
    m.pose = state.g_.t_;
    m.type = this->params_.type_;

    if(this->params_.type_ == MeasurementTypes::SEN_POS_VEL)
        m.twist = state.g_.R_*state.u_.p_;

    return m;
} 

//-----------------------------------------------------------------
template<class tState>
Eigen::Matrix<typename tState::DataType,Eigen::Dynamic,Eigen::Dynamic> SourceSENPosVel<tState>::DerivedOMinus(const Meas<DataType>& m1, const Meas<DataType>& m2) {

    if (m1.type == MeasurementTypes::SEN_POS && m1.type == m2.type) {
        return m1.pose - m2.pose;
    } else if (m1.type == MeasurementTypes::SEN_POS_VEL && m1.type == m2.type){
        Eigen::Matrix<DataType, meas_space_dim_*2,1> error;
        error.block(0,0,meas_space_dim_,1) = m1.pose - m2.pose;
        error.block(meas_space_dim_,0,meas_space_dim_,1) = m1.twist - m2.twist;
        return error;
    } else {
        throw std::runtime_error("SourceSENPosVel::OMinus Measurement type not supported.");
    }


}

//----------------------------------------------------------------------------------------
template<class tState>
Meas<typename tState::DataType> SourceSENPosVel<tState>::DerivedGenerateRandomMeasurement(const tState& state, const MatXd& meas_std){
    Meas<DataType> m;
    m.source_index = this->params_.source_index_;

    MatXd deviation = meas_std*this->GaussianRandomGenerator(meas_std.rows());

    switch (this->params_.type_)
    {
    case MeasurementTypes::SEN_POS:        
        m.pose = state.g_.t_ + deviation;
        m.type = MeasurementTypes::SEN_POS;
        break;
    case MeasurementTypes::SEN_POS_VEL:
        m.pose = state.g_.t_ + deviation.block(0,0,meas_space_dim_,1);   
        m.twist = state.g_.R_*state.u_.p_ + deviation.block(meas_space_dim_,0,meas_space_dim_,1);
        m.type = MeasurementTypes::SEN_POS_VEL;
        break;
    default:
        throw std::runtime_error("SourceSENPosVel::GenerateRandomMeasurement Measurement type not supported.");
        break;
    }

    return m;
}

//----------------------------------------------------------------------------------------

template<class tState>
Meas<typename tState::DataType> SourceSENPosVel<tState>::DerivedGetEstMeas(const tState& state, const MeasurementTypes type) {

    Meas<DataType> m;
    m.pose = state.g_.t_;
    m.type = type;

    if(type == MeasurementTypes::SEN_POS_VEL)
        m.twist = state.g_.R_*state.u_.p_;

    return m;


}

//-------------------------------------------------------------------------------------------------------------------------

template<class tState>
Eigen::Matrix<typename tState::DataType,Eigen::Dynamic,Eigen::Dynamic> SourceSENPosVel<tState>::DerivedGetLinObsMatState(const tState& state, const MeasurementTypes type) {

    const unsigned int sizeg = tState::Group::dim_;
    const unsigned int sizeu = tState::Algebra::dim_;

    // Construct the Jacobians
    MatXd H;
    switch (type)
    {
    case MeasurementTypes::SEN_POS:
        H = Eigen::Matrix<DataType, meas_space_dim_, cov_dim_>::Zero();
        H.block(0,0,meas_space_dim_, meas_space_dim_) = state.g_.R_;
        break;
    case MeasurementTypes::SEN_POS_VEL:
        H = Eigen::Matrix<DataType,meas_space_dim_+tState::Algebra::dim_t_vel_, cov_dim_>::Zero();
        H.block(0,0,meas_space_dim_, meas_space_dim_) = state.g_.R_; //dt/dt
        H.block(meas_space_dim_, tState::Group::dim_, meas_space_dim_,1) = state.g_.R_.block(0,0,tState::Algebra::dim_t_vel_,1); //dtd/rho_x
        if ( meas_space_dim_== 2) { 
            H.block(meas_space_dim_,meas_space_dim_, meas_space_dim_, tState::Group::RotAlgebra::dim_) = state.g_.R_ * tState::Group::RotAlgebra::Wedge(Eigen::Matrix<DataType,tState::Group::RotAlgebra::dim_,1>::Ones()) * state.u_.p_;
        } else {
            H.block(meas_space_dim_,meas_space_dim_, meas_space_dim_, tState::Group::RotAlgebra::dim_) = -state.g_.R_ * tState::Group::RotAlgebra::Wedge(state.u_.p_.block(0,0,tState::Group::RotAlgebra::dim_,1));
        }
        break;
    default:
        throw std::runtime_error("SourceSENPosVel::Init Measurement type not supported.");
        break;
    }
    return H;

}

//-------------------------------------------------------------------------------------------------------------------------

template<class tState>
Eigen::Matrix<typename tState::DataType,Eigen::Dynamic,Eigen::Dynamic> SourceSENPosVel<tState>::DerivedGetLinObsMatSensorNoise(const tState& state, const MeasurementTypes type) {

    const unsigned int sizeg = tState::Group::dim_;
    const unsigned int sizeu = tState::Algebra::dim_;

    MatXd V;

    // Construct the Jacobians
    switch (type)
    {
    case MeasurementTypes::SEN_POS:
        V = Eigen::Matrix<DataType,meas_space_dim_,meas_space_dim_>::Identity();
        break;
    case MeasurementTypes::SEN_POS_VEL:
        V = Eigen::Matrix<DataType,meas_space_dim_ + meas_space_dim_,meas_space_dim_+ tState::Algebra::dim_t_vel_>::Identity();
        break;
    default:
        throw std::runtime_error("SourceRN::Init Measurement type not supported.");
        break;
    }



    return V;

}

// Common Sources
typedef SourceSENPosVel<lie_groups::SE2_se2> SourceSE2PosVel;
typedef SourceSENPosVel<lie_groups::SE3_se3> SourceSE3PosVel;

} // namespace rransac

#endif // RRANSAC_COMMON_SOURCES_SOURCE_SEN_POS_VEL_H_
