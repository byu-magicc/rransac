#ifndef RRANSAC_COMMON_SOURCES_SOURCE_SEN_POS_VEL_H_
#define RRANSAC_COMMON_SOURCES_SOURCE_SEN_POS_VEL_H_
#pragma once


#include <typeinfo>

#include "rransac/common/sources/source_base.h"
#include "rransac/common/utilities.h"
#include "lie_groups/utilities.h"

namespace rransac {


/**
 * \class SourceSENPosVel
 * This source is meant to be used for a target that evolves on SEN and whose
 * measurements space is RN. It is assumed that the Jacobian of the observation
 * function w.r.t. the state is identity. Compatible with 
 * MeasurementType::SEN_POS and MeasurementType::SEN_POS_VEL 
 */ 
template <typename tState, MeasurementTypes tMeasurementType, template <typename > typename tTransformation>
class SourceSENPosVel: public SourceBase<tState,tMeasurementType, tTransformation, tState::Group::dim_pos_, SourceSENPosVel> {

public:


typedef tState State;                                                            /**< The state of the target. @see State. */
typedef typename tState::DataType DataType;                                      /**< The scalar object for the data. Ex. float, double, etc. */
typedef Eigen::Matrix<DataType,Eigen::Dynamic,Eigen::Dynamic> MatXd;             /**< The object type of the Jacobians. */
static constexpr unsigned int meas_pose_rows_ = tState::Group::dim_pos_;         /**< The number of rows in the pose measurement. */
static constexpr unsigned int meas_pose_cols_ = 1;                               /**< The number of columns in the pose measurement. */
static constexpr unsigned int meas_twist_rows_ = tState::Group::dim_pos_;        /**< The number of rows in the twist measurement. */
static constexpr unsigned int meas_twist_cols_ = 1;                              /**< The number of columns in the twist measurement. */
static constexpr MeasurementTypes measurement_type_ = tMeasurementType;          /**< The measurement type of the source. */
typedef utilities::CompatibleWithModelSENPosVel ModelCompatibility;              /**< Indicates which model this source is compatible with. */
typedef SourceBase<tState,tMeasurementType, tTransformation, tState::Group::dim_pos_, SourceSENPosVel> Base; /**< The Base source class */

static_assert(lie_groups::utilities::StateIsSEN_seN<tState>::value, "SourceSENPosVel: The state is not compatible with the model");
static_assert( tMeasurementType == MeasurementTypes::SEN_POS || tMeasurementType==MeasurementTypes::SEN_POS_VEL, "SourceSENPosVel: The measurement type is not compatible with the source."    );

static constexpr unsigned int l_dim_ =  tState::Algebra::dim_a_vel_ + 1;         /**< The dimension of the angular velocity of the target plus one. */
static constexpr unsigned int cov_dim_ = tState::Group::dim_ + tState::Algebra::dim_ - tState::Algebra::dim_t_vel_ + 1; /**< The dimension of the state covariance. */


/** 
 * Initializes the Jacobians
 */
SourceSENPosVel();

/** 
 * Initializes the measurement source. Currently does nothing. 
 * @param[in] params The source parameters.
 */
void DerivedInit(const SourceParameters& params){}     

/** 
 * Returns the jacobian of the observation function w.r.t. the states.
 * @param[in] state The state of a target at which the Jacobian is be evaluated.
 */
static MatXd DerivedGetLinObsMatState(const State& state);
                   

/** 
 * Returns the jacobian of the observation function w.r.t. the sensor noise.
 * @param[in] state The state of a target at which the Jacobian is be evaluated.
 */
static MatXd DerivedGetLinObsMatSensorNoise(const State& state);


/**
 *  Implements the observation function and returns an estimated measurement based on the state. 
 * Currently, the measurement is only given a pose, twist, and measurement type. 
 * @param[in] state A state of the target.
 */
static Meas<DataType> DerivedGetEstMeas(const State& state);

/**
 * Performs the OMinus operation between two measurement (m1 ominus m2) of the same type. In other words, this
 * method computes the geodesic distance between two measurements of the same type.
 * @param[in] m1 a measurement
 * @param[in] m2 a measurement
 */
static MatXd DerivedOMinus(const Meas<DataType>& m1, const Meas<DataType>& m2);


/**
 * Generates a random measurement from a Gaussian distribution with mean defined by the state and standard deviation of the measurement covariance
 * @param[in] state The state that serves as the mean in the Gaussian distribution
 * @param[in] meas_std The measurement standard deviation
 */ 
Meas<DataType> DerivedGenerateRandomMeasurement(const MatXd& meas_std, const tState& state) const;


};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                            Definitions
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template <typename tState, MeasurementTypes tMeasurementType, template <typename > typename tTransformation>
SourceSENPosVel<tState,tMeasurementType,tTransformation>::SourceSENPosVel() {


    if (Base::has_vel_) {
        Base::H_ = Eigen::Matrix<DataType,Base::meas_space_dim_+tState::Algebra::dim_t_vel_, cov_dim_>::Zero();
        Base::V_ = Eigen::Matrix<DataType,Base::meas_space_dim_ + Base::meas_space_dim_,Base::meas_space_dim_+ tState::Algebra::dim_t_vel_>::Identity();
    } else {
        Base::H_ = Eigen::Matrix<DataType, Base::meas_space_dim_, cov_dim_>::Zero();
        Base::V_ = Eigen::Matrix<DataType,Base::meas_space_dim_,Base::meas_space_dim_>::Identity();
    }


}

//-----------------------------------------------------------------
template <typename tState, MeasurementTypes tMeasurementType, template <typename > typename tTransformation>
Eigen::Matrix<typename tState::DataType,Eigen::Dynamic,Eigen::Dynamic> SourceSENPosVel<tState,tMeasurementType,tTransformation>::DerivedGetLinObsMatState(const tState& state)  {

    MatXd H = Base::H_;

    if (Base::has_vel_) {
        H.block(0,0,Base::meas_space_dim_, Base::meas_space_dim_) = state.g_.R_; //dt/dt
        H.block(Base::meas_space_dim_, tState::Group::dim_, Base::meas_space_dim_,1) = state.g_.R_.block(0,0,tState::Algebra::dim_t_vel_,1); //dtd/rho_x
        if ( Base::meas_space_dim_== 2) { 
            H.block(Base::meas_space_dim_,Base::meas_space_dim_, Base::meas_space_dim_, tState::Group::RotAlgebra::dim_) = state.g_.R_ * tState::Group::RotAlgebra::Wedge(Eigen::Matrix<DataType,tState::Group::RotAlgebra::dim_,1>::Ones()) * state.u_.p_;
        } else {
            H.block(Base::meas_space_dim_,Base::meas_space_dim_, Base::meas_space_dim_, tState::Group::RotAlgebra::dim_) = -state.g_.R_ * tState::Group::RotAlgebra::Wedge(state.u_.p_.block(0,0,tState::Group::RotAlgebra::dim_,1));
        }
    } else {
        H.block(0,0,Base::meas_space_dim_, Base::meas_space_dim_) = state.g_.R_;
    }

    return H;

}

//-----------------------------------------------------------------

template <typename tState, MeasurementTypes tMeasurementType, template <typename > typename tTransformation>
Eigen::Matrix<typename tState::DataType,Eigen::Dynamic,Eigen::Dynamic> SourceSENPosVel<tState,tMeasurementType,tTransformation>::DerivedGetLinObsMatSensorNoise(const tState& state) {

  return Base::V_;

}


//-----------------------------------------------------------------
template <typename tState, MeasurementTypes tMeasurementType, template <typename > typename tTransformation>
Meas<typename tState::DataType> SourceSENPosVel<tState,tMeasurementType,tTransformation>::DerivedGetEstMeas(const tState& state) {
    Meas<DataType> m;
    m.pose = state.g_.t_;
    m.type = tMeasurementType;

    if(Base::has_vel_)
        m.twist = state.g_.R_*state.u_.p_;

    return m;
} 

//-----------------------------------------------------------------
template <typename tState, MeasurementTypes tMeasurementType, template <typename > typename tTransformation>
Eigen::Matrix<typename tState::DataType,Eigen::Dynamic,Eigen::Dynamic> SourceSENPosVel<tState,tMeasurementType,tTransformation>::DerivedOMinus(const Meas<DataType>& m1, const Meas<DataType>& m2) {

    Eigen::Matrix<DataType, Base::meas_space_dim_*Base::meas_space_dim_mult_,1> error;

    if (Base::has_vel_) {
        error.block(0,0,Base::meas_space_dim_,1) = m1.pose - m2.pose;
        error.block(Base::meas_space_dim_,0,Base::meas_space_dim_,1) = m1.twist - m2.twist;
    } else {
        error = m1.pose - m2.pose;
    }

    return error;

}

//----------------------------------------------------------------------------------------
template <typename tState, MeasurementTypes tMeasurementType, template <typename > typename tTransformation>
Meas<typename tState::DataType> SourceSENPosVel<tState,tMeasurementType,tTransformation>::DerivedGenerateRandomMeasurement(const MatXd& meas_std, const tState& state) const {
    Meas<DataType> m;
    m.source_index = this->params_.source_index_;
    m.type = tMeasurementType;

    MatXd deviation = meas_std*utilities::GaussianRandomGenerator(meas_std.rows());

    m.pose = state.g_.t_ + deviation.block(0,0,Base::meas_space_dim_,1); 

    if (Base::has_vel_) {
        m.twist = state.g_.R_*state.u_.p_ + deviation.block(Base::meas_space_dim_,0,Base::meas_space_dim_,1);
    } 

    return m;
}


} // namespace rransac

#endif // RRANSAC_COMMON_SOURCES_SOURCE_SEN_POS_VEL_H_
