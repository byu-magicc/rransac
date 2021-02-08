#ifndef RRANSAC_COMMON_SOURCES_SOURCE_RN_H_
#define RRANSAC_COMMON_SOURCES_SOURCE_RN_H_
#pragma once

#include <typeinfo>

#include "rransac/common/sources/source_base.h"
#include "rransac/common/utilities.h"
#include "rransac/utilities.h"

namespace rransac {


/**
 * \class SourceRN
 * This source is meant to be used for a target that evolves on RN, whose
 * measurements space is RN. The Jacobian of the observation function w.r.t. the
 * noise is assumed identity.
 */ 

template<class tState>
class SourceRN : public SourceBase<tState, SourceRN<tState>> {

public:

typedef tState State;                                                       /**< The state of the target. @see State. */
typedef typename tState::DataType DataType;                                 /**< The scalar object for the data. Ex. float, double, etc. */
typedef Eigen::Matrix<DataType,Eigen::Dynamic,Eigen::Dynamic> MatXd;        /**< The object type of the Jacobians. */
static constexpr unsigned int meas_space_dim_ = tState::Group::dim_;        /**< The dimension of the measurement space. */
static constexpr unsigned int meas_pose_rows_ = tState::Group::dim_;        /**< The number of rows in the pose measurement. */
static constexpr unsigned int meas_pose_cols_ = 1;                          /**< The number of columns in the pose measurement. */
static constexpr unsigned int meas_twist_rows_ = tState::Group::dim_;       /**< The number of rows in the twist measurement. */
static constexpr unsigned int meas_twist_cols_ = 1;                         /**< The number of columns in the twist measurement. */
typedef utilities::CompatibleWithModelRN ModelCompatibility;                /**< Indicates which model this source is compatible with. */


static_assert(lie_groups::utilities::StateIsRN_rN<tState>::value, "The state is not compatible with the model");


SourceRN()=default;
~SourceRN()=default;

/** 
 * Initializes the measurement source by initializing the Jacobians. 
 * @param[in] params The source parameters.
 */
void DerivedInit(const SourceParameters& params);      

/** 
 * Returns the jacobian of the observation function w.r.t. the states. This is an optimized version. 
 * @param[in] state The state of a target at which the Jacobian is be evaluated.
 */
MatXd DerivedGetLinObsMatState(const tState& state) const {return this->H_;}   

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
void SourceRN<tState>::DerivedInit(const SourceParameters& params) {

    const unsigned int sizeg = tState::Group::dim_;
    const unsigned int sizeu = tState::Algebra::dim_;

    // Construct the Jacobians
    switch (this->params_.type_)
    {
    case MeasurementTypes::RN_POS:
        this->H_ = Eigen::Matrix<DataType,sizeg,sizeg+sizeu>::Zero();
        this->H_.block(0,0,sizeg,sizeg).setIdentity();
        this->V_ = Eigen::Matrix<DataType,sizeg,sizeg>::Identity();
        break;
    case MeasurementTypes::RN_POS_VEL:
        this->H_ = Eigen::Matrix<DataType,sizeg+sizeu,sizeg+sizeu>::Identity();
        this->V_ = Eigen::Matrix<DataType,sizeg+sizeu,sizeg+sizeu>::Identity();
        break;
    default:
        throw std::runtime_error("SourceRN::Init Measurement type not supported.");
        break;
    }
    
}

//---------------------------------------------------------------------------

template<class tState>
Meas<typename tState::DataType> SourceRN<tState>::DerivedGetEstMeas(const tState& state) const {
    Meas<DataType> m;
    m.type = this->params_.type_;
    m.pose = state.g_.data_;
    m.twist = state.u_.data_;
    return m;
} 

//---------------------------------------------------------------------------
template<class tState>
Eigen::Matrix<typename tState::DataType,Eigen::Dynamic,Eigen::Dynamic> SourceRN<tState>::DerivedOMinus(const Meas<DataType>& m1, const Meas<DataType>& m2) {

    if (m1.type == MeasurementTypes::RN_POS && m2.type == m1.type) {
        return m1.pose - m2.pose;
    } else if (m1.type == MeasurementTypes::RN_POS_VEL && m2.type == m1.type){
        Eigen::Matrix<DataType, tState::Group::dim_*2,1> error;
        error.block(0,0,tState::Group::dim_,1) = m1.pose - m2.pose;
        error.block(tState::Group::dim_,0,tState::Group::dim_,1) = m1.twist - m2.twist;
        return error;
    } else {
        throw std::runtime_error("SourceRN::OMinus Measurement type not supported.");
    }
}

//---------------------------------------------------------------------------------------------
template<class tState>
Meas<typename tState::DataType> SourceRN<tState>::DerivedGenerateRandomMeasurement(const tState& state, const MatXd& meas_std){
    Meas<DataType> m;
    m.source_index = this->params_.source_index_;

    MatXd deviation = meas_std*this->GaussianRandomGenerator(meas_std.rows());

    switch (this->params_.type_)
    {
    case MeasurementTypes::RN_POS:        
        m.pose = tState::Group::OPlus(state.g_.data_,deviation);
        m.type = MeasurementTypes::RN_POS;
        break;
    case MeasurementTypes::RN_POS_VEL:
        m.pose = tState::Group::OPlus(state.g_.data_, deviation.block(0,0,tState::Group::dim_,1));
        m.twist = state.u_.data_ + deviation.block(tState::Group::dim_,0,tState::Group::dim_,1);
        m.type = MeasurementTypes::RN_POS_VEL;
        break;
    default:
        throw std::runtime_error("SourceRN::GenerateRandomMeasurement Measurement type not supported.");
        break;
    }

    return m;
}

//-------------------------------------------------------------------------------------------------------------------------

template<class tState>
Meas<typename tState::DataType> SourceRN<tState>::DerivedGetEstMeas(const State& state, const MeasurementTypes type) {

    Meas<DataType> m;
    m.type = type;
    m.pose = state.g_.data_;
    if (MeasurementTypes::RN_POS_VEL == type)
        m.twist = state.u_.data_;
    return m;

}

//-------------------------------------------------------------------------------------------------------------------------

template<class tState>
Eigen::Matrix<typename tState::DataType,Eigen::Dynamic,Eigen::Dynamic> SourceRN<tState>::DerivedGetLinObsMatState(const State& state, const MeasurementTypes type) {

    const unsigned int sizeg = tState::Group::dim_;
    const unsigned int sizeu = tState::Algebra::dim_;

    // Construct the Jacobians
    MatXd H;
    switch (type)
    {
    case MeasurementTypes::RN_POS:
        H = Eigen::Matrix<DataType,sizeg,sizeg+sizeu>::Zero();
        H.block(0,0,sizeg,sizeg).setIdentity();
        break;
    case MeasurementTypes::RN_POS_VEL:
        H = Eigen::Matrix<DataType,sizeg+sizeu,sizeg+sizeu>::Identity();
        break;
    default:
        throw std::runtime_error("SourceRN::Init Measurement type not supported.");
        break;
    }
    return H;

}

//-------------------------------------------------------------------------------------------------------------------------

template<class tState>
Eigen::Matrix<typename tState::DataType,Eigen::Dynamic,Eigen::Dynamic> SourceRN<tState>::DerivedGetLinObsMatSensorNoise(const State& state, const MeasurementTypes type) {

    const unsigned int sizeg = tState::Group::dim_;
    const unsigned int sizeu = tState::Algebra::dim_;

    MatXd V;

    // Construct the Jacobians
    switch (type)
    {
    case MeasurementTypes::RN_POS:
        V = Eigen::Matrix<DataType,sizeg,sizeg>::Identity();
        break;
    case MeasurementTypes::RN_POS_VEL:
        V = Eigen::Matrix<DataType,sizeg+sizeu,sizeg+sizeu>::Identity();
        break;
    default:
        throw std::runtime_error("SourceRN::Init Measurement type not supported.");
        break;
    }

    return V;

}

//Common Sources
typedef SourceRN<lie_groups::R2_r2> SourceR2;
typedef SourceRN<lie_groups::R3_r3> SourceR3;


} // namesapce rransac


#endif // RRANSAC_COMMON_SOURCES_SOURCE_RN_H_
