#ifndef RRANSAC_COMMON_SOURCES_SOURCE_RN_H_
#define RRANSAC_COMMON_SOURCES_SOURCE_RN_H_
#pragma once


#include "common/sources/source_base.h"
#include "common/utilities.h"
#include <typeinfo>

namespace rransac {


/**
 * \class SourceRN
 * This source is meant to be used for a target that evolves on RN and whose
 * measurements are on RN or on RN and its tangent space. Compatible with 
 * MeasurementType::RN_POS and MeasurementType::RN_VEL
 */ 

template<class tState>
class SourceRN : public SourceBase<tState, SourceRN<tState>> {

public:

typedef tState State;
typedef typename tState::DataType DataType;
typedef Eigen::Matrix<DataType,Eigen::Dynamic,Eigen::Dynamic> MatXd;
static constexpr unsigned int meas_dim_ = tState::Group::dim_;
static constexpr unsigned int meas_pose_rows_ = tState::Group::dim_;
static constexpr unsigned int meas_pose_cols_ = 1;
static constexpr unsigned int meas_twist_rows_ = tState::Group::dim_;
static constexpr unsigned int meas_twist_cols_ = 1;
typedef utilities::CompatibleWithModelRN ModelCompatibility;
static_assert(lie_groups::utilities::StateIsRN_rN<tState>::value, "The state is not compatible with the model");


SourceRN()=default;
~SourceRN()=default;

/** Initializes the measurement source. This function must set the parameters.  */
void DerivedInit(const SourceParameters& params);      

/** Returns the jacobian of the observation function w.r.t. the states */
MatXd DerivedGetLinObsMatState(const tState& state) const {return this->H_;}   

/** Returns the jacobian of the observation function w.r.t. the states. This method is not optimized. */
static MatXd DerivedGetLinObsMatState(const State& state, const MeasurementTypes type);

/** Returns the jacobian of the observation function w.r.t. the sensor noise */
MatXd DerivedGetLinObsMatSensorNoise(const tState& state) const {return this->V_;}                         

/** Returns the jacobian of the observation function w.r.t. the sensor noise. This method is not optimized */
static MatXd DerivedGetLinObsMatSensorNoise(const State& state, const MeasurementTypes type);

/** Computes the estimated measurement given a state */
Meas<DataType> DerivedGetEstMeas(const tState& state) const ;

/** Computes the estimated measurement given a state and measurement type. This method is not optimized.*/
static Meas<DataType> DerivedGetEstMeas(const State& state, const MeasurementTypes type);

/**
 * Returns the error between the estimated measurement and the measurement
 */
static MatXd DerivedOMinus(const Meas<DataType>& m1, const Meas<DataType>& m2);



/**
 * Generates a random measurement from a Gaussian distribution with mean defined by the state and covariance defined by meas_cov
 * @param state The state that serves as the mean in the Gaussian distribution
 * @param meas_std The measurement standard deviation
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

    // Make sure that the state is a valid type
    if (typeid(tState).name() != typeid(lie_groups::State<lie_groups::Rn,double,sizeg>).name())
    {
        throw std::runtime_error("SourceRNPos::Init State is not supported by this source type");
    }

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
