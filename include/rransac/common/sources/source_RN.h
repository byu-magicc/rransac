#ifndef RRANSAC_COMMON_SOURCES_SOURCE_RN_H_
#define RRANSAC_COMMON_SOURCES_SOURCE_RN_H_

#include "common/sources/source_base.h"
#include <typeinfo>

namespace rransac {


/**
 * \class SourceRN
 * This source is meant to be used for a target that evolves on RN and whose
 * measurements are on RN or on RN and its tangent space. Compatible with 
 * MeasurementType::RN_POS and MeasurementType::RN_VEL
 */ 

template<class tState, typename tDataType = double>
class SourceRN : public SourceBase<tState, SourceRN<tState,tDataType>, tDataType> {

public:

typedef tState State;
typedef Eigen::Matrix<tDataType,Eigen::Dynamic,Eigen::Dynamic> MatXd;
static constexpr unsigned int meas_dim_ = tState::Group::dim_;

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
Meas<tDataType> DerivedGetEstMeas(const tState& state) const ;

/** Computes the estimated measurement given a state and measurement type. This method is not optimized.*/
static Meas<tDataType> DerivedGetEstMeas(const State& state, const MeasurementTypes type);

/**
 * Returns the error between the estimated measurement and the measurement
 */
static MatXd DerivedOMinus(const Meas<tDataType>& m1, const Meas<tDataType>& m2);



/**
 * Generates a random measurement from a Gaussian distribution with mean defined by the state and covariance defined by meas_cov
 * @param state The state that serves as the mean in the Gaussian distribution
 * @param meas_std The measurement standard deviation
 */ 
Meas<tDataType> DerivedGenerateRandomMeasurement(const tState& state, const MatXd& meas_std);

};




//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                            Definitions
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class tState, typename tDataType>
void SourceRN<tState, tDataType>::DerivedInit(const SourceParameters& params) {

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
        this->H_ = Eigen::Matrix<tDataType,sizeg,sizeg+sizeu>::Zero();
        this->H_.block(0,0,sizeg,sizeg).setIdentity();
        this->V_ = Eigen::Matrix<tDataType,sizeg,sizeg>::Identity();
        break;
    case MeasurementTypes::RN_POS_VEL:
        this->H_ = Eigen::Matrix<tDataType,sizeg+sizeu,sizeg+sizeu>::Identity();
        this->V_ = Eigen::Matrix<tDataType,sizeg+sizeu,sizeg+sizeu>::Identity();
        break;
    default:
        throw std::runtime_error("SourceRN::Init Measurement type not supported.");
        break;
    }
    
}

//---------------------------------------------------------------------------

template<class tState, typename tDataType>
Meas<tDataType> SourceRN<tState, tDataType>::DerivedGetEstMeas(const tState& state) const {
    Meas<tDataType> m;
    m.type = this->params_.type_;
    m.pose = state.g_.data_;
    m.twist = state.u_.data_;
    return m;
} 

//---------------------------------------------------------------------------
template<class tState, typename tDataType>
Eigen::Matrix<tDataType,Eigen::Dynamic,Eigen::Dynamic> SourceRN<tState, tDataType>::DerivedOMinus(const Meas<tDataType>& m1, const Meas<tDataType>& m2) {

    if (m1.type == MeasurementTypes::RN_POS && m2.type == m1.type) {
        return m1.pose - m2.pose;
    } else if (m1.type == MeasurementTypes::RN_POS_VEL && m2.type == m1.type){
        Eigen::Matrix<double, tState::Group::dim_*2,1> error;
        error.block(0,0,tState::Group::dim_,1) = m1.pose - m2.pose;
        error.block(tState::Group::dim_,0,tState::Group::dim_,1) = m1.twist - m2.twist;
        return error;
    } else {
        throw std::runtime_error("SourceRN::OMinus Measurement type not supported.");
    }
}

//---------------------------------------------------------------------------------------------
template<class tState, typename tDataType >
Meas<tDataType> SourceRN<tState, tDataType>::DerivedGenerateRandomMeasurement(const tState& state, const MatXd& meas_std){
    Meas<tDataType> m;
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

template<class tState, typename tDataType >
Meas<tDataType> SourceRN<tState, tDataType>::DerivedGetEstMeas(const State& state, const MeasurementTypes type) {

    Meas<tDataType> m;
    m.type = type;
    m.pose = state.g_.data_;
    if (MeasurementTypes::RN_POS_VEL == type)
        m.twist = state.u_.data_;
    return m;

}

//-------------------------------------------------------------------------------------------------------------------------

template<class tState, typename tDataType >
Eigen::Matrix<tDataType,Eigen::Dynamic,Eigen::Dynamic> SourceRN<tState, tDataType>::DerivedGetLinObsMatState(const State& state, const MeasurementTypes type) {

    const unsigned int sizeg = tState::Group::dim_;
    const unsigned int sizeu = tState::Algebra::dim_;

    // Construct the Jacobians
    MatXd H;
    switch (type)
    {
    case MeasurementTypes::RN_POS:
        H = Eigen::Matrix<tDataType,sizeg,sizeg+sizeu>::Zero();
        H.block(0,0,sizeg,sizeg).setIdentity();
        break;
    case MeasurementTypes::RN_POS_VEL:
        H = Eigen::Matrix<tDataType,sizeg+sizeu,sizeg+sizeu>::Identity();
        break;
    default:
        throw std::runtime_error("SourceRN::Init Measurement type not supported.");
        break;
    }
    return H;

}

//-------------------------------------------------------------------------------------------------------------------------

template<class tState, typename tDataType >
Eigen::Matrix<tDataType,Eigen::Dynamic,Eigen::Dynamic> SourceRN<tState, tDataType>::DerivedGetLinObsMatSensorNoise(const State& state, const MeasurementTypes type) {

    const unsigned int sizeg = tState::Group::dim_;
    const unsigned int sizeu = tState::Algebra::dim_;

    MatXd V;

    // Construct the Jacobians
    switch (type)
    {
    case MeasurementTypes::RN_POS:
        V = Eigen::Matrix<tDataType,sizeg,sizeg>::Identity();
        break;
    case MeasurementTypes::RN_POS_VEL:
        V = Eigen::Matrix<tDataType,sizeg+sizeu,sizeg+sizeu>::Identity();
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
