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

template<class tState>
class SourceRN : public SourceBase<tState, SourceRN<tState>> {

public:

typedef tState State;
static constexpr unsigned int meas_dim_ = tState::g_type_::dim_;

SourceRN()=default;
~SourceRN()=default;

/** Initializes the measurement source. This function must set the parameters.  */
void DerivedInit(const SourceParameters& params);      

/** Returns the jacobian of the observation function w.r.t. the states */
Eigen::MatrixXd DerivedGetLinObsMatState(const tState& state) {return this->H_;}                        

/** Returns the jacobian of the observation function w.r.t. the sensor noise */
Eigen::MatrixXd DerivedGetLinObsMatSensorNoise(const tState& state) {return this->V_;}                         

/** Computes the estimated measurement given a state */
Meas DerivedGetEstMeas(const tState& state);

/**
 * Returns the error between the estimated measurement and the measurement
 */
Eigen::MatrixXd DerivedOMinus(const Meas& m1, const Meas& m2);

/**
 * Maps the pose to Euclidean space. In this case, it just returns the pose.
 * @param Meas The measurement whose pose needs to be transformed
 */
Eigen::MatrixXd DerivedToEuclidean(const Meas& m)  {
    return m.pose;
}

/**
 * Generates a random measurement from a Gaussian distribution with mean defined by the state and covariance defined by meas_cov
 * @param state The state that serves as the mean in the Gaussian distribution
 * @param meas_std The measurement standard deviation
 */ 
Meas DerivedGenerateRandomMeasurement(const tState& state, const Eigen::MatrixXd& meas_std);

};




//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                            Definitions
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class tState>
void SourceRN<tState>::DerivedInit(const SourceParameters& params) {

    const unsigned int sizeg = tState::g_type_::dim_;
    const unsigned int sizeu = tState::u_type_::dim_;

    // Make sure that the state is a valid type
    if (typeid(tState).name() != typeid(lie_groups::State<lie_groups::Rn<sizeg>>).name())
    {
        throw std::runtime_error("SourceRNPos::Init State is not supported by this source type");
    }

    // Construct the Jacobians
    switch (this->params_.type_)
    {
    case MeasurementTypes::RN_POS:
        this->H_ = Eigen::Matrix<double,sizeg,sizeg+sizeu>::Zero();
        this->H_.block(0,0,sizeg,sizeg).setIdentity();
        this->V_ = Eigen::Matrix<double,sizeg,sizeg>::Identity();
        break;
    case MeasurementTypes::RN_POS_VEL:
        this->H_ = Eigen::Matrix<double,sizeg+sizeu,sizeg+sizeu>::Identity();
        this->V_ = Eigen::Matrix<double,sizeg+sizeu,sizeg+sizeu>::Identity();
        break;
    default:
        throw std::runtime_error("SourceRN::Init Measurement type not supported.");
        break;
    }
    
}

//---------------------------------------------------------------------------

template<class tState>
Meas SourceRN<tState>::DerivedGetEstMeas(const tState& state) {
    Meas m;
    m.pose = state.g_.data_;
    m.twist = state.u_.data_;
    return m;
} 

//---------------------------------------------------------------------------
template<class tState>
Eigen::MatrixXd SourceRN<tState>::DerivedOMinus(const Meas& m1, const Meas& m2) {

    if (this->params_.type_ == MeasurementTypes::RN_POS) {
        return m1.pose - m2.pose;
    } else if (this->params_.type_ == MeasurementTypes::RN_POS_VEL){
        Eigen::Matrix<double, tState::g_type_::dim_*2,1> error;
        error.block(0,0,tState::g_type_::dim_,1) = m1.pose - m2.pose;
        error.block(tState::g_type_::dim_,0,tState::g_type_::dim_,1) = m1.twist - m2.twist;
        return error;
    } else {
        throw std::runtime_error("SourceRN::OMinus Measurement type not supported.");
    }
}

//---------------------------------------------------------------------------------------------
template<class tState>
Meas SourceRN<tState>::DerivedGenerateRandomMeasurement(const tState& state, const Eigen::MatrixXd& meas_std){
    Meas m;
    m.source_index = this->params_.source_index_;

    Eigen::MatrixXd deviation = meas_std*this->GaussianRandomGenerator(meas_std.rows());

    switch (this->params_.type_)
    {
    case MeasurementTypes::RN_POS:        
        m.pose = tState::g_type_::OPlus(state.g_.data_,deviation);
        m.type = MeasurementTypes::RN_POS;
        break;
    case MeasurementTypes::RN_POS_VEL:
        m.pose = tState::g_type_::OPlus(state.g_.data_, deviation.block(0,0,tState::g_type_::dim_,1));
        m.twist = state.u_.data_ + deviation.block(tState::g_type_::dim_,0,tState::g_type_::dim_,1);
        m.type = MeasurementTypes::RN_POS_VEL;
        break;
    default:
        throw std::runtime_error("SourceRN::GenerateRandomMeasurement Measurement type not supported.");
        break;
    }

    return m;
}

//Common Sources
typedef SourceRN<lie_groups::R2_r2> SourceR2;
typedef SourceRN<lie_groups::R3_r3> SourceR3;


} // namesapce rransac


#endif // RRANSAC_COMMON_SOURCES_SOURCE_RN_H_
