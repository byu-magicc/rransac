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

template<class S>
class SourceRN : public SourceBase<S,SourceRN<S>> {

public:

typedef S type_;
static constexpr unsigned int dim = S::g_type_::dim_;

SourceRN()=default;
~SourceRN()=default;

/** Initializes the measurement source. This function must set the parameters.  */
void Init(const SourceParameters& params);      

/** Returns the jacobian of the observation function w.r.t. the states */
Eigen::MatrixXd GetLinObsMatState(const S& state) {return this->H_;}                        

/** Returns the jacobian of the observation function w.r.t. the sensor noise */
Eigen::MatrixXd GetLinObsMatSensorNoise(const S& state) {return this->V_;}                         

/** Computes the estimated measurement given a state */
Meas GetEstMeas(const S& state) {
    Meas m;
    m.pose = state.g_.data_;
    m.twist = state.u_.data_;
    return m;
    } 

/**
 * Returns the error between the estimated measurement and the measurement
 */
Eigen::MatrixXd OMinus(const Meas& m1, const Meas& m2) {

    if (this->params_.type_ == MeasurementTypes::RN_POS) {
        return m1.pose - m2.pose;
    } else if (this->params_.type_ == MeasurementTypes::RN_POS_VEL){
        Eigen::Matrix<double, S::g_type_::dim_*2,1> error;
        error.block(0,0,S::g_type_::dim_,1) = m1.pose - m2.pose;
        error.block(S::g_type_::dim_,0,S::g_type_::dim_,1) = m1.twist - m2.twist;
        return error;
    } else {
        throw std::runtime_error("SourceRN::OMinus Measurement type not supported.");
    }


}

/**
 * Maps the pose to Euclidean space. In this case, it just returns the pose.
 * @param Meas The measurement whose pose needs to be transformed
 */
Eigen::MatrixXd ToEuclidean(const Meas& m)  {
    return m.pose;
}

/**
 * Generates a random measurement from a Gaussian distribution with mean defined by the state and covariance defined by meas_cov
 * @param state The state that serves as the mean in the Gaussian distribution
 * @param meas_std The measurement standard deviation
 */ 
Meas GenerateRandomMeasurement(const S& state, const Eigen::MatrixXd& meas_std){
    Meas m;
    m.source_index = this->params_.source_index_;



    Eigen::MatrixXd deviation = meas_std*this->GaussianRandomGenerator(meas_std.rows());
    // Eigen::MatrixXd deviation = meas_std*this->GaussianRandomGenerator(5);

 

    switch (this->params_.type_)
    {
    case MeasurementTypes::RN_POS:        
        m.pose = S::g_type_::OPlus(state.g_.data_,deviation);
        m.type = MeasurementTypes::RN_POS;
        break;
    case MeasurementTypes::RN_POS_VEL:
        m.pose = S::g_type_::OPlus(state.g_.data_, deviation.block(0,0,S::g_type_::dim_,1));
        m.twist = state.u_.data_ + deviation.block(S::g_type_::dim_,0,S::g_type_::dim_,1);
        m.type = MeasurementTypes::RN_POS_VEL;
        break;
    default:
        throw std::runtime_error("SourceRN::GenerateRandomMeasurement Measurement type not supported.");
        break;
    }

    return m;
}

};


//--------------------------------------------------

template<class S>
void SourceRN<S>::Init(const SourceParameters& params) {

    const unsigned int sizeg = S::g_type_::dim_;
    const unsigned int sizeu = S::u_type_::dim_;

    // Make sure that the state is a valid type
    if (typeid(S).name() != typeid(lie_groups::State<lie_groups::Rn<sizeg>>).name())
    {
        throw std::runtime_error("SourceRNPos::Init State is not supported by this source type");
    }

    // Construct the Jacobians
    switch (params.type_)
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

    this->params_ = params;

}

// template<class S>
// Eigen::MatrixXd SourceRN<S>::GetLinObsMatState(S const& state) {return this->H_;}


} // namesapce rransac


#endif // RRANSAC_COMMON_SOURCES_SOURCE_RN_H_
