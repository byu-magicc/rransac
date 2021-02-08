#ifndef RRANSAC_COMMON_SOURCES_SOURCE_SEN_POSE_TWIST_H_
#define RRANSAC_COMMON_SOURCES_SOURCE_SEN_POSE_TWIST_H_
#pragma once


#include <typeinfo>

#include "rransac/common/sources/source_base.h"
#include "rransac/common/utilities.h"
#include "lie_groups/utilities.h"


namespace rransac {

/**
 * \class SourceSENPoseTwist
 * This source is meant to be used for a target that evolves on SEN and whose
 * measurements is also SEN. Compatible with 
 * MeasurementType::SEN_POSE and MeasurementType::SEN_POSE_TWIST.
 * It is assumed that the partial derivative of the observation function w.r.t.
 * the measurement noise is identity.
 */ 

template <class tState>
class SourceSENPoseTwist : public SourceBase<tState,SourceSENPoseTwist<tState>> {

public:

typedef tState State;                                                               /**< The state of the target. @see State. */
typedef typename tState::DataType DataType;                                         /**< The scalar object for the data. Ex. float, double, etc. */
typedef Eigen::Matrix<DataType,Eigen::Dynamic,Eigen::Dynamic> MatXd;                /**< The object type of the Jacobians. */
static constexpr unsigned int meas_space_dim_ = State::Group::dim_;                 /**< The dimension of the measurement space. */
static constexpr unsigned int meas_pose_rows_ = State::Group::size1_;               /**< The number of rows in the pose measurement. */
static constexpr unsigned int meas_pose_cols_ = State::Group::size2_;               /**< The number of columns in the pose measurement. */
static constexpr unsigned int meas_twist_rows_ = State::Group::dim_;                /**< The number of rows in the twist measurement. */
static constexpr unsigned int meas_twist_cols_ = 1;                                 /**< The number of columns in the twist measurement. */
typedef utilities::CompatibleWithModelSENPoseTwist ModelCompatibility;              /**< Indicates which model this source is compatible with. */

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
MatXd DerivedGetLinObsMatState(tState const& state) const {return this->H_;};  

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
static MatXd DerivedGetLinObsMatSensorNoise(const tState& state, const MeasurementTypes type);

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
static MatXd DerivedOMinus(const Meas<DataType>& m1, const Meas<DataType>& m2) ;


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

//-----------------------------------------------------------------

template <class tState>
void SourceSENPoseTwist<tState>::DerivedInit(const SourceParameters& params) {


    // Verify measurement type
    switch (this->params_.type_)
    {
    case MeasurementTypes::SEN_POSE:
        this->V_ = Eigen::Matrix<DataType,meas_space_dim_,meas_space_dim_>::Identity();
        this->H_ = Eigen::Matrix<DataType, meas_space_dim_, State::dim_>::Zero();
        this->H_.block(0,0,meas_space_dim_,meas_space_dim_).setIdentity();
        break;
    case MeasurementTypes::SEN_POSE_TWIST:
        this->V_ = Eigen::Matrix<DataType,meas_space_dim_*2,meas_space_dim_*2>::Identity();
        this->H_ = Eigen::Matrix<DataType,meas_space_dim_*2,State::dim_>::Identity();
        break;
    default:
        throw std::runtime_error("SourceSENPoseTwist::Init Measurement type not supported.");
        break;
    }

}

//----------------------------------------------------------------------------------------
template <class tState>
Meas<typename tState::DataType> SourceSENPoseTwist<tState>::DerivedGetEstMeas(const tState& state) const {
    Meas<DataType> m;
    m.type = this->params_.type_;
    m.pose = state.g_.data_;
    m.twist = state.u_.data_;
    return m;
} 

//----------------------------------------------------------------------------------------
template <class tState>
Eigen::Matrix<typename tState::DataType,Eigen::Dynamic,Eigen::Dynamic> SourceSENPoseTwist<tState>::DerivedOMinus(const Meas<DataType>& m1, const Meas<DataType>& m2) {

    if (m1.type == MeasurementTypes::SEN_POSE && m1.type == m2.type) {
        return State::Group::OMinus(m1.pose,m2.pose);
    } else if (m1.type == MeasurementTypes::SEN_POSE_TWIST && m1.type == m2.type){
        Eigen::Matrix<DataType, meas_space_dim_*2,1> error;
        error.block(0,0,meas_space_dim_,1) = State::Group::OMinus(m1.pose,m2.pose);
        error.block(meas_space_dim_,0,meas_space_dim_,1) = m1.twist - m2.twist;
        return error;
    } else {
        throw std::runtime_error("SourceSENPoseTwist::OMinus Measurement type not supported.");
    }
}

//----------------------------------------------------------------------------------------
template <class tState>
Meas<typename tState::DataType> SourceSENPoseTwist<tState>::DerivedGenerateRandomMeasurement(const tState& state, const MatXd& meas_std){
    Meas<DataType> m;
    m.source_index = this->params_.source_index_;

    MatXd deviation = meas_std*this->GaussianRandomGenerator(meas_std.rows());

    switch (this->params_.type_)
    {
    case MeasurementTypes::SEN_POSE:        
        m.pose = State::Group::OPlus(state.g_.data_,deviation);
        m.type = MeasurementTypes::SEN_POSE;
        break;
    case MeasurementTypes::SEN_POSE_TWIST:
        m.pose = State::Group::OPlus(state.g_.data_, deviation.block(0,0,State::Group::dim_,1));
        m.twist = state.u_.data_ + deviation.block(State::Group::dim_,0,State::Group::dim_,1);
        m.type = MeasurementTypes::SEN_POSE_TWIST;

        break;
    default:
        throw std::runtime_error("SourceSENPoseTwist::GenerateRandomMeasurement Measurement type not supported.");
        break;
    }

    return m;
}

//----------------------------------------------------------------------------------------

template<class tState >
Meas<typename tState::DataType> SourceSENPoseTwist<tState>::DerivedGetEstMeas(const tState& state, const MeasurementTypes type) {

    Meas<DataType> m;
    m.type = type;
    m.pose = state.g_.data_;
    if (MeasurementTypes::SEN_POSE_TWIST == type)
        m.twist = state.u_.data_;
    return m;

}

//----------------------------------------------------------------------------------------

template<class tState>
Eigen::Matrix<typename tState::DataType,Eigen::Dynamic,Eigen::Dynamic> SourceSENPoseTwist<tState>::DerivedGetLinObsMatState(const tState& state, const MeasurementTypes type) {

    const unsigned int sizeg = tState::Group::dim_;
    const unsigned int sizeu = tState::Algebra::dim_;

    // Construct the Jacobians
    MatXd H;
    switch (type)
    {
    case MeasurementTypes::SEN_POSE:
        H = Eigen::Matrix<DataType, meas_space_dim_, State::dim_>::Zero();
        H.block(0,0,meas_space_dim_,meas_space_dim_).setIdentity();
        break;
    case MeasurementTypes::SEN_POSE_TWIST:
        H = Eigen::Matrix<DataType,meas_space_dim_*2,State::dim_>::Identity();
        break;
    default:
        throw std::runtime_error("SourceSENPoseTwist::Init Measurement type not supported.");
        break;
    }
    return H;

}

//-------------------------------------------------------------------------------------------------------------------------

template<class tState >
Eigen::Matrix<typename tState::DataType,Eigen::Dynamic,Eigen::Dynamic> SourceSENPoseTwist<tState>::DerivedGetLinObsMatSensorNoise(const tState& state, const MeasurementTypes type) {

    const unsigned int sizeg = tState::Group::dim_;
    const unsigned int sizeu = tState::Algebra::dim_;

    MatXd V;

    // Construct the Jacobians
    switch (type)
    {
    case MeasurementTypes::SEN_POSE:
        V = Eigen::Matrix<DataType,meas_space_dim_,meas_space_dim_>::Identity();
        break;
    case MeasurementTypes::SEN_POSE_TWIST:
        V = Eigen::Matrix<DataType,meas_space_dim_*2,meas_space_dim_*2>::Identity();
        break;
    default:
        throw std::runtime_error("SourceSENPoseTwist::Init Measurement type not supported.");
        break;
    }

    return V;

}

// Common Sources
typedef SourceSENPoseTwist<lie_groups::SE2_se2> SourceSE2PoseTwist;
typedef SourceSENPoseTwist<lie_groups::SE3_se3> SourceSE3PoseTwist;



} // namespace rransac
#endif // RRANSAC_COMMON_SOURCES_SOURCE_SEN_POSE_TWIST_H_