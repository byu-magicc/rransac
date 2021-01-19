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

template <class tState, typename tDataType = double>
class SourceSENPoseTwist : public SourceBase<tState,SourceSENPoseTwist<tState>> {

public:

typedef tState State;
typedef tDataType DataType;
typedef Eigen::Matrix<tDataType,Eigen::Dynamic,Eigen::Dynamic> MatXd;
static constexpr unsigned int meas_dim_ = tState::Group::dim_;


/** Initializes the measurement source. This function must set the parameters.  */
void DerivedInit(const SourceParameters& params);      

/** Returns the jacobian of the observation function w.r.t. the states 
 * @param state A state of a model.
*/
MatXd DerivedGetLinObsMatState(tState const& state) const {return this->H_;};  

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
static MatXd DerivedOMinus(const Meas<tDataType>& m1, const Meas<tDataType>& m2) ;

/**
 * Maps the pose to Euclidean space. The translation is unchanged; however, the rotation is transformed using Cayley coordinates of the first kind.
 * @param m The measurement whose pose needs to be transformed
 */
MatXd DerivedToEuclidean(const Meas<tDataType>& m);

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

//-----------------------------------------------------------------

template <class tState, typename tDataType >
void SourceSENPoseTwist<tState, tDataType>::DerivedInit(const SourceParameters& params) {

    // Verify state
    if (typeid(State).name() != typeid(lie_groups::SE2_se2).name() && typeid(State).name() != typeid(lie_groups::SE3_se3).name())
    {
        throw std::runtime_error("SourceSENPoseTwist::Init State is not supported by this source type");
    }

    // Verify measurement type
    switch (this->params_.type_)
    {
    case MeasurementTypes::SEN_POSE:
        this->V_ = Eigen::Matrix<tDataType,meas_dim_,meas_dim_>::Identity();
        this->H_ = Eigen::Matrix<tDataType, meas_dim_, State::dim_>::Zero();
        this->H_.block(0,0,meas_dim_,meas_dim_).setIdentity();
        break;
    case MeasurementTypes::SEN_POSE_TWIST:
        this->V_ = Eigen::Matrix<tDataType,meas_dim_*2,meas_dim_*2>::Identity();
        this->H_ = Eigen::Matrix<tDataType,meas_dim_*2,State::dim_>::Identity();
        break;
    default:
        throw std::runtime_error("SourceSENPoseTwist::Init Measurement type not supported.");
        break;
    }

}

//----------------------------------------------------------------------------------------
template <class tState, typename tDataType >
Meas<tDataType> SourceSENPoseTwist<tState, tDataType>::DerivedGetEstMeas(const tState& state) const {
    Meas<tDataType> m;
    m.type = this->params_.type_;
    m.pose = state.g_.data_;
    m.twist = state.u_.data_;
    return m;
} 

//----------------------------------------------------------------------------------------
template <class tState, typename tDataType >
Eigen::Matrix<tDataType,Eigen::Dynamic,Eigen::Dynamic> SourceSENPoseTwist<tState, tDataType>::DerivedOMinus(const Meas<tDataType>& m1, const Meas<tDataType>& m2) {

    if (m1.type == MeasurementTypes::SEN_POSE && m1.type == m2.type) {
        return State::Group::OMinus(m1.pose,m2.pose);
    } else if (m1.type == MeasurementTypes::SEN_POSE_TWIST && m1.type == m2.type){
        Eigen::Matrix<tDataType, meas_dim_*2,1> error;
        error.block(0,0,meas_dim_,1) = State::Group::OMinus(m1.pose,m2.pose);
        error.block(meas_dim_,0,meas_dim_,1) = m1.twist - m2.twist;
        return error;
    } else {
        throw std::runtime_error("SourceSENPoseTwist::OMinus Measurement type not supported.");
    }
}

//----------------------------------------------------------------------------------------
template <class tState, typename tDataType >
Eigen::Matrix<tDataType,Eigen::Dynamic,Eigen::Dynamic> SourceSENPoseTwist<tState, tDataType>::DerivedToEuclidean(const Meas<tDataType>& m)  {
    typename State::Group g(m.pose);
    Eigen::Matrix<tDataType,  State::Group::dim_pos_,State::Group::dim_pos_> I = Eigen::Matrix<tDataType,  State::Group::dim_pos_,State::Group::dim_pos_>::Identity();
    Eigen::Matrix<tDataType, State::Group::dim_pos_,State::Group::dim_pos_>&& u = 2.0*(g.R_-I)*(g.R_+I).inverse();
    Eigen::Matrix<tDataType, State::Group::dim_,1> pose_euclidean;
    pose_euclidean.block(0,0,State::Group::dim_pos_,1) = g.t_;
    pose_euclidean.block(State::Group::dim_pos_,0,State::Group::dim_rot_,1) = State::Group::rot_algebra::Vee(u);


    return pose_euclidean;
}

//----------------------------------------------------------------------------------------
template <class tState, typename tDataType >
Meas<tDataType> SourceSENPoseTwist<tState, tDataType>::DerivedGenerateRandomMeasurement(const tState& state, const MatXd& meas_std){
    Meas<tDataType> m;
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

template<class tState, typename tDataType >
Meas<tDataType> SourceSENPoseTwist<tState, tDataType>::DerivedGetEstMeas(const tState& state, const MeasurementTypes type) {

    Meas<tDataType> m;
    m.type = type;
    m.pose = state.g_.data_;
    if (MeasurementTypes::SEN_POSE_TWIST == type)
        m.twist = state.u_.data_;
    return m;

}

//----------------------------------------------------------------------------------------

template<class tState, typename tDataType >
Eigen::Matrix<tDataType,Eigen::Dynamic,Eigen::Dynamic> SourceSENPoseTwist<tState, tDataType>::DerivedGetLinObsMatState(const tState& state, const MeasurementTypes type) {

    const unsigned int sizeg = tState::Group::dim_;
    const unsigned int sizeu = tState::Algebra::dim_;

    // Construct the Jacobians
    MatXd H;
    switch (type)
    {
    case MeasurementTypes::SEN_POSE:
        H = Eigen::Matrix<tDataType, meas_dim_, State::dim_>::Zero();
        H.block(0,0,meas_dim_,meas_dim_).setIdentity();
        break;
    case MeasurementTypes::SEN_POSE_TWIST:
        H = Eigen::Matrix<tDataType,meas_dim_*2,State::dim_>::Identity();
        break;
    default:
        throw std::runtime_error("SourceSENPoseTwist::Init Measurement type not supported.");
        break;
    }
    return H;

}

//-------------------------------------------------------------------------------------------------------------------------

template<class tState, typename tDataType >
Eigen::Matrix<tDataType,Eigen::Dynamic,Eigen::Dynamic> SourceSENPoseTwist<tState, tDataType>::DerivedGetLinObsMatSensorNoise(const tState& state, const MeasurementTypes type) {

    const unsigned int sizeg = tState::Group::dim_;
    const unsigned int sizeu = tState::Algebra::dim_;

    MatXd V;

    // Construct the Jacobians
    switch (type)
    {
    case MeasurementTypes::SEN_POSE:
        V = Eigen::Matrix<tDataType,meas_dim_,meas_dim_>::Identity();
        break;
    case MeasurementTypes::SEN_POSE_TWIST:
        V = Eigen::Matrix<tDataType,meas_dim_*2,meas_dim_*2>::Identity();
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