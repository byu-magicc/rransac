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

template <class tState>
class SourceSENPoseTwist : public SourceBase<tState,SourceSENPoseTwist<tState>> {

public:

typedef tState State;
static constexpr unsigned int meas_dim_ = tState::g_type_::dim_;


/** Initializes the measurement source. This function must set the parameters.  */
void DerivedInit(const SourceParameters& params);      

/** Returns the jacobian of the observation function w.r.t. the states 
 * @param state A state of a model.
*/
Eigen::MatrixXd DerivedGetLinObsMatState(tState const& state) const {return this->H_;};                        

/** Returns the jacobian of the observation function w.r.t. the sensor noise */
Eigen::MatrixXd DerivedGetLinObsMatSensorNoise(const tState& state) const {return this->V_;}                        

/** Computes the estimated measurement given a state */
Meas DerivedGetEstMeas(const tState& state) const ;

/**
 * Returns the error between the estimated measurement and the measurement
 */
static Eigen::MatrixXd DerivedOMinus(const Meas& m1, const Meas& m2) ;

/**
 * Maps the pose to Euclidean space. The translation is unchanged; however, the rotation is transformed using Cayley coordinates of the first kind.
 * @param m The measurement whose pose needs to be transformed
 */
Eigen::MatrixXd DerivedToEuclidean(const Meas& m);

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

//-----------------------------------------------------------------

template<class tState>
void SourceSENPoseTwist<tState>::DerivedInit(const SourceParameters& params) {

    // Verify state
    if (typeid(State).name() != typeid(lie_groups::SE2_se2).name() && typeid(State).name() != typeid(lie_groups::SE3_se3).name())
    {
        throw std::runtime_error("SourceSENPoseTwist::Init State is not supported by this source type");
    }

    // Verify measurement type
    switch (this->params_.type_)
    {
    case MeasurementTypes::SEN_POSE:
        this->V_ = Eigen::Matrix<double,meas_dim_,meas_dim_>::Identity();
        this->H_ = Eigen::Matrix<double, meas_dim_, State::dim_>::Zero();
        this->H_.block(0,0,meas_dim_,meas_dim_).setIdentity();
        break;
    case MeasurementTypes::SEN_POSE_TWIST:
        this->V_ = Eigen::Matrix<double,meas_dim_*2,meas_dim_*2>::Identity();
        this->H_ = Eigen::Matrix<double,meas_dim_*2,State::dim_>::Identity();
        break;
    default:
        throw std::runtime_error("SourceSENPoseTwist::Init Measurement type not supported.");
        break;
    }

}

//----------------------------------------------------------------------------------------
template<class tState>
Meas SourceSENPoseTwist<tState>::DerivedGetEstMeas(const tState& state) const {
    Meas m;
    m.pose = state.g_.data_;
    m.twist = state.u_.data_;
    return m;
} 

//----------------------------------------------------------------------------------------
template<class tState>
Eigen::MatrixXd SourceSENPoseTwist<tState>::DerivedOMinus(const Meas& m1, const Meas& m2) {

    if (m1.type == MeasurementTypes::SEN_POSE && m1.type == m2.type) {
        return State::g_type_::OMinus(m1.pose,m2.pose);
    } else if (m1.type == MeasurementTypes::SEN_POSE_TWIST && m1.type == m2.type){
        Eigen::Matrix<double, meas_dim_*2,1> error;
        error.block(0,0,meas_dim_,1) = State::g_type_::OMinus(m1.pose,m2.pose);
        error.block(meas_dim_,0,meas_dim_,1) = m1.twist - m2.twist;
        return error;
    } else {
        throw std::runtime_error("SourceSENPoseTwist::OMinus Measurement type not supported.");
    }
}

//----------------------------------------------------------------------------------------
template<class tState>
Eigen::MatrixXd SourceSENPoseTwist<tState>::DerivedToEuclidean(const Meas& m)  {
    typename State::g_type_ g(m.pose);
    Eigen::Matrix<double,  State::g_type_::dim_pos_,State::g_type_::dim_pos_> I = Eigen::Matrix<double,  State::g_type_::dim_pos_,State::g_type_::dim_pos_>::Identity();
    Eigen::Matrix<double, State::g_type_::dim_pos_,State::g_type_::dim_pos_>&& u = 2.0*(g.R_-I)*(g.R_+I).inverse();
    Eigen::Matrix<double, State::g_type_::dim_,1> pose_euclidean;
    pose_euclidean.block(0,0,State::g_type_::dim_pos_,1) = g.t_;
    pose_euclidean.block(State::g_type_::dim_pos_,0,State::g_type_::dim_rot_,1) = State::g_type_::rot_algebra::Vee(u);


    return pose_euclidean;
}

//----------------------------------------------------------------------------------------
template<class tState>
Meas SourceSENPoseTwist<tState>::DerivedGenerateRandomMeasurement(const tState& state, const Eigen::MatrixXd& meas_std){
    Meas m;
    m.source_index = this->params_.source_index_;

    Eigen::MatrixXd deviation = meas_std*this->GaussianRandomGenerator(meas_std.rows());

    switch (this->params_.type_)
    {
    case MeasurementTypes::SEN_POSE:        
        m.pose = State::g_type_::OPlus(state.g_.data_,deviation);
        m.type = MeasurementTypes::SEN_POSE;
        break;
    case MeasurementTypes::SEN_POSE_TWIST:
        m.pose = State::g_type_::OPlus(state.g_.data_, deviation.block(0,0,State::g_type_::dim_,1));
        m.twist = state.u_.data_ + deviation.block(State::g_type_::dim_,0,State::g_type_::dim_,1);
        m.type = MeasurementTypes::SEN_POSE_TWIST;

        break;
    default:
        throw std::runtime_error("SourceSENPoseTwist::GenerateRandomMeasurement Measurement type not supported.");
        break;
    }

    return m;
}

// Common Sources
typedef SourceSENPoseTwist<lie_groups::SE2_se2> SourceSE2PoseTwist;
typedef SourceSENPoseTwist<lie_groups::SE3_se3> SourceSE3PoseTwist;



} // namespace rransac
#endif // RRANSAC_COMMON_SOURCES_SOURCE_SEN_POSE_TWIST_H_