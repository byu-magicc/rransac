#ifndef RRANSAC_COMMON_SOURCES_SOURCE_SEN_POS_VEL_H_
#define RRANSAC_COMMON_SOURCES_SOURCE_SEN_POS_VEL_H_
#pragma once


#include <typeinfo>

#include "rransac/common/sources/source_base.h"
#include "rransac/common/utilities.h"
#include "lie_groups/utilities.h"

namespace rransac {


/**
 * \class SourceSE3CamDepth
 * This source is meant to be used for a target that evolves on SE3 and that is observed 
 * by a camera and another source that measures depth. Thus, the measurement space is 
 * (R2)x(R1). The camera measurement is expressed on the normalized image sphere. This
 * source is compatible with MeasurementType::CamDepth
 */ 
template<class tState>
class SourceSE3CamDepth: public SourceBase<tState,SourceSE3CamDepth<tState>> {

public:


typedef tState State;                                                            /**< The state of the target. @see State. */
typedef typename tState::DataType DataType;                                      /**< The scalar object for the data. Ex. float, double, etc. */
typedef Eigen::Matrix<DataType,Eigen::Dynamic,Eigen::Dynamic> MatXd;             /**< The object type of the Jacobians. */
static constexpr unsigned int meas_space_dim_ = 5;                               /**< The dimension of the measurement space. */
// static constexpr unsigned int meas_pose_rows_ = tState::Group::dim_pos_;         /**< The number of rows in the pose measurement. */
// static constexpr unsigned int meas_pose_cols_ = 1;                               /**< The number of columns in the pose measurement. */
// static constexpr unsigned int meas_twist_rows_ = tState::Group::dim_pos_;        /**< The number of rows in the twist measurement. */
// static constexpr unsigned int meas_twist_cols_ = 1;                              /**< The number of columns in the twist measurement. */
typedef utilities::CompatibleWithModelSENPosVel ModelCompatibility;              /**< Indicates which model this source is compatible with. */


// static constexpr unsigned int l_dim_ =  tState::Algebra::dim_a_vel_ + 1;         /**< The dimension of the angular velocity of the target plus one. */
static constexpr unsigned int cov_dim_ = tState::Group::dim_ + tState::Algebra::dim_ - tState::Algebra::dim_t_vel_ + 1; /**< The dimension of the state covariance. */


static_assert(std::is_same<tState,liegroups::SE3_se3>, "The state is not compatible with the source model");



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
void SourceSE3CamDepth<tState>::DerivedInit(const SourceParameters& params) {

    // Verify measurement type
    this->V_ = Eigen::Matrix<DataType,meas_space_dim_,meas_space_dim_>::Identity();
    this->H_ = Eigen::Matrix<DataType, meas_space_dim_, this->cov_dim_>::Zero();
}

//-----------------------------------------------------------------
template<class tState>
Eigen::Matrix<typename tState::DataType,Eigen::Dynamic,Eigen::Dynamic> SourceSE3CamDepth<tState>::DerivedGetLinObsMatState(const tState& state) const {

    using dim_pos = tState::Group::dim_pos_;
    using dim_rot = tSTate::Group::dim_rot_;
    using dim_t_vel = tState::Algebra::dim_t_vel_;
    using dim_a_vel = tState::Algebra::dim_a_vel_;

    MatXd H = this->H_;
    Eigen::Matrix<DataType,dim_pos,1>& t = state.g_.t_;
    Eigen::Matrix<DataType,dim_t_vel,1>& p = state.u_.p_;
    Eigen::Matrix<DataType,dim_rot,dim_rot>& R = state.g_.R_;
    double d = t.norm();
    double d3 = powf(d,3);
    double d5 = powf(d,5);
    Eigen::Matrix<DataType,dim_pos,dim_pos> ttR = t*t.transpose()*R;
    Eigen::Matrix<DataType,3,3> tmp = R/d - ttR/d3;

    H.block(0,0,dim_pos,dim_pos) = t.transpose()*R/d;
    H.block(dim_pos,0,dim_rot,dim_pos) = R/d - ttR/d3;

    H.block(dim_pos+dim_rot,0,dim_pos,dim_pos) = 3*ttR*p*t.transpose()*R/d5 - (R*p*t.transpose()*R + t.transpose()*p + R*t.transpose()*R*p)/d3;
    H.block(dim_pos+dim_rot,dim_pos,dim_rot,dim_rot) = -tmp*liegroups::se3::SSM(p);
    H.block(6,6,3,1) = tmp.block(0,0,3,1);


}

//-----------------------------------------------------------------
template<class tState>
Meas<typename tState::DataType> SourceSE3CamDepth<tState>::DerivedGetEstMeas(const tState& state) const{
    
    Eigen::Matrix<DataType,dim_pos,1>& t = state.g_.t_;
    Eigen::Matrix<DataType,dim_t_vel,1>& p = state.u_.p_;
    Eigen::Matrix<DataType,dim_rot,dim_rot>& R = state.g_.R_;

    Meas<DataType> m;
    m.type = this->params_.type_;
    m.source_index = this->params_.source_index_;
    m.pose = Eigen::Matrix<double,3,1>::Zero();
    m.twist = Eigen::Matrix<double,2,1>::Zero();
    double d = t.norm();
    m.pose(0,0) = d;
    if (d != 0) {
        m.pose.block(1,0,2,1) = state.g_.t_/d;
        m.twist = R*p/d - t*t.transpose()*R*p/powf(d,3);
    }

    
    return m;
} 

//-----------------------------------------------------------------
template<class tState>
Eigen::Matrix<typename tState::DataType,Eigen::Dynamic,Eigen::Dynamic> SourceSE3CamDepth<tState>::DerivedOMinus(const Meas<DataType>& m1, const Meas<DataType>& m2) {

    Eigen::Matrix<DataType, meas_space_dim_,1> error;
    error.block(0,0,3,1) = m1.pose - m2.pose;
    error.block(0,3,2,1) = m1.twist - m2.twist;

}

//----------------------------------------------------------------------------------------
template<class tState>
Meas<typename tState::DataType> SourceSE3CamDepth<tState>::DerivedGenerateRandomMeasurement(const tState& state, const MatXd& meas_std){
    Meas<DataType> m = this->DerivedGetEstMeas(state);
    m.source_index = this->params_.source_index_;

    MatXd deviation = meas_std*this->GaussianRandomGenerator(meas_std.rows());

    m.pose += deviation.block(0,0,3,1);
    m.twist += deviation.block(3,0,2,1);
    

    return m;
}

//----------------------------------------------------------------------------------------

template<class tState>
Meas<typename tState::DataType> SourceSE3CamDepth<tState>::DerivedGetEstMeas(const tState& state, const MeasurementTypes type) {

    Eigen::Matrix<DataType,dim_pos,1>& t = state.g_.t_;
    Eigen::Matrix<DataType,dim_t_vel,1>& p = state.u_.p_;
    Eigen::Matrix<DataType,dim_rot,dim_rot>& R = state.g_.R_;

    Meas<DataType> m;
    m.type = type;
    m.pose = Eigen::Matrix<double,3,1>::Zero();
    m.twist = Eigen::Matrix<double,2,1>::Zero();
    double d = t.norm();
    m.pose(0,0) = d;
    if (d != 0) {
        m.pose.block(1,0,2,1) = state.g_.t_/d;
        m.twist = R*p/d - t*t.transpose()*R*p/powf(d,3);
    }

    return m;


}

//-------------------------------------------------------------------------------------------------------------------------

template<class tState>
Eigen::Matrix<typename tState::DataType,Eigen::Dynamic,Eigen::Dynamic> SourceSE3CamDepth<tState>::DerivedGetLinObsMatState(const tState& state, const MeasurementTypes type) {

    using dim_pos = tState::Group::dim_pos_;
    using dim_rot = tSTate::Group::dim_rot_;
    using dim_t_vel = tState::Algebra::dim_t_vel_;
    using dim_a_vel = tState::Algebra::dim_a_vel_;
    MatXd H = Eigen::Matrix<double,meas_space_dim_, cov_dim_>::Zero();
    
    Eigen::Matrix<DataType,dim_pos,1>& t = state.g_.t_;
    Eigen::Matrix<DataType,dim_t_vel,1>& p = state.u_.p_;
    Eigen::Matrix<DataType,dim_rot,dim_rot>& R = state.g_.R_;
    double d = t.norm();
    double d3 = powf(d,3);
    double d5 = powf(d,5);
    Eigen::Matrix<DataType,dim_pos,dim_pos> ttR = t*t.transpose()*R;
    Eigen::Matrix<DataType,3,3> tmp = R/d - ttR/d3;

    H.block(0,0,dim_pos,dim_pos) = t.transpose()*R/d;
    H.block(dim_pos,0,dim_rot,dim_pos) = R/d - ttR/d3;

    H.block(dim_pos+dim_rot,0,dim_pos,dim_pos) = 3*ttR*p*t.transpose()*R/d5 - (R*p*t.transpose()*R + t.transpose()*p + R*t.transpose()*R*p)/d3;
    H.block(dim_pos+dim_rot,dim_pos,dim_rot,dim_rot) = -tmp*liegroups::se3::SSM(p);
    H.block(6,6,3,1) = tmp.block(0,0,3,1);
   
    return H;

}

//-------------------------------------------------------------------------------------------------------------------------

template<class tState>
Eigen::Matrix<typename tState::DataType,Eigen::Dynamic,Eigen::Dynamic> SourceSE3CamDepth<tState>::DerivedGetLinObsMatSensorNoise(const tState& state, const MeasurementTypes type) {

    const unsigned int sizeg = tState::Group::dim_;
    const unsigned int sizeu = tState::Algebra::dim_;

    MatXd V= Eigen::Matrix<DataType,meas_space_dim_,meas_space_dim_>::Identity();


    return V;

}

// Common Sources
typedef SourceSE3CamDepth<lie_groups::SE2_se2> SourceSE2PosVel;
typedef SourceSE3CamDepth<lie_groups::SE3_se3> SourceSE3PosVel;

} // namespace rransac

#endif // RRANSAC_COMMON_SOURCES_SOURCE_SEN_POS_VEL_H_
