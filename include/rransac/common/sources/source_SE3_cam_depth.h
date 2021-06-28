#ifndef RRANSAC_COMMON_SOURCES_SOURCE_SE3_CAM_DEPTH_H_
#define RRANSAC_COMMON_SOURCES_SOURCE_SE3_CAM_DEPTH_H_
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
template <typename tState, MeasurementTypes tMeasurementType, template <typename > typename tTransformation>
class SourceSE3CamDepth: public SourceBase<tState,tMeasurementType, tTransformation<tState>, SourceSE3CamDepth<tState,tMeasurementType,tTransformation>> {

public:


typedef tState State;                                                            /**< The state of the target. @see State. */
typedef typename tState::DataType DataType;                                      /**< The scalar object for the data. Ex. float, double, etc. */
typedef Eigen::Matrix<DataType,Eigen::Dynamic,Eigen::Dynamic> MatXd;             /**< The object type of the Jacobians. */
static constexpr unsigned int meas_space_dim_ = 7;                               /**< The dimension of the measurement space. */
static constexpr unsigned int meas_pose_rows_ = 4;                               /**< The number of rows in the pose measurement. */
static constexpr unsigned int meas_pose_cols_ = 1;                               /**< The number of columns in the pose measurement. */
static constexpr unsigned int meas_twist_rows_ = 3;                              /**< The number of rows in the twist measurement. */
static constexpr unsigned int meas_twist_cols_ = 1;                              /**< The number of columns in the twist measurement. */
static constexpr MeasurementTypes measurement_type_ = tMeasurementType;          /**< The measurement type of the source. */
typedef utilities::CompatibleWithModelSENPosVel ModelCompatibility;              /**< Indicates which model this source is compatible with. */
typedef SourceBase<tState,tMeasurementType, tTransformation<tState>, SourceSE3CamDepth<tState,tMeasurementType,tTransformation>> Base;                          /**< The source base class. */


static constexpr int dim_mult_ = 1;    /**< a constant used when the measurement contains velocity. */
static constexpr int has_vel_ = true; /**< Indicates if the measurement contains velocity.  */

// static constexpr unsigned int l_dim_ =  tState::Algebra::dim_a_vel_ + 1;         /**< The dimension of the angular velocity of the target plus one. */
static constexpr unsigned int cov_dim_ = tState::Group::dim_ + tState::Algebra::dim_ - tState::Algebra::dim_t_vel_ + 1; /**< The dimension of the state covariance. */


static_assert(std::is_same<tState,lie_groups::SE3_se3>::value, "SourceSE3CamDepth: The state is not compatible with the source model");
static_assert( tMeasurementType == MeasurementTypes::SE3_CAM_DEPTH, "SourceSE3CamDepth: The measurement type is not compatible with the source."    );



/** 
 * Initializes the measurement source by initializing the Jacobians. 
 * @param[in] params The source parameters.
 */
void DerivedInit(const SourceParameters& params);      


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
void SourceSE3CamDepth<tState,tMeasurementType,tTransformation>::DerivedInit(const SourceParameters& params) {

    // Verify measurement type
    Base::V_ = Eigen::Matrix<DataType,meas_space_dim_,meas_space_dim_>::Identity();
    Base::H_ = Eigen::Matrix<DataType, meas_space_dim_, this->cov_dim_>::Zero();
}

//-----------------------------------------------------------------
template <typename tState, MeasurementTypes tMeasurementType, template <typename > typename tTransformation>
Eigen::Matrix<typename tState::DataType,Eigen::Dynamic,Eigen::Dynamic> SourceSE3CamDepth<tState,tMeasurementType,tTransformation>::DerivedGetLinObsMatState(const tState& state) {

    static constexpr unsigned int dim_pos =  tState::Group::dim_pos_;
    static constexpr unsigned int dim_rot =  tState::Group::dim_rot_;
    static constexpr unsigned int dim_t_vel =  tState::Algebra::dim_t_vel_;
    static constexpr unsigned int dim_a_vel =  tState::Algebra::dim_a_vel_;


    MatXd H = Base::H_;
    const Eigen::Matrix<DataType,dim_pos,1>& t = state.g_.t_;
    const Eigen::Matrix<DataType,dim_t_vel,1>& p = state.u_.p_;
    const Eigen::Matrix<DataType,dim_rot,dim_rot>& R = state.g_.R_;
    double d = t.norm();
    double d3 = pow(d,3);
    double d5 = pow(d,5);
    Eigen::Matrix<DataType,3,3> ttR = t*t.transpose()*R;
    Eigen::Matrix<DataType,3,3> tmp = R/d - ttR/d3;

    H.block(0,0,1,3) = t.transpose()*R/d;
    H.block(1,0,3,3) = R/d - ttR/d3;

    H.block(4,0,3,3) = 3*ttR*(p*t.transpose())*R/d5 - (R*(p*t.transpose())*R + t*p.transpose() + R*(t.transpose()*R*p))/d3;
    H.block(4,3,3,3) = -tmp*lie_groups::se3<DataType>::SSM(p);
    H.block(4,6,3,1) = tmp.block(0,0,3,1);

    return H;


}

//-------------------------------------------------------------------------------------------------------------------------

template <typename tState, MeasurementTypes tMeasurementType, template <typename > typename tTransformation>
Eigen::Matrix<typename tState::DataType,Eigen::Dynamic,Eigen::Dynamic> SourceSE3CamDepth<tState,tMeasurementType,tTransformation>::DerivedGetLinObsMatSensorNoise(const tState& state) {

    return Base::V_;

}

//-----------------------------------------------------------------
template <typename tState, MeasurementTypes tMeasurementType, template <typename > typename tTransformation>
Meas<typename tState::DataType> SourceSE3CamDepth<tState,tMeasurementType,tTransformation>::DerivedGetEstMeas(const tState& state) {
    
    static constexpr unsigned int dim_pos =  tState::Group::dim_pos_;
    static constexpr unsigned int dim_rot =  tState::Group::dim_rot_;
    static constexpr unsigned int dim_t_vel =  tState::Algebra::dim_t_vel_;
    static constexpr unsigned int dim_a_vel =  tState::Algebra::dim_a_vel_;

    const Eigen::Matrix<DataType,dim_pos,1>& t = state.g_.t_;
    const Eigen::Matrix<DataType,dim_t_vel,1>& p = state.u_.p_;
    const Eigen::Matrix<DataType,dim_rot,dim_rot>& R = state.g_.R_;

    Meas<DataType> m;
    m.type = tMeasurementType;
    m.pose = Eigen::Matrix<double,4,1>::Zero();
    m.twist = Eigen::Matrix<double,3,1>::Zero();
    double d = t.norm();
    m.pose(0,0) = d;
    if (d != 0) {
        m.pose.block(1,0,3,1) = state.g_.t_/d;
        m.twist = R*p/d - t*t.transpose()*R*p/pow(d,3);
    }

    
    return m;
} 

//-----------------------------------------------------------------
template <typename tState, MeasurementTypes tMeasurementType, template <typename > typename tTransformation>
Eigen::Matrix<typename tState::DataType,Eigen::Dynamic,Eigen::Dynamic> SourceSE3CamDepth<tState,tMeasurementType,tTransformation>::DerivedOMinus(const Meas<DataType>& m1, const Meas<DataType>& m2) {

    Eigen::Matrix<DataType, meas_space_dim_,1> error;
    error.block(0,0,4,1) = m1.pose - m2.pose;
    error.block(4,0,3,1) = m1.twist - m2.twist;


    return error;

}

//----------------------------------------------------------------------------------------
template <typename tState, MeasurementTypes tMeasurementType, template <typename > typename tTransformation>
Meas<typename tState::DataType> SourceSE3CamDepth<tState,tMeasurementType,tTransformation>::DerivedGenerateRandomMeasurement(const MatXd& meas_std,const tState& state) const {
    Meas<DataType> m = this->DerivedGetEstMeas(state);
    m.source_index = this->params_.source_index_;
    m.type = tMeasurementType;

    MatXd deviation = meas_std*utilities::GaussianRandomGenerator(meas_std.rows());

    m.pose += deviation.block(0,0,4,1);
    m.pose.block(1,0,3,1)/= m.pose.block(1,0,3,1).norm();
    m.twist += deviation.block(4,0,3,1);
    

    return m;
}



} // namespace rransac

#endif // RRANSAC_COMMON_SOURCES_SOURCE_SE3_CAM_DEPTH_H_
