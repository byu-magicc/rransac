#ifndef RRANSAC_COMMON_SOURCES_SOURCE_SE3_CAM_DEPTH_H_
#define RRANSAC_COMMON_SOURCES_SOURCE_SE3_CAM_DEPTH_H_
#pragma once


#include <typeinfo>

#include "rransac/common/sources/source_base.h"
#include "rransac/common/utilities.h"
#include "lie_groups/utilities.h"
#include "lie_groups/state.h"

namespace rransac {
using namespace utilities;

/**
 * \class SourceSE3CamDepth
 * This source is meant to be used for a target that evolves on SE3 and that is observed 
 * by a camera and another source that measures depth. Thus, the measurement space is 
 * (R2)x(R1). The camera measurement is expressed on the normalized image sphere. This
 * source is compatible with MeasurementType::CamDepth
 */ 
template <typename _State, MeasurementTypes _MeasurementType, template <typename > typename _Transformation>
class SourceSE3CamDepth: public SourceBase<SourceDerivedTraits<_State,MatXT<typename _State::DataType>, MatXT<typename _State::DataType>,_Transformation,4,1,3,1,4,3,true,_MeasurementType,CompatibleWithModelSENPosVel>, SourceSE3CamDepth> {

public:

typedef SourceBase<SourceDerivedTraits<_State,MatXT<typename _State::DataType>, MatXT<typename _State::DataType>,_Transformation,4,1,3,1,4,3,true,_MeasurementType,CompatibleWithModelSENPosVel>, SourceSE3CamDepth> Base;


typedef typename Base::State State;                                            /**< The state of the target. @see State. */
typedef typename Base::DataType DataType;                                      /**< The scalar object for the data. Ex. float, double, etc. */
typedef typename Base::MatH MatH;                                              /**< The object type of the Jacobian H. */
typedef typename Base::MatV MatV;                                              /**< The object type of the Jacobians V. */
typedef typename Base::Transformation Transformation;                          /**< The transformation used to transform the measurements and tracks. */
typedef typename Base::MatMeasCov MatMeasCov;                                  /**< The data type of the measurement covariance. */
typedef typename Base::VecMeas VecMeas;                                        /**< The data type of the measurement covariance. */
static constexpr unsigned int meas_pose_rows_  = Base::meas_pose_rows_;        /**< The number of rows in the pose measurement. */
static constexpr unsigned int meas_pose_cols_  = Base::meas_pose_cols_;        /**< The number of columns in the pose measurement. */
static constexpr unsigned int meas_twist_rows_ = Base::meas_twist_rows_;       /**< The number of rows in the twist measurement. */
static constexpr unsigned int meas_twist_cols_ = Base::meas_twist_cols_;       /**< The number of columns in the twist measurement. */
static constexpr unsigned int meas_pose_dim_   = Base::meas_pose_dim_;         /**< The measurement pose dimension. */
static constexpr unsigned int meas_twist_dim_  = Base::meas_twist_dim_;        /**< The measurement twist dimension. */
static constexpr unsigned int total_meas_dim_  = Base::total_meas_dim_;        /**< The total measurement dimension. */
static constexpr bool has_vel_ = Base::has_vel_;                               /**< Indicates if the measurement contains velocity.  */
static constexpr MeasurementTypes measurement_type_ = Base::measurement_type_; /**< The measurement type of the source. */
typedef typename Base::ModelCompatibility ModelCompatibility;                  /**< Indicates which model the source is compatible with. */
typedef typename Base::TransformDataType TransformDataType;                    /**< The error type of the difference between two measurements. */
typedef typename Base::Measurement Measurement;                                /**< The measurement data type. */

// static constexpr unsigned int l_dim_ =  State::Algebra::dim_a_vel_ + 1;         /**< The dimension of the angular velocity of the target plus one. */
static constexpr unsigned int cov_dim_ = _State::Group::dim_ + _State::Algebra::dim_ - _State::Algebra::dim_t_vel_ + 1; /**< The dimension of the state covariance. */


static_assert(std::is_same<_State,lie_groups::SE3_se3>::value, "SourceSE3CamDepth: The state is not compatible with the source model");
static_assert( measurement_type_ == MeasurementTypes::SE3_CAM_DEPTH, "SourceSE3CamDepth: The measurement type is not compatible with the source."    );


/**
 * Initializes the Jacobians
 */ 
SourceSE3CamDepth();

/** 
 * Initializes the measurement source. Currently it does nothing. 
 * @param[in] params The source parameters.
 */
void DerivedInit(const SourceParameters& params){}   


/** 
 * Returns the jacobian of the observation function w.r.t. the states.
 * @param[in] state The state of a target at which the Jacobian is be evaluated.
 */
static MatH DerivedGetLinObsMatState(const State& state);

                      

/** 
 * Returns the jacobian of the observation function w.r.t. the sensor noise.
 * @param[in] state The state of a target at which the Jacobian is be evaluated.
 */
static MatV DerivedGetLinObsMatSensorNoise(const State& state);



/**
 *  Implements the observation function and returns an estimated measurement based on the state. 
 * Currently, the measurement is only given a pose, twist, and measurement type. 
 * @param[in] state A state of the target.
 */
static Measurement DerivedGetEstMeas(const State& state);

/**
 * Performs the OMinus operation between two measurement (m1 ominus m2) of the same type. In other words, this
 * method computes the geodesic distance between two measurements of the same type.
 * @param[in] m1 a measurement
 * @param[in] m2 a measurement
 */
static VecMeas DerivedOMinus(const Measurement& m1, const Measurement& m2);



/**
 * Generates a random measurement from a Gaussian distribution with mean defined by the state and standard deviation of the measurement covariance
 * @param[in] state The state that serves as the mean in the Gaussian distribution
 * @param[in] meas_std The measurement standard deviation
 */ 
Measurement DerivedGenerateRandomMeasurement(const MatMeasCov& meas_std, const State& state) const;


};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                            Definitions
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template <typename _State, MeasurementTypes _MeasurementType, template <typename > typename _Transformation>
SourceSE3CamDepth<_State,_MeasurementType,_Transformation>::SourceSE3CamDepth() {

    // Verify measurement type
    Base::V_ = Eigen::Matrix<DataType, total_meas_dim_,total_meas_dim_>::Identity();
    Base::H_ = Eigen::Matrix<DataType, total_meas_dim_, cov_dim_>::Zero();
}

//-----------------------------------------------------------------
template <typename _State, MeasurementTypes _MeasurementType, template <typename > typename _Transformation>
typename SourceSE3CamDepth<_State,_MeasurementType,_Transformation>::MatH SourceSE3CamDepth<_State,_MeasurementType,_Transformation>::DerivedGetLinObsMatState(const State& state) {

    static constexpr unsigned int dim_pos =  State::Group::dim_pos_;
    static constexpr unsigned int dim_rot =  State::Group::dim_rot_;
    static constexpr unsigned int dim_t_vel =  State::Algebra::dim_t_vel_;
    static constexpr unsigned int dim_a_vel =  State::Algebra::dim_a_vel_;


    MatH H = Base::H_;
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

template <typename _State, MeasurementTypes _MeasurementType, template <typename > typename _Transformation>
typename SourceSE3CamDepth<_State,_MeasurementType,_Transformation>::MatV SourceSE3CamDepth<_State,_MeasurementType,_Transformation>::DerivedGetLinObsMatSensorNoise(const State& state) {

    return Base::V_;

}

//-----------------------------------------------------------------
template <typename _State, MeasurementTypes _MeasurementType, template <typename > typename _Transformation>
typename SourceSE3CamDepth<_State,_MeasurementType,_Transformation>::Measurement SourceSE3CamDepth<_State,_MeasurementType,_Transformation>::DerivedGetEstMeas(const State& state) {
    
    static constexpr unsigned int dim_pos =  State::Group::dim_pos_;
    static constexpr unsigned int dim_rot =  State::Group::dim_rot_;
    static constexpr unsigned int dim_t_vel =  State::Algebra::dim_t_vel_;
    static constexpr unsigned int dim_a_vel =  State::Algebra::dim_a_vel_;

    const Eigen::Matrix<DataType,dim_pos,1>& t = state.g_.t_;
    const Eigen::Matrix<DataType,dim_t_vel,1>& p = state.u_.p_;
    const Eigen::Matrix<DataType,dim_rot,dim_rot>& R = state.g_.R_;

    Measurement m;
    m.type = measurement_type_;
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
template <typename _State, MeasurementTypes _MeasurementType, template <typename > typename _Transformation>
typename SourceSE3CamDepth<_State,_MeasurementType,_Transformation>::VecMeas SourceSE3CamDepth<_State,_MeasurementType,_Transformation>::DerivedOMinus(const Measurement& m1, const Measurement& m2) {

    VecMeas error;
    error.block(0,0,4,1) = m1.pose - m2.pose;
    error.block(4,0,3,1) = m1.twist - m2.twist;


    return error;

}

//----------------------------------------------------------------------------------------
template <typename _State, MeasurementTypes _MeasurementType, template <typename > typename _Transformation>
typename SourceSE3CamDepth<_State,_MeasurementType,_Transformation>::Measurement SourceSE3CamDepth<_State,_MeasurementType,_Transformation>::DerivedGenerateRandomMeasurement(const MatMeasCov& meas_std,const State& state) const {
    Measurement m = this->DerivedGetEstMeas(state);
    m.source_index = this->params_.source_index_;
    m.type = measurement_type_;

    VecMeas deviation = meas_std*utilities::GaussianRandomGenerator(meas_std.rows());

    m.pose += deviation.block(0,0,4,1);
    m.pose.block(1,0,3,1)/= m.pose.block(1,0,3,1).norm();
    m.twist += deviation.block(4,0,3,1);
    

    return m;
}



} // namespace rransac

#endif // RRANSAC_COMMON_SOURCES_SOURCE_SE3_CAM_DEPTH_H_
