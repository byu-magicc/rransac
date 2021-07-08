#ifndef RRANSAC_COMMON_SOURCES_SOURCE_RN_H_
#define RRANSAC_COMMON_SOURCES_SOURCE_RN_H_
#pragma once

#include <typeinfo>

#include "rransac/common/sources/source_base.h"
#include "rransac/common/utilities.h"
#include "lie_groups/utilities.h"
#include "lie_groups/state.h"

namespace rransac {

using namespace utilities;

/**
 * \class SourceRN
 * This source is meant to be used for a target that evolves on RN, whose
 * measurements space is RN. The Jacobian of the observation function w.r.t. the
 * noise is assumed identity.
 */ 

template <typename _State, MeasurementTypes _MeasurementType, template <typename > typename _Transformation>
class SourceRN : public SourceBase<SourceDerivedTraits<_State,MatXT<typename _State::DataType>,MatXT<typename _State::DataType>,_Transformation,_State::Group::dim_,1,MeasHasVelocity<_MeasurementType>::value ? _State::Algebra::dim_ : 0 , MeasHasVelocity<_MeasurementType>::value ? 1 : 0,_State::Group::dim_,MeasHasVelocity<_MeasurementType>::value ? _State::Algebra::dim_ : 0, MeasHasVelocity<_MeasurementType>::value, _MeasurementType,CompatibleWithModelRN>, SourceRN> {

public:


typedef SourceBase<SourceDerivedTraits<_State,MatXT<typename _State::DataType>,MatXT<typename _State::DataType>,_Transformation,_State::Group::dim_,1,MeasHasVelocity<_MeasurementType>::value ? _State::Algebra::dim_ : 0 , MeasHasVelocity<_MeasurementType>::value ? 1 : 0,_State::Group::dim_,MeasHasVelocity<_MeasurementType>::value ? _State::Algebra::dim_ : 0, MeasHasVelocity<_MeasurementType>::value, _MeasurementType,CompatibleWithModelRN>, SourceRN> Base;

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
typedef typename Base::Measurement Measurement;                                            /**< The measurement data type. */



// Perform compatibility checks
static_assert(lie_groups::utilities::StateIsRN_rN<_State>::value, "SourceRN: The state is not compatible with the model");
static_assert( _MeasurementType == MeasurementTypes::RN_POS || _MeasurementType==MeasurementTypes::RN_POS_VEL, "SourceRN: The measurement type is not compatible with the source."    );


/**
 * Initializes the Jacobians
 */ 
SourceRN();
~SourceRN()=default;

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
SourceRN<_State, _MeasurementType, _Transformation>::SourceRN() {


    Base::H_ = Eigen::Matrix<DataType, total_meas_dim_, State::dim_>::Zero();
    Base::H_.block(0,0, total_meas_dim_, total_meas_dim_).setIdentity();
    Base::V_ = Eigen::Matrix<DataType, total_meas_dim_,total_meas_dim_>::Identity();
    
}

//-------------------------------------------------------------------------------------------------------------------------

template <typename _State, MeasurementTypes _MeasurementType, template <typename > typename _Transformation>
typename SourceRN<_State, _MeasurementType, _Transformation>::MatH SourceRN<_State, _MeasurementType, _Transformation>::DerivedGetLinObsMatState(const State& state) {

    return Base::H_;
}

//-------------------------------------------------------------------------------------------------------------------------

template <typename _State, MeasurementTypes _MeasurementType, template <typename > typename _Transformation>
typename SourceRN<_State, _MeasurementType, _Transformation>::MatV SourceRN<_State, _MeasurementType, _Transformation>::DerivedGetLinObsMatSensorNoise(const State& state) {

    return Base::V_;
}

//---------------------------------------------------------------------------

template <typename _State, MeasurementTypes _MeasurementType, template <typename > typename _Transformation>
typename SourceRN<_State,_MeasurementType, _Transformation>::Measurement SourceRN<_State,_MeasurementType, _Transformation>::DerivedGetEstMeas(const State& state) {
    Measurement m;
    m.type = measurement_type_;
    m.pose = state.g_.data_;
    if (has_vel_) {
        m.twist = state.u_.data_.block(0,0,meas_pose_dim_,1);
    }
    return m;
} 

//---------------------------------------------------------------------------
template <typename _State, MeasurementTypes _MeasurementType, template <typename > typename _Transformation>
typename SourceRN<_State,_MeasurementType, _Transformation>::VecMeas SourceRN<_State,_MeasurementType, _Transformation>::DerivedOMinus(const Measurement& m1, const Measurement& m2) {


    VecMeas error;
    error.block(0,0,meas_pose_dim_,1) = m1.pose - m2.pose;

    if(has_vel_) {
        error.block(meas_pose_dim_,0,meas_pose_dim_,1) = m1.twist - m2.twist;
    }
    
    return error;
}

//---------------------------------------------------------------------------------------------
template <typename _State, MeasurementTypes _MeasurementType, template <typename > typename _Transformation>
typename SourceRN<_State,_MeasurementType, _Transformation>::Measurement SourceRN<_State,_MeasurementType, _Transformation>::DerivedGenerateRandomMeasurement(const MatMeasCov& meas_std, const State& state) const {
    Measurement m;

    VecMeas deviation = meas_std*utilities::GaussianRandomGenerator(meas_std.rows());

    m = DerivedGetEstMeas(state);
    m.source_index = this->params_.source_index_;

    m.pose += deviation.block(0,0,meas_pose_dim_,1);

    if(has_vel_) {
        m.twist += deviation.block(meas_pose_dim_,0,meas_pose_dim_,1);
    }

    return m;
}









} // namespace rransac


#endif // RRANSAC_COMMON_SOURCES_SOURCE_RN_H_
