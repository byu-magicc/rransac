#ifndef RRANSAC_COMMON_SOURCES_SOURCE_R2_R3_RADAR_H_
#define RRANSAC_COMMON_SOURCES_SOURCE_R2_R3_RADAR_H_
#pragma once

#include <typeinfo>

#include "rransac/common/sources/source_base.h"
#include "rransac/common/utilities.h"
#include "lie_groups/utilities.h"
#include "lie_groups/state.h"

namespace rransac {

using namespace utilities;


/**
 * \class SourceRadarR2_R3
 * This source is meant to be used for a target that evolves on R2 or R3 and is observed by 
 * a radar sensor with measurements being in polar coordinates for R2 or spherical coordinates for R3.
 * If the measurement type is R2_R3_RADAR_DEPTH_DERIV, then the source includes the derivative of the range
 * in the measurement model. The measurement pose is given in a 3x1 vector in the order range, azimuth, zenith.
 * \latexinclude radar_r2_r3.tex
 *
 */ 

template <typename _State, MeasurementTypes _MeasurementType, template <typename > typename _Transformation>
class SourceRadarR2_R3 : public SourceBase<SourceDerivedTraits<_State, MatXT<typename _State::DataType>,MatXT<typename _State::DataType>,_Transformation,_State::Group::size1_,_State::Group::size2_, MeasHasVelocity<_MeasurementType>::value ? 1 : 0, MeasHasVelocity<_MeasurementType>::value ? 1:0, _State::Group::dim_, MeasHasVelocity<_MeasurementType>::value ? 1:0, MeasHasVelocity<_MeasurementType>::value, _MeasurementType, utilities::CompatibleWithModelRN>, SourceRadarR2_R3> {

public:

typedef SourceBase<SourceDerivedTraits<_State, MatXT<typename _State::DataType>,MatXT<typename _State::DataType>,_Transformation,_State::Group::size1_,_State::Group::size2_, MeasHasVelocity<_MeasurementType>::value ? 1 : 0, MeasHasVelocity<_MeasurementType>::value ? 1:0, _State::Group::dim_, MeasHasVelocity<_MeasurementType>::value ? 1:0, MeasHasVelocity<_MeasurementType>::value, _MeasurementType, utilities::CompatibleWithModelRN>, SourceRadarR2_R3> Base;

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
static constexpr MeasurementTypes measurement_type_ = Base::measurement_type_;             /**< The measurement type of the source. */
typedef typename Base::ModelCompatibility ModelCompatibility;                  /**< Indicates which model the source is compatible with. */
typedef typename Base::TransformDataType TransformDataType;                    /**< The error type of the difference between two measurements. */
typedef typename Base::Measurement Measurement;                                /**< The measurement data type. */

static constexpr unsigned int state_dim_ = _State::dim_;
static constexpr bool polar_coordinates_ = meas_pose_dim_ == 2 ? true : false;      /**< This variable is true if the measurement is in polar coordinates. */
static constexpr bool spherical_coordinates_ = meas_pose_dim_ == 3 ? true : false;  /**< This value is true if the measurement is in spherical coordinates. */

// Perform compatibility checks
static_assert(lie_groups::utilities::StateIsRN_rN<_State>::value, "SourceRadarR2_R3: The state is not compatible with the model");
static_assert(_State::Group::dim_ != 2 || _State::Group::dim_ != 3, "SourceRadarR2_R3: The group element of the state must be either 2 or 3.");
static_assert( measurement_type_ == MeasurementTypes::R2_R3_RADAR || measurement_type_==MeasurementTypes::R2_R3_RADAR_DEPTH_DERIV, "SourceRadarR2_R3: The measurement type is not compatible with the source."    );


/**
 * Initializes the Jacobians
 */ 
SourceRadarR2_R3();
~SourceRadarR2_R3()=default;

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
SourceRadarR2_R3<_State, _MeasurementType, _Transformation>::SourceRadarR2_R3() {


    Base::H_ = Eigen::Matrix<DataType, total_meas_dim_, state_dim_>::Zero();
    Base::V_ = Eigen::Matrix<DataType, total_meas_dim_, total_meas_dim_>::Identity();
    
}

//-------------------------------------------------------------------------------------------------------------------------

template <typename _State, MeasurementTypes _MeasurementType, template <typename > typename _Transformation>
typename SourceRadarR2_R3<_State, _MeasurementType, _Transformation>::MatH SourceRadarR2_R3<_State, _MeasurementType, _Transformation>::DerivedGetLinObsMatState(const State& state) {

    MatH H =Base::H_;

    DataType r = state.g_.data_.norm();

    // Make sure that the range is sufficiently big; otherwise, H will be the zero matrix. 
    if(r > 1e-4) {

        
        const DataType x = state.g_.data_(0);
        const DataType y = state.g_.data_(1);
        
        // polar coordinates
        if(polar_coordinates_) {

            // DataType azimuth = atan2(y,x);
            DataType tmp = x*x+y*y;
            DataType&& r_sqrd = pow(r,2);
            DataType&& r_3 = pow(r,3.0);

            H(0,0) = x/r;
            H(0,1) = y/r;
            H(1,0) = -y/tmp;
            H(1,1) =  x/tmp;



            if(has_vel_) {
                DataType tmp2 = (state.g_.data_.transpose()*state.u_.data_.block(0,0,meas_pose_rows_,1))(0,0);

                H(2,0) = state.u_.data_(0,0)/r - x*tmp2/r_3;
                H(2,1) = state.u_.data_(1,0)/r - y*tmp2/r_3;
                H(2,2) = x/r;
                H(2,3) = y/r;
            }

        } // spherical  coordinates 
        else if (spherical_coordinates_) {

            const DataType z = state.g_.data_(2);
            // DataType&& zenith = atan2(sqrt(x*x + y*y),z);
            // DataType&& azimuth = atan2(y,x);
            DataType&& tmp = sqrt(x*x+y*y);
            DataType&& r_sqrd = pow(r,2);
            DataType&& r_3 = pow(r,3.0);

            H(0,0) = x/r;
            H(0,1) = y/r;
            H(0,2) = z/r;
            H(1,0) = -y/pow(tmp,2);
            H(1,1) = x/pow(tmp,2);
            H(2,0) = x*z/(r_sqrd*tmp);
            H(2,1) = y*z/(r_sqrd*tmp);
            H(2,2) = -tmp/r_sqrd;
            
            if(has_vel_) {

                DataType tmp2 = (state.g_.data_.transpose()*state.u_.data_.block(0,0,meas_pose_rows_,1))(0,0);

                H(3,0) = state.u_.data_(0,0)/r - tmp2*x/r_3;
                H(3,1) = state.u_.data_(1,0)/r - tmp2*y/r_3;
                H(3,2) = state.u_.data_(2,0)/r - tmp2*z/r_3;
                H(3,3) = x/r;
                H(3,4) = y/r;
                H(3,5) = z/r;
            }

        } else {
            throw std::runtime_error("SourceRadarR2_R3::DerivedGetLinObsMatState: The measurement dimensions are not valid.");
        }
    } else {
        std::cerr << "SourceRadarR2_R3::DerivedGetLinObsMatState: The range of target is less than 1e-4. Setting the Jacobian H to zero. ";
    }





    return H;
}

//-------------------------------------------------------------------------------------------------------------------------

template <typename _State, MeasurementTypes _MeasurementType, template <typename > typename _Transformation>
typename SourceRadarR2_R3<_State, _MeasurementType, _Transformation>::MatV SourceRadarR2_R3<_State, _MeasurementType, _Transformation>::DerivedGetLinObsMatSensorNoise(const State& state) {

    return Base::V_;
}

//---------------------------------------------------------------------------

template <typename _State, MeasurementTypes _MeasurementType, template <typename > typename _Transformation>
typename SourceRadarR2_R3<_State, _MeasurementType, _Transformation>::Measurement SourceRadarR2_R3<_State, _MeasurementType, _Transformation>::DerivedGetEstMeas(const State& state) {
    
    Measurement meas;
    meas.pose.resize(meas_pose_rows_,1);
    meas.twist.resize(meas_twist_rows_,1);
    const DataType&& r = state.g_.data_.norm();
    const DataType x = state.g_.data_(0);
    const DataType y = state.g_.data_(1);
    DataType&& azimuth = atan2(y,x);

    // if (azimuth < 0.0) {
    //     azimuth += 2*M_PI;
    // }
    
    // polar coordinates
    if(polar_coordinates_) {

        
        meas.pose << r, azimuth;
        utilities::WrapAngle(meas.pose(1));
        if(has_vel_) {
            meas.twist << (x*state.u_.data_(0) + y*state.u_.data_(1))/r;
        }

    } // spherical  coordinates 
    else if (spherical_coordinates_) {

        const DataType z = state.g_.data_(2);
        const DataType&& zenith = atan2(sqrt(x*x + y*y),z);
        meas.pose << r, azimuth, zenith;
        utilities::WrapAngle(meas.pose(1));
        utilities::WrapAngle(meas.pose(2));

        if(has_vel_) {
            meas.twist << (x*state.u_.data_(0) + y*state.u_.data_(1) + z*state.u_.data_(2))/r;
        }

    } else {
        throw std::runtime_error("SourceRadarR2_R3::DerivedGetEstMeas: The measurement dimensions are not valid.");
    }
    
    meas.type = measurement_type_;

    return meas;
} 

//---------------------------------------------------------------------------
template <typename _State, MeasurementTypes _MeasurementType, template <typename > typename _Transformation>
typename SourceRadarR2_R3<_State, _MeasurementType, _Transformation>::VecMeas SourceRadarR2_R3<_State,_MeasurementType, _Transformation>::DerivedOMinus(const Measurement& m1, const Measurement& m2) {


    VecMeas error;
    error(0) = m1.pose(0) - m2.pose(0);
    error(1) = m1.pose(1) - m2.pose(1); 
    utilities::WrapAngle(error(1));
    if (spherical_coordinates_) {
        error(2) = m1.pose(2)-m2.pose(2);
        utilities::WrapAngle(error(2));
    }



    if(Base::has_vel_) {
        error.block(meas_pose_dim_,0,1,1) = m1.twist - m2.twist;
    }
    
    return error;
}

//---------------------------------------------------------------------------------------------
template <typename _State, MeasurementTypes _MeasurementType, template <typename > typename _Transformation>
typename SourceRadarR2_R3<_State, _MeasurementType, _Transformation>::Measurement SourceRadarR2_R3<_State,_MeasurementType, _Transformation>::DerivedGenerateRandomMeasurement(const MatMeasCov& meas_std, const State& state) const {
    Measurement m;

    VecMeas deviation = meas_std*utilities::GaussianRandomGenerator(meas_std.rows());

    m = DerivedGetEstMeas(state);
    m.source_index = this->params_.source_index_;

    m.pose += deviation.block(0,0,meas_pose_dim_,1);

    if(has_vel_) {
        m.twist += deviation.block(meas_pose_dim_,0,1,1);
    }

    return m;
}








} // namespace rransac


#endif // RRANSAC_COMMON_SOURCES_SOURCE_R2_R3_RADAR_H_
