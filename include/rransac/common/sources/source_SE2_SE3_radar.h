#ifndef RRANSAC_COMMON_SOURCES_SOURCE_SE2_SE3_RADAR_H_
#define RRANSAC_COMMON_SOURCES_SOURCE_SE2_SE3_RADAR_H_
#pragma once

#include <typeinfo>

#include "rransac/common/sources/source_base.h"
#include "rransac/common/utilities.h"
#include "lie_groups/utilities.h"
#include "lie_groups/state.h"
#include "lie_groups/lie_algebras/se2.h"

namespace rransac
{
using namespace utilities;

/**
 * \class SourceRadarSE2SE3
 * This source is meant to be used for a target that evolves on SE2 or SE3 and is observed by 
 * a radar sensor with measurements being in polar coordinates for SE2 or spherical coordinates for SE3.
 * For polar coordinates, the measurement pose is given in a 2x1 vector in the order: range, azimuth.
 * For spherical coordinates, the measurement pose is given in a 3x1 vector in the order: range, azimuth, zenith.
 * \latexinclude radar_r2_r3.tex
 *
 */ 

template <typename _State, MeasurementTypes _MeasurementType, template <typename > typename _Transformation>
class SourceRadarSE2SE3 : public SourceBase<SourceDerivedTraits<_State, MatXT<typename _State::DataType>,MatXT<typename _State::DataType>,_Transformation,_State::Group::dim_pos_,1, 0, 0, _State::Group::dim_pos_, 0, MeasHasVelocity<_MeasurementType>::value, _MeasurementType, utilities::CompatibleWithModelSENPosVel>, SourceRadarSE2SE3> {

public:

typedef SourceBase<SourceDerivedTraits<_State, MatXT<typename _State::DataType>,MatXT<typename _State::DataType>,_Transformation,_State::Group::dim_pos_,1, 0, 0, _State::Group::dim_pos_, 0, MeasHasVelocity<_MeasurementType>::value, _MeasurementType, utilities::CompatibleWithModelSENPosVel>, SourceRadarSE2SE3> Base;

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

typedef Eigen::Matrix<DataType,meas_pose_rows_,meas_pose_cols_> MeasPoseDataType; /**< The data type of the measurement pose. */

static constexpr unsigned int state_dim_ = _State::Group::dim_ + _State::Algebra::dim_ - _State::Algebra::dim_t_vel_ + 1; /**< The dimension of the state. */
static constexpr bool polar_coordinates_ = meas_pose_dim_ == 2 ? true : false;      /**< This variable is true if the measurement is in polar coordinates. */
static constexpr bool spherical_coordinates_ = meas_pose_dim_ == 3 ? true : false;  /**< This value is true if the measurement is in spherical coordinates. */

// Perform compatibility checks
static_assert(lie_groups::utilities::StateIsSEN_seN<_State>::value, "SourceRadarSE2SE3: The state is not compatible with the model");
static_assert(_State::Group::dim_ != 3 || _State::Group::dim_ != 6, "SourceRadarSE2SE3: The dimensions of the group element of the state must be either 2 or 3.");
static_assert( measurement_type_ == MeasurementTypes::SE2_SE3_RADAR, "SourceRadarSE2SE3: The measurement type is not compatible with the source."    );


/** 
 * Initializes the Jacobians
 */
SourceRadarSE2SE3();

/** 
 * Initializes the measurement source. Currently does nothing. 
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
SourceRadarSE2SE3<_State,_MeasurementType,_Transformation>::SourceRadarSE2SE3() {

    Base::H_ = Eigen::Matrix<DataType,total_meas_dim_,state_dim_>::Zero();
    Base::V_ = MatMeasCov::Identity();

}

//-----------------------------------------------------------------
template <typename _State, MeasurementTypes _MeasurementType, template <typename > typename _Transformation>
typename SourceRadarSE2SE3<_State,_MeasurementType,_Transformation>::MatH SourceRadarSE2SE3<_State,_MeasurementType,_Transformation>::DerivedGetLinObsMatState(const State& state)  {

    MatH H = Base::H_;

    DataType d = state.g_.t_.norm();


    Eigen::Matrix<DataType,1,State::Group::dim_pos_> PT_R = state.g_.t_.transpose()*state.g_.R_;

    if (polar_coordinates_) {
        
        Eigen::Matrix<DataType,1,2> r1, r2;
        r1 = state.g_.R_.block(0,0,1,2);
        r2 = state.g_.R_.block(1,0,1,2);
        DataType Px = state.g_.t_(0);
        DataType Py = state.g_.t_(1);
        DataType tmp = Px*Px + Py*Py;


        H.block(0,0,1,State::Group::dim_pos_) = PT_R/d;
        H.block(1,0,1,State::Group::dim_pos_) = (Px*r2 - Py*r1)/tmp;
        
    } else {

        Eigen::Matrix<DataType,1,3> r1, r2, r3;
        r1 = state.g_.R_.block(0,0,1,3);
        r2 = state.g_.R_.block(1,0,1,3);
        r3 = state.g_.R_.block(2,0,1,3);

        DataType Px = state.g_.t_(0);
        DataType Py = state.g_.t_(1);
        DataType Pz = state.g_.t_(2);
        DataType tmp = Px*Px + Py*Py;

        H.block(0,0,1,State::Group::dim_pos_) = PT_R/d;
        H.block(1,0,1,State::Group::dim_pos_) = (Px*r2 - Py*r1)/tmp;
        H.block(2,0,1,State::Group::dim_pos_) = (Pz*(Px*r1 + Py*r2)/sqrt(tmp) - sqrt(tmp)*r3)/(d*d);


    }


    return H;

}

//-----------------------------------------------------------------

template <typename _State, MeasurementTypes _MeasurementType, template <typename > typename _Transformation>
typename SourceRadarSE2SE3<_State,_MeasurementType,_Transformation>::MatV SourceRadarSE2SE3<_State,_MeasurementType,_Transformation>::DerivedGetLinObsMatSensorNoise(const State& state) {

  return Base::V_;

}


//-----------------------------------------------------------------
template <typename _State, MeasurementTypes _MeasurementType, template <typename > typename _Transformation>
typename SourceRadarSE2SE3<_State,_MeasurementType,_Transformation>::Measurement SourceRadarSE2SE3<_State,_MeasurementType,_Transformation>::DerivedGetEstMeas(const State& state) {
    Measurement m;
    m.pose = MeasPoseDataType::Zero();
    m.type = measurement_type_;

    m.pose(0) = state.g_.t_.norm();
    DataType Px = state.g_.t_(0);
    DataType Py = state.g_.t_(1);
    m.pose(1) = atan2(Py,Px);
    utilities::WrapAngle(m.pose(1));


    if(spherical_coordinates_) {
        DataType Pz = state.g_.t_(2);

        m.pose(2) = atan2( sqrt(Px*Px + Py*Py), Pz);
        utilities::WrapAngle(m.pose(2));


    } 



    return m;
} 

//-----------------------------------------------------------------
template <typename _State, MeasurementTypes _MeasurementType, template <typename > typename _Transformation>
typename SourceRadarSE2SE3<_State,_MeasurementType,_Transformation>::VecMeas SourceRadarSE2SE3<_State,_MeasurementType,_Transformation>::DerivedOMinus(const Measurement& m1, const Measurement& m2) {

    VecMeas error;

    error(0) = m1.pose(0) - m2.pose(0);
    error(1) = m1.pose(1) - m2.pose(1);
    utilities::WrapAngle(error(1));

    if(spherical_coordinates_) {
        error(2) = m1.pose(2) - m2.pose(2);
        utilities::WrapAngle(error(2));
    }



    return error;

}

//----------------------------------------------------------------------------------------
template <typename _State, MeasurementTypes _MeasurementType, template <typename > typename _Transformation>
typename SourceRadarSE2SE3<_State,_MeasurementType,_Transformation>::Measurement SourceRadarSE2SE3<_State,_MeasurementType,_Transformation>::DerivedGenerateRandomMeasurement(const MatMeasCov& meas_std, const State& state) const {
    Measurement m;
    m.source_index = this->params_.source_index_;
    m.type = measurement_type_;

    VecMeas deviation = meas_std*utilities::GaussianRandomGenerator(meas_std.rows());

    m = DerivedGetEstMeas(state);
    m.pose += deviation.block(0,0,meas_pose_dim_,1); 

    utilities::WrapAngle(m.pose(1));


    if(spherical_coordinates_) {
        utilities::WrapAngle(m.pose(2));
    } 


    return m;
}





} // namespace rransac




#endif // RRANSAC_COMMON_SOURCES_SOURCE_SE2_SE3_RADAR_H_