#ifndef RRANSAC_COMMON_TRANSFORMATIONS_TRANS_RADAR_SE2_SE3_WITH_SE2_SE3_H_
#define RRANSAC_COMMON_TRANSFORMATIONS_TRANS_RADAR_SE2_SE3_WITH_SE2_SE3_H_
#pragma once

#include <Eigen/Core>

#include "lie_groups/state.h"
#include "rransac/common/measurement/measurement_base.h"
#include "rransac/common/transformations/transformation_base.h"
#include "lie_groups/lie_groups/SE2.h"
#include "lie_groups/lie_groups/SE3.h"
#include "rransac/common/utilities.h"

namespace rransac
{

template<typename _State, typename TransformDataType>
struct IsSE2orSE3 {
    static bool Value(const TransformDataType& transform_data);
};

//---------------------------------------------------------------------------------------------------------------------------

template<typename _DataType, int _NumTangentSpaces, typename TransformDataType>
struct IsSE2orSE3<lie_groups::State<lie_groups::SE2,_DataType,3,_NumTangentSpaces>,TransformDataType> {
    static bool Value(const TransformDataType& transform_data){ return lie_groups::SE2<_DataType>::isElement(transform_data);}
};

//---------------------------------------------------------------------------------------------------------------------------

template<typename _DataType, int _NumTangentSpaces, typename TransformDataType>
struct IsSE2orSE3<lie_groups::State<lie_groups::SE3,_DataType,6,_NumTangentSpaces>,TransformDataType> {
    static bool Value(const TransformDataType& transform_data){ return lie_groups::SE3<_DataType>::isElement(transform_data);}
};



//---------------------------------------------------------------------------------------------------------------------------


template<typename _State>
class TransRadarSE2SE3WithSE2SE3 : public TransformBase<TransformDerivedTraits<_State,Eigen::Matrix<typename _State::DataType,_State::Group::size1_,_State::Group::size2_>,Eigen::Matrix<typename _State::DataType,_State::Group::dim_+_State::Group::dim_rot_+1,_State::Group::dim_+_State::Group::dim_rot_+1>,false>, TransRadarSE2SE3WithSE2SE3> {

public:

typedef TransformBase<TransformDerivedTraits<_State,Eigen::Matrix<typename _State::DataType,_State::Group::size1_,_State::Group::size2_>,Eigen::Matrix<typename _State::DataType,_State::Group::dim_+_State::Group::dim_rot_+1,_State::Group::dim_+_State::Group::dim_rot_+1>,false>, TransRadarSE2SE3WithSE2SE3> Base;
typedef typename Base::State State;                                      /**< The State type being used. */
typedef typename Base::DataType DataType;                                /**< The scalar data type. */
typedef typename Base::TransformDataType TransformDataType;              /**< The transform data type being used. It is either an element of SE2 for R2 or SE3 for R3. */
typedef typename Base::MatCov MatCov;                                    /**< The covariance type of the track, and the transform jacobian type. */
typedef typename Base::Measurement Measurement;                          /**< The measurement type. */


static constexpr unsigned int group_dim_ = State::Group::dim_;                  /**< The dimension of the group. */
static constexpr unsigned int num_tangent_spaces_ = State::NumTangentSpaces;    /**< The number of tangent spaces in the group. Ex 1 tangent space is velocity and 2 is velocity and acceleration. */
static constexpr unsigned int state_dim_ = State::dim_;                         /**< The dimension of the state. */
static constexpr bool polar_coordinates_ = group_dim_ == 3 ? true : false;      /**< This variable is true if the measurement is in polar coordinates. */
static constexpr bool spherical_coordinates_ = group_dim_ == 6 ? true : false;  /**< This value is true if the measurement is in spherical coordinates. */
static constexpr unsigned int meas_dim_ = polar_coordinates_ ? 2 : 3; /**< The dimension of the measurement. */
/**
 * Used to initialize the transformation object.
 */ 
void DerivedInit() {}

/** 
 * Sets the transformation data member variable. This is data is used to transform
 * the tracks and measurements to the current tracking frame. The tracking frame is assumed to be 
 * stationary and so this function does nothing.
 * @param[in] data The data required to transform the measurements, states, and error covariance
 */ 
void DerivedSetData(const TransformDataType& data) {
    throw std::runtime_error("TransRadarR2R3WithSE2SE3::DerivedSetData: The tracking frame shouldn't be moving, so this function should not be called. ");
}


/** 
 * Transforms the measurement from the previous tracking frame to the current one.
 * Since the tracking frame is assumed fixed, this function shouldn't be called. 
 * @param[in] meas The measurement to be transformed.
 */ 
void DerivedTransformMeasurement(Measurement& meas) const {
    throw std::runtime_error("TransRadarR2R3WithSE2SE3::DerivedTransformMeasurement: The tracking frame shouldn't be moving, so this function should not be called. ");
}

/** 
 * Transforms the measurement into a different frame. It does this by
 * converting it to cartesian coordinants and then transforming it using the
 * data into another frame. The measurement is left in cartesian coordinates, 
 * and only the pose is transformed. 
 * This method should only be used to compare two measurements in a common frame. 
 * This is used to compare measurements from different radar frames.
 * @param[in] meas The measurement to be transformed.
 */ 
static Measurement DerivedTransformMeasurement(const Measurement& meas, const TransformDataType& transform_data);

/** 
 * Transforms the track using the transform data from the previous tracking frame
 * to the current one. Since the tracking frame is stationary, this function should 
 * not be called. 
 * @param[in] state The track's state to be transformed.
 * @param[in] cov   The track's error covariance to be transformed.
 */ 
void DerivedTransformTrack(State& state, MatCov& cov) const {
        throw std::runtime_error("TransRadarR2R3WithSE2SE3::DerivedTransformTrack: The tracking frame shouldn't be moving, so this function should not be called. ");

}

/** 
 * Transforms the state and error covariance using user provided transform data.
 * This function is useful when you need to transform a track to a measurement frame. 
 * @param[in] state The track's state to be transformed.
 * @param[in] cov   The track's error covariance to be transformed.
 * @param[in] transform_data The data used to transform the state and error covariance.
 * This is an element of SE2 if the state is R2, or SE3 if the state is R3;
 */ 
static void DerivedTransformTrack(State& state, MatCov& cov, const TransformDataType& transform_data);

/** 
 * Transforms the state using user provided transform data.
 * This function is useful when you need to transform a track to a measurement frame. 
 * @param[in] state The track's state to be transformed.
 * @param[in] transform_data The data used to transform the state and error covariance.
 * This is an element of SE2 if the state is R2, or SE3 if the state is R3;
 */ 
static State DerivedTransformState(const State& state, const TransformDataType& transform_data);

/** 
 * Returns the Jacobian of the transformation
 * @param[in] state The state of the target after it has been transformed using transform_data
 * @param[in] transform_data The data used in the transformation
 */ 
static MatCov DerivedGetTransformationJacobian(const State& state, const TransformDataType& transform_data);
/**
 * Verifies that the transform data provided by the user is in the requested from.
 * @param transform_data The data used to transform the state and error covariance.
 * This is an element of SE2 if the state is R2, or SE3 if the state is R3;
 */
static bool DerivedIsAcceptableTransformData(const TransformDataType& transform_data){return IsSE2orSE3<State,TransformDataType>::Value(transform_data);}

/**
 * Generates random transform data. The function can use the parameter scalar in order to 
 * generate a larger distribution of random transformations.
 * @param scalar A scalar used to generate a larger distribution. 
 */ 
static TransformDataType DerivedGetRandomTransform(const DataType scalar = static_cast<DataType>(1.0)){

    Eigen::Matrix<DataType,Eigen::Dynamic,Eigen::Dynamic> random_transform;

    if (polar_coordinates_) {
        random_transform = lie_groups::SE2<DataType>::Random(scalar);
    } else {
        random_transform = lie_groups::SE3<DataType>::Random(scalar);
    }

    return random_transform;
}

};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Definitions
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<class _State>
typename TransRadarSE2SE3WithSE2SE3<_State>::Measurement TransRadarSE2SE3WithSE2SE3<_State>::DerivedTransformMeasurement(const Measurement& meas, const TransformDataType& transform_data) {

    static constexpr size_t rot_dim = meas_dim_;
    static constexpr size_t pos_dim = meas_dim_;
    Eigen::Matrix<DataType,pos_dim,1> position;
    Measurement m;
    m.pose.resize(meas_dim_,1);
    DataType r = meas.pose(0);
    DataType azimuth = meas.pose(1);
    const Eigen::Matrix<DataType,rot_dim,rot_dim>& R = transform_data.block(0,0,rot_dim,rot_dim);
    const Eigen::Matrix<DataType,pos_dim,1>& position_offset = transform_data.block(0,rot_dim,pos_dim,1); 

    if(polar_coordinates_) {
        position << r*cos(azimuth), r*sin(azimuth);
        position = R*position+position_offset;
        m.pose(0) = position.norm();
        m.pose(1) = atan2(position(1),position(0));
        utilities::WrapAngle(m.pose(1));

    } else if(spherical_coordinates_) {
        DataType zenith = meas.pose(2);
        position << r*cos(azimuth)*sin(zenith),r*sin(azimuth)*sin(zenith),r*cos(zenith);
        position = R*position+position_offset;
        m.pose(0) = position.norm();
        m.pose(1) = atan2(position(1),position(0));
        m.pose(2) = atan2( position.block(0,0,2,1).norm(), position(2));
        utilities::WrapAngle(m.pose(1));
        utilities::WrapAngle(m.pose(2));

    } else {
        throw std::runtime_error("TransRadarSE2SE3WithSE2SE3::DerivedTransformMeasurement: Measurement should be in spherical or polar coordinates. Is the state right R2 or R3? ");
    }


    return m;
}

//---------------------------------------------------------------------------------------------------------------------------

template<class _State>
void TransRadarSE2SE3WithSE2SE3<_State>::DerivedTransformTrack(State& state, MatCov& cov, const TransformDataType& transform_data) {
        MatCov TJ = DerivedGetTransformationJacobian(state,transform_data);
        state = DerivedTransformState(state,transform_data);
        cov = TJ*cov*TJ.transpose();
}

//---------------------------------------------------------------------------------------------------------------------------

template<class _State>
typename TransRadarSE2SE3WithSE2SE3<_State>::State TransRadarSE2SE3WithSE2SE3<_State>::DerivedTransformState(const State& state, const TransformDataType& transform_data) {

    State transformed_state;
    transformed_state.g_.data_ = transform_data*state.g_.data_;
    transformed_state.u_.data_ = state.u_.data_;


    return transformed_state;
}

//---------------------------------------------------------------------------------------------------------------------------

template<class _State>
typename TransRadarSE2SE3WithSE2SE3<_State>::MatCov TransRadarSE2SE3WithSE2SE3<_State>::DerivedGetTransformationJacobian(const State& state, const TransformDataType& transform_data) {
    MatCov TJ = MatCov::Identity();
    return TJ;
}

    
} // namespace rransac


#endif // RRANSAC_COMMON_TRANSFORMATIONS_TRANS_RADAR_SE2_SE3_WITH_SE2_SE3_H_