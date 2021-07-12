#ifndef RRANSAC_COMMON_TRANSFORATIONS_SE3_CAM_DEPTH
#define RRANSAC_COMMON_TRANSFORATIONS_SE3_CAM_DEPTH
#pragma once


#include <Eigen/Core>
#include "rransac/common/measurement/measurement_base.h"
#include "rransac/common/transformations/transformation_base.h"
#include "lie_groups/state.h"
#include "lie_groups/lie_groups/SE3.h"

namespace rransac
{


/** \class TransformSE3CamDepth
 * This transformation is used when target tracking is done on SE3 and the target is observed by a camera that measures the depth to the 
 * target and the line of sight vector to the target. It should only be used with the source SourceSE3CamDepth, and the model ModelSENPosVel.
 * If the tracking frame changes, the measurement is not transformed into the current tracking frame since it cannot be. Instead, the transformation
 * from the tracking frame to the measurement frame is stored in the measurement. 
*/

template<typename _State>
class TransformSE3CamDepth : public TransformBase<TransformDerivedTraits<_State,Eigen::Matrix<typename _State::DataType,4,4>,Eigen::Matrix<typename _State::DataType,10,10>,false>, TransformSE3CamDepth> {

public:

typedef TransformBase<TransformDerivedTraits<_State,Eigen::Matrix<typename _State::DataType,4,4>,Eigen::Matrix<typename _State::DataType,10,10>,false>, TransformSE3CamDepth> Base;
typedef typename Base::State State;                                      /**< The State type being used. */
typedef typename Base::DataType DataType;                                /**< The scalar data type. */
typedef typename Base::TransformDataType TransformDataType;              /**< The transform data type being used. It is either an element of SE2 for R2 or SE3 for R3. */
typedef typename Base::MatCov MatCov;                                    /**< The covariance type of the track, and the transform jacobian type. */
typedef typename Base::Measurement Measurement;                          /**< The measurement type. */
static_assert(std::is_same<_State,lie_groups::State< lie_groups::SE3,typename _State::DataType,_State::N>>::value, "TransformSE3CamDepth: The state is not supported");

typedef Eigen::Matrix<DataType,3,1> VecPos;
typedef Eigen::Matrix<DataType,4,1> VecMeas;
/**
 * Used to initialize the object, but it doesn't need to initialize anyting.
 */ 
void DerivedInit() {};

/** 
 * Sets the data
 * @param data The data required to transform the measurements, states, and error covariance. This should be a transformation on SE(3) from the previous tracking
 * frame to the current tracking frame.
 */ 
void DerivedSetData(const TransformDataType data) {
    this->data_ = data;
}

/** 
 * Transforms the measurement from the previous tracking frame to the current one.
 * @param[in,out] meas The measurement to be transformed.
 */ 
void DerivedTransformMeasurement(Measurement& meas) const ;

/** 
 * Transforms the measurement into a different frame. It does this by
 * converting it to cartesian coordinants and then transforming it using the
 * data into another frame. The measurement is left in cartesian coordinates, 
 * and only the pose is transformed. 
 * This method should only be used to compare two measurements in a common frame. 
 * This is used to compare measurements from different radar frames.
 * @param[in] meas The measurement to be transformed.
 */ 
static Measurement DerivedTransformMeasurement(const Measurement& meas, const TransformDataType& transform_data) {
    VecPos pos;
    VecMeas vec_meas;
    Measurement m;
    const Eigen::Matrix<DataType,3,3>& rotation = transform_data.block(0,0,3,3);
    const Eigen::Matrix<DataType,3,1>& position_offset = transform_data.block(0,3,3,1);
    pos =rotation*meas.pose(0,0)*meas.pose.block(1,0,3,1) + position_offset;
    vec_meas(0) = pos.norm();
    vec_meas.block(1,0,3,1) = pos.normalized();
    m.pose = vec_meas;
    return m;
}


/** 
 * Transforms the track provided that the state is SE3_se3.
 * Since the velocity is expressed in the body frame, it doesn't need to be transformed.
 * Since the error covariance is on the error state, it doesn't need to be transformed. 
 * @param[in,out] state The track's state to be transformed.
 * @param[in,out] cov   The track's error covariance.
 */ 
void DerivedTransformTrack(State& state, MatCov& cov) const {
    state.g_.data_ = this->data_*state.g_.data_;
}


/** 
 * Transforms the state and error covariance using user provided transform data.
 * The transform data is an element of SE(3) and transforms the state from the current
 * tracking frame to a measurement frame.
 * @param[in] state The track's state to be transformed.
 * @param[in] cov   The track's error covariance to be transformed.
 * @param[in] transform_data The data used to transform the state and error covariance
 */ 
static void DerivedTransformTrack(State& state, MatCov& cov, const TransformDataType& transform_data) {
   state.g_.data_ = transform_data*state.g_.data_;
}

/** 
 * Transforms the state using user provided transform data.
 * The transform data is an element of SE(3) and transforms the state from the current
 * tracking frame to a measurement frame.
 * @param[in] state The track's state to be transformed.
 * @param[in] transform_data The data used to transform the state and error covariance
 */ 
static State DerivedTransformState(const State& state, const TransformDataType& transform_data) {
    State transformed_state;
    transformed_state.u_.data_= state.u_.data_;
    transformed_state.g_.data_ = transform_data*state.g_.data_;
    return transformed_state;
}

/** 
 * Returns the Jacobian of the transformation
 * The transform data is an element of SE(3) and transforms the state from the current
 * tracking frame to a measurement frame.
 * @param[in] state The state of the target after it has been transformed using transform_data
 * @param[in] transform_data The data used in the transformation
 */ 
static MatCov DerivedGetTransformationJacobian(const State& state, const TransformDataType& transform_data) {
   return MatCov::Identity();
}


/**
 * Verifies that the transform data provided by the user is in the requested from.
 * @param transform_data The tranformation data to be tested. 
 */
static bool DerivedIsAcceptableTransformData(const TransformDataType& transform_data) {

    return lie_groups::SE3<DataType>::isElement(transform_data);
} 

/**
 * Generates random transform data. The function can use the parameter scalar in order to 
 * generate a larger distribution of random transformations.
 * @param scalar A scalar used to generate a larger distribution. 
 */ 
static TransformDataType GetRandomTransform(const DataType scalar = static_cast<DataType>(1.0)){
    return lie_groups::SE3<DataType>::Random(scalar);
}


};

template<class _State>
void TransformSE3CamDepth<_State>::DerivedTransformMeasurement(Measurement& meas) const {

    if(!meas.transform_state) {
        meas.transform_state = true;
        meas.transform_meas = true;
        meas.transform_data_m_t = this->data_;
        meas.transform_data_t_m = this->data_.inverse();
    } else {
        meas.transform_data_m_t = this->data_ * meas.transform_data_m_t;
        meas.transform_data_t_m = meas.transform_data_t_m * this->data_.inverse();
    }

    // DataType d_oa_a = meas.pose(0,0);
    // Eigen::Matrix<DataType,3,1> s_oa_a = meas.pose.block(1,0,3,1);
    // Eigen::Matrix<DataType,3,1> sd_oa_a = meas.twist;
    // Eigen::Matrix<DataType,3,3> I = Eigen::Matrix<DataType,3,3>::Identity();

    // Eigen::Matrix<DataType,3,3> R_ab = this->data_.block(0,0,3,3);
    // Eigen::Matrix<DataType,3,1> t_ab_b = this->data_.block(0,3,3,1);

    // Eigen::Matrix<DataType,3,1> t_oa_a = d_oa_a*s_oa_a;
    // Eigen::Matrix<DataType,3,1> t_ob_b = R_ab*t_oa_a + t_ab_b;
    // Eigen::Matrix<DataType,3,1> v_oa_a = (I*pow(d_oa_a,2) - (t_oa_a*t_oa_a.transpose())).inverse()*pow(d_oa_a,3)*sd_oa_a;

    // std::cout << "I: " << std::endl << I << std::endl;
    // std::cout << "d_oa_a: " << std::endl << d_oa_a << std::endl;
    // std::cout << "t_oa_a: " << std::endl << t_oa_a << std::endl;
    // std::cout << "sd_oa_a:  " << std::endl << sd_oa_a << std::endl;


    // DataType d_ob_b = sqrt(  (t_ob_b.transpose()*t_ob_b)(0,0)   );
    // Eigen::Matrix<DataType,3,1> s_ob_b = t_ob_b/d_ob_b;
    // Eigen::Matrix<DataType,3,1> sd_ob_b = R_ab*v_oa_a/d_ob_b - (t_ob_b*t_ob_b.transpose()*R_ab*v_oa_a)/pow(d_ob_b,3);

    // meas.pose(0,0) = d_ob_b;
    // meas.pose.block(1,0,3,1) = s_ob_b;
    // meas.twist = sd_ob_b;

    // std::cout << "v: " << std::endl << v_oa_a << std::endl;


}





} // namespace rransac
#endif // RRANSAC_COMMON_TRANSFORATIONS_SE3_CAM_DEPTH