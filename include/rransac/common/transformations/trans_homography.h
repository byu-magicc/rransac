#ifndef RRANSAC_COMMON_TRANS_HOMOGRAPHY_R2_H_
#define RRANSAC_COMMON_TRANS_HOMOGRAPHY_R2_H_
#pragma once


#include <Eigen/Core>
#include "rransac/common/measurement/measurement_base.h"
#include "rransac/common/transformations/transformation_base.h"
#include "lie_groups/state.h"
#include "lie_groups/lie_groups/SE2.h"

namespace rransac
{

constexpr int TransformHomographyCovDim(int state_dim) {
    return (state_dim == 2) ? 4: 5;
} 

/** \class TransHomography
 * This transformation is used when target tracking is done on an image plane and the measurement source is a camera. 
 * It transforms the measurements and tracks using the homography and is compatible with SourceSENPosVel and ModelSENPosVel when the
 * target's configuration manifold is SE2 and the measurement space is R2, and it is also compatible with SourceRN and ModelRN when
 * the target's configuration manifold is R2 and the measurement space is R2. 
*/

template<class _State>
class TransformHomography : public TransformBase<TransformDerivedTraits<_State,Eigen::Matrix<typename _State::DataType,3,3>,Eigen::Matrix<typename _State::DataType,TransformHomographyCovDim(_State::Group::dim_),TransformHomographyCovDim(_State::Group::dim_)>,false> ,TransformHomography> {

public:

typedef TransformBase<TransformDerivedTraits<_State,Eigen::Matrix<typename _State::DataType,3,3>,Eigen::Matrix<typename _State::DataType,TransformHomographyCovDim(_State::Group::dim_),TransformHomographyCovDim(_State::Group::dim_)>,false> ,TransformHomography> Base;
typedef typename Base::State State;                                      /**< The State type being used. */
typedef typename Base::DataType DataType;                                /**< The scalar data type. */
typedef typename Base::TransformDataType TransformDataType;              /**< The transform data type being used. It is either an element of SE2 for R2 or SE3 for R3. */
typedef typename Base::MatCov MatCov;                                    /**< The covariance type of the track, and the transform jacobian type. */
typedef typename Base::Measurement Measurement;                          /**< The measurement type. */



typedef Eigen::Matrix<DataType,2,1> Vec2d;  /**< The object type of the measurement. */
typedef Eigen::Matrix<DataType,2,2> Mat2d;  
typedef Eigen::Matrix<DataType,4,4> Mat4d;


// Components of the Homograpy H = [H1, h2; h3_T^T, h4] where T stands for transpose.
Eigen::Matrix<DataType,2,2> H1_;               /**< The homography is represented as a 3x3 matrix and can be segmented as H = [H1, h2; h3_T^T, h4]. */  
Eigen::Matrix<DataType,2,1> h2_;               /**< The homography is represented as a 3x3 matrix and can be segmented as H = [H1, h2; h3_T^T, h4]. */  
Eigen::Matrix<DataType,1,2> h3_T_;             /**< The homography is represented as a 3x3 matrix and can be segmented as H = [H1, h2; h3_T^T, h4]. */  
Eigen::Matrix<DataType,1,1> h4_;               /**< The homography is represented as a 3x3 matrix and can be segmented as H = [H1, h2; h3_T^T, h4]. */  
Eigen::Matrix<DataType,2,2> H1p_;              /**< In the case the homography is being used to transform a track on SE2, the matrix H1_ represents 
                                                    a rotation about the z axis. This matrix could be noisy and  not on the manifold SO2. IN that case,
                                                    we must project it back onto SO2. This projection is the matrix H1p_. */

/**
 * Used to initialize the object, but it doesn't need to initialize anyting.
 */ 
void DerivedInit() {};

/** 
 * Sets the data and the block components of the homography.
 * @param data The data required to transform the measurements, states, and error covariance
 */ 
void DerivedSetData(const TransformDataType data) {
    this->data_ = data;
    H1_ = data.block(0,0,2,2);
    h2_ = data.block(0,2,2,1);
    h3_T_ = data.block(2,0,1,2);
    h4_ = data.block(2,2,1,1);
}

/** 
 * Transforms the measurement from the previous tracking frame to the current one.
 * @param[in] meas The measurement to be transformed.
 */ 
void DerivedTransformMeasurement(Measurement& meas) const {

    if (meas.twist.rows() != 0) {
        Mat2d&& tmp = ConstructTranslationalVelTransform(meas.pose);
        meas.twist = tmp*meas.twist;
    }

    meas.pose = TransformPosition(meas.pose);
}

static Measurement DerivedTransformMeasurement(const Measurement& meas, const TransformDataType& transform_data) {
    throw std::runtime_error("TransformHomography::DerivedTransformMeasurement The method is not implemented.");
}

/** 
 * Transforms the track provided that the state is SE2_se2 or R2_r2.
 * @param[in] state The track's state to be transformed.
 * @param[in] cov   The track's error covariance to be transformed.
 */ 
void DerivedTransformTrack(State& state, MatCov& cov) const {
    throw std::runtime_error("TransformHomography::DerivedTransformTrack The method is not implemented.");
}

/** 
 * Transforms the state and error covariance using user provided transform data.
 * This function is useful when you need to transform a track to a measurement frame. 
 * @param[in] state The track's state to be transformed.
 * @param[in] cov   The track's error covariance to be transformed.
 * @param[in] transform_data The data used to transform the state and error covariance
 */ 
static void DerivedTransformTrack(State& state, MatCov& cov, const TransformDataType& transform_data) {
   throw std::runtime_error("TransformHomography::DerivedTransformTrack This method is not implemented.");
}

/** 
 * Transforms the state using user provided transform data.
 * This function is useful when you need to transform a track to a measurement frame. 
 * @param[in] state The track's state to be transformed.
 * @param[in] transform_data The data used to transform the state and error covariance
 */ 
static State DerivedTransformState(const State& state, const TransformDataType& transform_data) {
    throw std::runtime_error("TransformHomography::DerivedTransformState This method is not implemented.");
}

/** 
 * Returns the Jacobian of the transformation
 * @param[in] state The state of the target after it has been transformed using transform_data
 * @param[in] transform_data The data used in the transformation
 */ 
static MatCov DerivedGetTransformationJacobian(const State& state, const TransformDataType& transform_data) {
   throw std::runtime_error("TransformHomography::DerivedTransformState This method is not implemented.");
}

/**
 * Generates random transform data. The function can use the parameter scalar in order to 
 * generate a larger distribution of random transformations.
 * @param scalar A scalar used to generate a larger distribution. 
 */ 
static TransformDataType DerivedGetRandomTransform(const DataType scalar = static_cast<DataType>(1.0)){
    lie_groups::SE3<DataType> relative_pose(lie_groups::SE3<DataType>::Random(scalar));
    TransformDataType random_transform;
    Eigen::Matrix<DataType,3,1> normal_vec;
    normal_vec << 0,0,1;
    random_transform = relative_pose.R_ + relative_pose.t_*normal_vec.transpose()/ (normal_vec.transpose()*relative_pose.t_);

    return random_transform.normalized();
}

// private: 

/**
 * Transforms the position of the pixel from the previous frame to the current frame
 * @param[in] pos The pixel location
 */ 
Vec2d TransformPosition(const Vec2d& pos) const{
    DataType tmp = (h3_T_*pos + h4_)(0,0);
    return (H1_*pos + h2_)/tmp;
}

/**
 * Construct the linear transformation required to transform the pixel translational velocity and error covariance using 
 * the current pixel position
 * @param[in] pos The pixel location
 */ 
Mat2d ConstructTranslationalVelTransform(const Vec2d& pos) const {
    DataType&& tmp = (h3_T_*pos + h4_)(0,0);
    return (tmp*H1_ - (H1_*pos +h2_)*h3_T_)/(tmp*tmp);
}

/**
 * The covariance transform for R2_r2 is a 4x4 matrix. This function constructs the lower left block of the transform.
 * @param[in] pos The pixel location
 * @param[in] vel The pixel velocity
 */ 
Mat2d ConstructCovTrans12(const Vec2d& pos, const Vec2d& vel) const {
    DataType&& tmp = (h3_T_*pos + h4_)(0,0);
    DataType&& tmp2 = tmp*tmp;
    DataType&& tmp3 = tmp2*tmp;
    return static_cast<DataType>(2.0)*(H1_*pos + h2_)*(h3_T_*vel)*h3_T_/tmp3 - H1_*(vel*h3_T_)/tmp2 - H1_*(h3_T_*vel)/tmp2;
}

/**
 * Constructs the new rotation matrix using the transformed velocity
 * @param vel_transformed The velocity transformed into the current surveillance frame
 */ 
Mat2d TransformRotation(const Vec2d& vel_transformed) const {

    // The velocity is small so set the rotation to identity.
    if (vel_transformed.norm() == 0) {
        return Mat2d::Identity();
    }

    Vec2d tmp = vel_transformed.normalized();
    Mat2d R_transformed;
    R_transformed << tmp(0), -tmp(1), tmp(1), tmp(0);
    return R_transformed;

}


/**
 * Verifies that the transform data provided by the user is in the requested from.
 * @param transform_data The tranformation data to be tested. 
 */
static bool DerivedIsAcceptableTransformData(const TransformDataType& transform_data) {

    bool correct;

    // transform data should be a 3x3 matrix
    if(transform_data.cols() != 3 || transform_data.rows() != 3) {
        correct = false;
    } else {
        correct = true;
    }


    return correct;
} 



};

/////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////   R2_r2                                    /////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

/** 
 * Transforms the track provided that the state is R2_r2.
 * @param[in] state The track's state to be transformed.
 * @param[in] cov   The track's error covariance to be transformed.
 */ 
template<>
inline void TransformHomography<lie_groups::R2_r2>::DerivedTransformTrack(lie_groups::R2_r2& state, Eigen::Matrix<typename lie_groups::R2_r2::DataType,4,4>& cov) const {


    Mat2d&& G = ConstructTranslationalVelTransform(state.g_.data_);
    
    // Transform the covariance
    Mat4d cov_trans;    // matrix used to transform the covariance
    cov_trans.block(0,0,2,2) = G;
    cov_trans.block(2,2,2,2) = G;
    cov_trans.block(0,2,2,2).setZero();
    cov_trans.block(2,0,2,2) = ConstructCovTrans12(state.g_.data_, state.u_.data_);

    cov = cov_trans * cov *cov_trans.transpose();

    // Transform the state
    state.g_.data_ = TransformPosition(state.g_.data_);
    state.u_.data_ = G*state.u_.data_;

}

/////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////   SE2-se2                                  /////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

/** 
 * Sets the data and the block components of the homography. The matrix H1p_ is the projection
 * of the matrix H1_ onto SO2 which represents the rotation about the z axis. 
 * @param data The data required to transform the measurements, states, and error covariance
 */ 
template<>
inline void TransformHomography<lie_groups::SE2_se2>::DerivedSetData(const TransformDataType data) {
    this->data_ = data;
    H1_ = data.block(0,0,2,2);
    h2_ = data.block(0,2,2,1);
    h3_T_ = data.block(2,0,1,2);
    h4_ = data.block(2,2,1,1);

    if ( fabs((H1_.transpose()*H1_).trace() -2.0) == 0) {
        H1p_ = H1_;
    } else {
        Eigen::JacobiSVD<Eigen::Matrix<DataType,2,2>> svd(H1_, Eigen::ComputeFullU | Eigen::ComputeFullV);
        H1p_ = svd.matrixU()*(svd.matrixV().transpose());
    }
}


/** 
 * Transforms the track provided that the state is SE2_se2.
 * @param[in] state The track's state to be transformed.
 * @param[in] cov   The track's error covariance to be transformed.
 */ 
template<>
inline void TransformHomography<lie_groups::SE2_se2>::DerivedTransformTrack(lie_groups::SE2_se2& state, Eigen::Matrix<typename lie_groups::SE2_se2::DataType,5,5>& cov) const {
        
    // Constraining the Homomgraphy to SE2 greatly simplifies the transformation.

    // Transform the position
    state.g_.t_= TransformPosition(state.g_.t_);

    // Transform the rotation
    state.g_.R_ = H1p_*state.g_.R_;

    // Transform the velocity
    state.u_.p_ = state.u_.p_/h4_(0);

    // The angular velocity doesn't change. 


    
    // Transform the covariance
    DataType&& alpha = static_cast<DataType>(1.0)/h4_(0);
    Eigen::Matrix<DataType,5,5> T;
    T.setZero();
    T.diagonal() << alpha, alpha, 1, alpha, 1;
    cov = T*cov*T.transpose();

}

} // namespace rransac

#endif // RRANSAC_COMMON_TRANS_HOMOGRAPHY_R2_H_
