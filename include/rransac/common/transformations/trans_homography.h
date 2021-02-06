#ifndef RRANSAC_COMMON_TRANS_HOMOGRAPHY_R2_H_
#define RRANSAC_COMMON_TRANS_HOMOGRAPHY_R2_H_
#pragma once


#include <Eigen/Core>
#include "common/measurement/measurement_base.h"
#include "common/transformations/transformation_base.h"

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

template<class tState>
class TransformHomography : public TransformBase<Eigen::Matrix<typename tState::DataType,3,3>, tState, Eigen::Matrix<typename tState::DataType,TransformHomographyCovDim(tState::g_type_::dim_),TransformHomographyCovDim(tState::g_type_::dim_)>, TransformHomography<tState>> {

public:
typedef typename tState::DataType DataType; /**< The scalar object for the data. Ex. float, double, etc. */
typedef Eigen::Matrix<DataType,2,1> Vec2d;  /**< The object type of the measurement. */
typedef Eigen::Matrix<DataType,2,2> Mat2d;  
typedef Eigen::Matrix<DataType,3,3> Mat3d;  /**< The data type of the homography. */
typedef Eigen::Matrix<DataType,4,4> Mat4d;
typedef Eigen::Matrix<DataType,TransformHomographyCovDim(tState::g_type_::dim_),TransformHomographyCovDim(tState::g_type_::dim_)> MatCov; /**< The object type of the track's error covariance. */
typedef Mat3d MatData;

// Components of the Homograpy H = [H1, h2; h3_T^T, h4] where T stands for transpose.
Eigen::Matrix<DataType,2,2> H1_;               /**< The homography is represented as a 3x3 matrix and can be segmented as H = [H1, h2; h3_T^T, h4]. */  
Eigen::Matrix<DataType,2,1> h2_;               /**< The homography is represented as a 3x3 matrix and can be segmented as H = [H1, h2; h3_T^T, h4]. */  
Eigen::Matrix<DataType,1,2> h3_T_;             /**< The homography is represented as a 3x3 matrix and can be segmented as H = [H1, h2; h3_T^T, h4]. */  
Eigen::Matrix<DataType,1,1> h4_;               /**< The homography is represented as a 3x3 matrix and can be segmented as H = [H1, h2; h3_T^T, h4]. */  
 

/**
 * Used to initialize the object, but it doesn't need to initialize anyting.
 */ 
void DerivedInit() {};

/** 
 * Sets the data and the block components of the homography.
 * @param data The data required to transform the measurements, states, and error covariance
 */ 
void DerivedSetData(const Mat3d data) {
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
void DerivedTransformMeasurement(Meas<DataType>& meas) const {

    if (meas.twist.rows() != 0) {
        Mat2d&& tmp = ConstructTranslationalVelTransform(meas.pose);
        meas.twist = tmp*meas.twist;
    }

    meas.pose = TransformPosition(meas.pose);
}

/** 
 * Transforms the track provided that the state is SE2_se2 or R2_r2.
 * @param[in] state The track's state to be transformed.
 * @param[in] cov   The track's error covariance to be transformed.
 */ 
void DerivedTransformTrack(tState& state, MatCov& cov) const {
    throw std::runtime_error("TransformHomography::DerivedTransformTrack The state is not supported.");
}

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
    if (vel_transformed.norm() == 0)
        return Mat2d::Identity();

    Vec2d tmp = vel_transformed.normalized();
    Mat2d R_transformed;
    R_transformed << tmp(0), -tmp(1), tmp(1), tmp(0);
    return R_transformed;

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
void TransformHomography<lie_groups::R2_r2>::DerivedTransformTrack(lie_groups::R2_r2& state, Eigen::Matrix<typename lie_groups::R2_r2::DataType,4,4>& cov) const {


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
 * Transforms the track provided that the state is SE2_se2.
 * @param[in] state The track's state to be transformed.
 * @param[in] cov   The track's error covariance to be transformed.
 */ 
template<>
void TransformHomography<lie_groups::SE2_se2>::DerivedTransformTrack(lie_groups::SE2_se2& state, Eigen::Matrix<typename lie_groups::SE2_se2::DataType,5,5>& cov) const {
        
    // Constraining the Homomgraphy to SE2 greatly simplifies the transformation.

    // Transform the position
    state.g_.t_= TransformPosition(state.g_.t_);

    // Transform the rotation
    state.g_.R_ = H1_*state.g_.R_;

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

}

#endif // RRANSAC_COMMON_TRANS_HOMOGRAPHY_R2_H_
