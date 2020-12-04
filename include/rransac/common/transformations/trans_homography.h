#ifndef RRANSAC_COMMON_TRANS_HOMOGRAPHY_R2_H_
#define RRANSAC_COMMON_TRANS_HOMOGRAPHY_R2_H_

#include <Eigen/Core>
#include "common/measurement/measurement_base.h"
#include "common/transformations/transformation_base.h"

namespace rransac
{
/** \class TransHomography
 * Transforms the measurements, states and error covariance using the homography. The homography is Eigen::Matrix3d
 * the state is R2 or SE2 and the measurement is of type MeasurementType::RN_POS, MeasurementType::RN_POS_VEL, MeasurementType::RN_POS.
 * or MeasurementType::SEN_POS_VEL and is of TWO dimensions. When the state is SE2, it is assumed that tracking is done in the virtual image frame; 
 * otherwise it won't work. 
*/

template<class tState, class tMatCov>
class TransformHomography : public TransformBase<Eigen::Matrix3d, tState, tMatCov, TransformHomography<tState,tMatCov>> {

public:

// Components of the Homograpy H = [H1, h2; h3_T^T, h4] where T stands for transpose.
Eigen::Matrix2d H1_;
Eigen::Matrix<double,2,1> h2_;
Eigen::Matrix<double,1,2> h3_T_;
Eigen::Matrix<double,1,1> h4_;

/** 
 * The parent class sets the transformation data member variable. This derived class has the 
 * opportunity to perform other calculations using the data. 
 * @param data The data required to transform the measurements, states, and error covariance
 */ 
void DerivedSetData(Eigen::Matrix3d data) {
    this->data_ = data;
    H1_ = data.block(0,0,2,2);
    h2_ = data.block(0,2,2,1);
    h3_T_ = data.block(2,0,1,2);
    h4_ = data.block(2,2,1,1);
}

/** 
 * Transforms the measurement using data_ from the previous surveillance frame to the current one.
 * @param meas The measurement to be transformed.
 */ 
void DerivedTransformMeasurement(Meas& meas) const {

    if (meas.type == MeasurementTypes::RN_POS_VEL) {
        Eigen::Matrix2d&& tmp = ConstructTranslationalVelTransform(meas.pose);
        meas.twist = tmp*meas.twist;
    }

    meas.pose = TransformPosition(meas.pose);
}

/** 
 * Transforms the track using the transform data. i.e. transform the estimated 
 * state and error covariance.
 * @param state The track's state to be transformed.
 * @param cov   The track's error covariance to be transformed.
 */ 
void DerivedTransformTrack(tState& state, tMatCov& cov) const {
// void TransformTrack(S& state, Eigen::MatrixXd& cov) {
    throw std::runtime_error("TransformHomography: The state is not supported.");
}

/**
 * Transforms the position of the pixel from the previous frame to the current frame
 * @param pos The pixel location
 */ 
Eigen::Matrix<double,2,1> TransformPosition(const Eigen::Matrix<double,2,1>& pos) const{
    double tmp = (h3_T_*pos + h4_)(0,0);
    return (H1_*pos + h2_)/tmp;
}

/**
 * Construct the linear transformation required to transform the pixel translational velocity and error covariance using 
 * the current pixel position
 * @param pos The pixel location
 */ 
Eigen::Matrix2d ConstructTranslationalVelTransform(const Eigen::Matrix<double,2,1>& pos) const {
    double&& tmp = (h3_T_*pos + h4_)(0,0);
    return (tmp*H1_ - (H1_*pos +h2_)*h3_T_)/(tmp*tmp);
}

/**
 * The covariance transform for R2_r2 is a 4x4 matrix. This function constructs the lower left block of the transform.
 * @param pos The pixel location
 * @param vel The pixel velocity
 */ 
Eigen::Matrix2d ConstructCovTrans12(const Eigen::Matrix<double,2,1>& pos, const Eigen::Matrix<double,2,1>& vel) const {
    double&& tmp = (h3_T_*pos + h4_)(0,0);
    double&& tmp2 = tmp*tmp;
    double&& tmp3 = tmp2*tmp;
    return 2.0*(H1_*pos + h2_)*(h3_T_*vel)*h3_T_/tmp3 - H1_*(vel*h3_T_)/tmp2 - H1_*(h3_T_*vel)/tmp2;
}

/**
 * Constructs the new rotation matrix using the transformed velocity
 * @param vel_transformed The velocity transformed into the current surveillance frame
 */ 
Eigen::Matrix2d TransformRotation(const Eigen::Matrix<double,2,1>& vel_transformed) const {

    // The velocity is small so set the rotation to identity.
    if (vel_transformed.norm() == 0)
        return Eigen::Matrix2d::Identity();

    Eigen::Matrix<double,2,1> tmp = vel_transformed.normalized();
    Eigen::Matrix2d R_transformed;
    R_transformed << tmp(0), -tmp(1), tmp(1), tmp(0);
    return R_transformed;

}

};

/////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////   R2_r2                                    /////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

template<>
void TransformHomography<lie_groups::R2_r2, Eigen::Matrix<double,4,4>>::DerivedTransformTrack(lie_groups::R2_r2& state, Eigen::Matrix<double,4,4>& cov) const {


    Eigen::Matrix2d&& G = ConstructTranslationalVelTransform(state.g_.data_);
    
    // Transform the covariance
    Eigen::Matrix4d cov_trans;    // matrix used to transform the covariance
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

template<>
void TransformHomography<lie_groups::SE2_se2, Eigen::Matrix<double,5,5>>::DerivedTransformTrack(lie_groups::SE2_se2& state, Eigen::Matrix<double,5,5>& cov) const {
        
    // Constraining the Homomgraphy to SE2 greatly simplifies the transformation.

    // Transform the position
    state.g_.t_= TransformPosition(state.g_.t_);

    // Transform the rotation
    state.g_.R_ = H1_*state.g_.R_;

    // Transform the velocity
    state.u_.p_ = state.u_.p_/h4_(0);

    // The angular velocity doesn't change. 


    
    // Transform the covariance
    double&& alpha = 1/h4_(0);
    Eigen::Matrix<double,5,5> T;
    T.setZero();
    T.diagonal() << alpha, alpha, 1, alpha, 1;
    cov = T*cov*T.transpose();

}

}

#endif // RRANSAC_COMMON_TRANS_HOMOGRAPHY_R2_H_
