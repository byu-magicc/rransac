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
 * or MeasurementType::SEN_POS_VEL and is of TWO dimensions.
*/

template<class S>
class TransformHomography : public TransformBase<Eigen::Matrix3d, S, TransformHomography<S>> {

public:


/** 
 * The parent class sets the transformation data member variable. This derived class has the 
 * opportunity to perform other calculations using the data. 
 * @param data The data required to transform the measurements, states, and error covariance
 */ 
void SetData(Eigen::Matrix3d data) {
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
void TransformMeasurement(Meas& meas) {
    throw std::runtime_error("TransformHomography: The state is not supported.");
}

/** 
 * Transforms the track using the transform data. i.e. transform the estimated 
 * state and error covariance.
 * @param state The track's state to be transformed.
 * @param cov   The track's error covariance to be transformed.
 */ 
// void TransformTrack(S& state, Eigen::Matrix<double,S::dim_,S::dim_> cov) {
void TransformTrack(S& state, Eigen::MatrixXd& cov) {
    throw std::runtime_error("TransformHomography: The state is not supported.");
}

// protected:

// Components of the Homograpy H = [H1, h2; h3_T^T, h4] where T stands for transpose.
Eigen::Matrix2d H1_;
Eigen::Matrix<double,2,1> h2_;
Eigen::Matrix<double,1,2> h3_T_;
Eigen::Matrix<double,1,1> h4_;

// private:

/**
 * Transforms the position of the pixel from the previous frame to the current frame
 * @param pos The pixel location
 */ 
Eigen::Matrix<double,2,1> TransformPosition(const Eigen::Matrix<double,2,1>& pos){
    double tmp = (h3_T_*pos + h4_)(0,0);
    return (H1_*pos + h2_)/tmp;
}

/**
 * Construct the linear transformation required to transform the pixel velocity and error covariance using 
 * the current pixel position
 * @param pos The pixel location
 */ 
Eigen::Matrix2d ConstructVelTransform(const Eigen::Matrix<double,2,1>& pos){
    double&& tmp = (h3_T_*pos + h4_)(0,0);
    return (tmp*H1_ - (H1_*pos +h2_)*h3_T_)/(tmp*tmp);
    // return H1_;
}

/**
 * The covariance transform for R2_r2 is a 4x4 matrix. This function constructs the lower left block of the transform.
 * @param pos The pixel location
 * @param vel The pixel velocity
 */ 
Eigen::Matrix2d ConstructCovTrans12(const Eigen::Matrix<double,2,1>& pos, const Eigen::Matrix<double,2,1>& vel) {
    double&& tmp = (h3_T_*pos + h4_)(0,0);
    double&& tmp2 = tmp*tmp;
    double&& tmp3 = tmp2*tmp;
    return 2.0*(H1_*pos + h2_)*(h3_T_*vel)*h3_T_/tmp3 - H1_*(vel*h3_T_)/tmp2 - H1_*(h3_T_*vel)/tmp2;
}

};

//-------------------------------------------------------------------------------------------

template<>
void TransformHomography<lie_groups::R2_r2>::TransformMeasurement(Meas& meas) {
    
    if (meas.type == MeasurementTypes::RN_POS_VEL) {
        Eigen::Matrix2d&& tmp = ConstructVelTransform(meas.pose);
        meas.twist = tmp*meas.twist;
    }

    meas.pose = TransformPosition(meas.pose);

}

//-------------------------------------------------------------------------------------------

template<>
// void TransformHomography<lie_groups::R2_r2>::TransformTrack(lie_groups::R2_r2& state, Eigen::Matrix<double,lie_groups::R2_r2::dim_,lie_groups::R2_r2::dim_>& cov) {
// void TransformHomography<lie_groups::R2_r2>::TransformTrack(lie_groups::R2_r2& state, Eigen::MatrixXd cov)
void TransformHomography<lie_groups::R2_r2>::TransformTrack(lie_groups::R2_r2& state, Eigen::MatrixXd& cov) 
{
    
    Eigen::Matrix2d&& G = ConstructVelTransform(state.g_.data_);
    
    // Transform the state
    state.g_.data_ = TransformPosition(state.g_.data_);
    state.u_.data_ = G*state.u_.data_;

    // Transform the covariance
    Eigen::Matrix4d cov_trans;    // matrix used to transform the covariance
    cov_trans.block(0,0,2,2) = G;
    cov_trans.block(2,2,2,2) = G;
    cov_trans.block(0,2,2,2).setZero();
    cov_trans.block(2,0,2,2) = ConstructCovTrans12(state.g_.data_, state.u_.data_);

    cov = cov_trans * cov *cov_trans.transpose();

}

// template<>
// void TransformHomography<lie_groups::R2_r2>::TransformState(lie_groups::R2_r2& state){

// }

// template<>
// void TransformHomography<lie_gorups::SE2_se2>::TransformState(lie_gorups::SE2_se2& state){
    
// }

}

#endif // RRANSAC_COMMON_TRANS_HOMOGRAPHY_R2_H_
