#include "common/models/model_cv_r2.h"

namespace rransac {



void ModelCvR2::Init(const Parameters& params) {

    identity_2d_.setIdentity();

    // Initialize F matrix
    F_.block(0,0,2,2).setIdentity();
    F_.block(2,0,2,2).setZero();
    F_.block(2,2,2,2).setIdentity();

    // Initialize the G matrix
    G_.block(2,0,2,2).setZero();

    // Initialize Q with parameters
    Q_ = params.process_noise_cov_;

    
    
}

//--------------------------------------------------------------------------------------------------------

// We assume a constant velocity model so we must update the pose, but the twist stays the same.
ModelCvR2::State ModelCvR2::PropagateState(const lie_groups::State<G,U>& state, const double dt) {
    ModelCvR2::State s(state);
    s.g_.BoxPlus(s.u_*dt);
    return s;

}

//--------------------------------------------------------------------------------------------------------

ModelCvR2::Mat ModelCvR2::GetLinTransFuncMatState(const ModelCvR2::State& state, const double dt) {
    F_.block(0,2,2,2) = identity_2d_*dt;
    return F_;
}

//--------------------------------------------------------------------------------------------------------

ModelCvR2::Mat ModelCvR2::GetLinTransFuncMatNoise(const ModelCvR2::State& state, const double dt) {

    G_.block(0,0,2,2) = identity_2d_*dt;
    G_.block(2,2,2,2) =  G_.block(0,0,2,2);
    G_.block(0,2,2,2) =  G_.block(0,0,2,2)*dt/2.0;

}

//--------------------------------------------------------------------------------------------------------

void ModelCvR2::PropagateModel(const double dt) {

    state_.g_.BoxPlus(state_.u_*dt);

    Eigen::Matrix4d&& F = GetLinTransFuncMatState(this->state_, dt);
    Eigen::Matrix4d&& G = GetLinTransFuncMatNoise(this->state_, dt);

    err_cov_ = F*err_cov_*F.transpose() + G*Q_*G.transpose();

}

//--------------------------------------------------------------------------------------------------------

void UpdateModel(const Parameters& param) {

}


} // namespace rransac
