#ifndef RRANSAC_COMMON_MODELS_SEN_POS_VEL_H_
#define RRANSAC_COMMON_MODELS_SEN_POS_VEL_H_

#include "common/models/model_base.h"


namespace rransac {

template <typename tState, typename tTransformation>
class ModelSENPosVel : public ModelBase<SourceSENPosVel<tState>, tTransformation,  tState::g_type_::dim_ + tState::u_type_::dim_ - tState::u_type_::dim_t_vel_ + 1, ModelSENPosVel<tState, tTransformation>> {

public:

typedef tState State;
typedef tTransformation Transformation;

static constexpr unsigned int g_dim_ = State::g_type_::dim_;
static constexpr unsigned int cov_dim_ = State::g_type_::dim_ + State::u_type_::dim_ - State::u_type_::dim_t_vel_ + 1;
static constexpr unsigned int l_dim_ =  State::u_type_::dim_a_vel_ + 1;
typedef Eigen::Matrix<double,cov_dim_,cov_dim_> Mat;

/**
 * Computes the Jacobian of the state transition function with respect to the state evaluated at the current state estimate.
 * @param[in] state The state to linearize about
 * @param[in] dt A time interval
 * @return The Jacobian \f$ F_k\f$. 
 */ 
Mat GetLinTransFuncMatState(const State& state, const double dt);

/**
 * Computes the Jacobian of the state transition function with respect to the noise evaluated at the current state estimate.
 * @param[in] state The state to linearize about
 * @param[in] dt  A time interval
 * @return Returns the Jacobian \f$ G_k \f$
 */
Mat GetLinTransFuncMatNoise(const State& state, const double dt);

/**
* Update the state of the model using the provided state_update. The state_update provided is augmented to account
* for the additional translationals velocities (which are zero) before being added to state.
* @param state_update An element of the lie algebra of the state used to update the state. 
*/
void UpdateState(const Eigen::Matrix<double,cov_dim_,1>& state_update);

/**
 * Returns a Random State
 */ 
static State GetRandomState();

};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                            Definitions
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename tState, typename tTransformation>
typename ModelSENPosVel<tState, tTransformation>::Mat  ModelSENPosVel<tState, tTransformation>::GetLinTransFuncMatState(const State& state, const double dt) {  
    Eigen::Matrix<double, g_dim_,g_dim_> tmp = (state.u_*dt).Jr()*dt;  
    this->F_.block(0,0,g_dim_,g_dim_) = typename State::g_type_(State::u_type_::Exp(state.u_.data_*dt)).Adjoint();
    this->F_.block(0,g_dim_, g_dim_, 1) = tmp.block(0,0,g_dim_,1); // Jacobian w.r.t. rho x
    this->F_.block(0,g_dim_+1, g_dim_, State::u_type_::dim_a_vel_) = tmp.block(0, State::u_type_::dim_t_vel_, g_dim_, State::u_type_::dim_a_vel_); // Jacobian w.r.t. angular velocities
    return this->F_;
}

//--------------------------------------------------------------------------------------------------------------------

template <typename tState, typename tTransformation>
typename ModelSENPosVel<tState, tTransformation>::Mat ModelSENPosVel<tState, tTransformation>::GetLinTransFuncMatNoise(const State& state, const double dt){
    Eigen::Matrix<double, g_dim_,g_dim_> tmp = (state.u_*dt).Jr()*dt; 
    this->G_.block(0,0,g_dim_, g_dim_) = tmp;
    this->G_.block(0,g_dim_, g_dim_, 1) = tmp.block(0,0,g_dim_,1)*dt/2.0;
    this->G_.block(0,g_dim_+1, g_dim_, State::u_type_::dim_a_vel_) = tmp.block(0, State::u_type_::dim_t_vel_, g_dim_, State::u_type_::dim_a_vel_)*dt/2.0; // Jacobian w.r.t. angular velocities

    this->G_.block(g_dim_,g_dim_,l_dim_,l_dim_)= Eigen::Matrix<double,l_dim_,l_dim_>::Identity()*dt;
    return this->G_;

}

//--------------------------------------------------------------------------------------------------------------------

template <typename tState, typename tTransformation>
void ModelSENPosVel<tState, tTransformation>::UpdateState(const Eigen::Matrix<double,cov_dim_,1>& state_update) {
    Eigen::Matrix<double,g_dim_,1> twist_update;
    twist_update.setZero();
    twist_update(0,0) = state_update(g_dim_,0);          // get rho_x
    twist_update.block(State::u_type_::dim_t_vel_,0, State::u_type_::dim_a_vel_,1) = state_update.block(g_dim_+1,0,State::u_type_::dim_a_vel_,1);
    this->state_.g_.OPlusEq(state_update.block(0,0,g_dim_,1));
    this->state_.u_.data_ += twist_update;
}

//--------------------------------------------------------------------------------------------------------------------

template <typename tState, typename tTransformation>
tState ModelSENPosVel<tState, tTransformation>::GetRandomState(){
    State state = State::Random();

    state.g_.R_.block(0,0,state.u_.p_.rows(),1) = state.u_.p_.normalized(); 

    if(state.g_.dim_ == 3) {
        state.g_.R_.block(0,1,state.u_.p_.rows(),1) << - state.g_.data_(1,0), state.g_.data_(0,0);
    } else {

        state.g_.R_.block(0,1,state.u_.p_.rows(),1) << -state.g_.data_(1,0), state.g_.data_(0,0), state.g_.data_(2,0);
        state.g_.R_.block(0,2,state.u_.p_.rows(),1) = State::g_type_::rot_algebra::Wedge(state.g_.data_.block(0,0,state.u_.p_.rows(),1))*state.g_.data_.block(0,1,state.u_.p_.rows(),1);

    }

    double px = state.u_.p_.norm();
    state.u_.p_.setZero();
    state.u_.p_(0,0) = px;
    return  state;
}

} // rransac


#endif // RRANSAC_COMMON_MODELS_SEN_POS_VEL_H_