#ifndef RRANSAC_COMMON_MODELS_SEN_POSE_TWIST_H_
#define RRANSAC_COMMON_MODELS_SEN_POSE_TWIST_H_

#include "common/models/model_base.h"


namespace rransac {

template <typename tState, typename tTransformation>
class ModelSENPoseTwist : public ModelBase< SourceSENPoseTwist<tState>, tTransformation, tState::g_type_::dim_*2, ModelSENPoseTwist<tState, tTransformation>> {

public:

typedef tState State;
typedef tTransformation Transformation;
static constexpr unsigned int cov_dim_ = State::g_type_::dim_*2;
static constexpr unsigned int g_dim_ = State::g_type_::dim_;
typedef Eigen::Matrix<double,2*g_dim_,2*g_dim_> Mat;

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
* Update the state of the model using the provided state_update
* @param state_update An element of the lie algebra of the state used to update the state. 
*/
void UpdateState(const Eigen::Matrix<double,2*g_dim_,1>& state_update);

/**
 * Returns a Random State
 */ 
static State GetRandomState(){ return State::Random(); }

};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                            Definitions
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename tState, typename tTransformation>
typename ModelSENPoseTwist<tState,tTransformation>::Mat ModelSENPoseTwist<tState,tTransformation>::GetLinTransFuncMatState(const State& state, const double dt) {    
    this->F_.block(0,0,g_dim_, g_dim_) = typename State::g_type_(State::u_type_::Exp(state.u_.data_*dt)).Adjoint();
    this->F_.block(0,g_dim_,g_dim_,g_dim_) = (state.u_*dt).Jr()*dt;
    return this->F_;
}

//--------------------------------------------------------------------------------------------------------------------

template <typename tState, typename tTransformation>
typename ModelSENPoseTwist<tState,tTransformation>::Mat ModelSENPoseTwist<tState,tTransformation>::GetLinTransFuncMatNoise(const State& state, const double dt){
    Eigen::Matrix<double, g_dim_, g_dim_> tmp = (state.u_*dt).Jr()*dt;
    this->G_.block(0,0,g_dim_, g_dim_) = tmp;
    this->G_.block(0,g_dim_,g_dim_,g_dim_) = tmp*dt/2.0;
    this->G_.block(g_dim_,g_dim_,g_dim_,g_dim_)= Eigen::Matrix<double,g_dim_,g_dim_>::Identity()*dt;
    return this->G_;

}

//--------------------------------------------------------------------------------------------------------------------

template <typename tState, typename tTransformation>
void ModelSENPoseTwist<tState,tTransformation>::UpdateState(const Eigen::Matrix<double,2*g_dim_,1>& state_update) {
    this->state_.g_.OPlusEq(state_update.block(0,0,g_dim_,1));
    this->state_.u_.data_ += state_update.block(g_dim_,0,g_dim_,1);
}

} // rransac


#endif // RRANSAC_COMMON_MODELS_SEN_POSE_TWIST_H_