#ifndef RRANSAC_COMMON_MODELS_RN_H_
#define RRANSAC_COMMON_MODELS_RN_H_

#include "common/models/model_base.h"


namespace rransac {

template <typename tState, typename tTransformation>
class ModelRN : public ModelBase<SourceRN<tState>, tTransformation, tState::g_type_::dim_*2, ModelRN<tState, tTransformation>> {

public:

typedef tState State;
typedef tTransformation Transformation;
static constexpr unsigned int cov_dim_ = tState::g_type_::dim_*2;
// static constexpr unsigned int cov_dim_ = tState::g_type_::dim_*2;
static constexpr unsigned int g_dim_ = State::g_type_::dim_;
typedef Eigen::Matrix<double,2*g_dim_,2*g_dim_> Mat;

/**
 * Computes the Jacobian of the state transition function with respect to the state evaluated at the current state estimate.
 * @param[in] state The state to linearize about
 * @param[in] dt A time interval
 * @return The Jacobian \f$ F_k\f$. 
 */ 
Mat DerivedGetLinTransFuncMatState(const State& state, const double dt);

/**
 * Computes the Jacobian of the state transition function with respect to the noise evaluated at the current state estimate.
 * @param[in] state The state to linearize about
 * @param[in] dt  A time interval
 * @return Returns the Jacobian \f$ G_k \f$
 */
Mat DerivedGetLinTransFuncMatNoise(const State& state, const double dt);

/**
* Update the state of the model using the provided state_update
* @param state_update An element of the lie algebra of the state used to update the state. 
*/
void DerivedUpdateState(const Eigen::Matrix<double,2*g_dim_,1>& state_update);

/**
 * Returns a Random State
 */ 
static State DerivedGetRandomState(){ return State::Random();}

/**
 * 
 */
 static Eigen::Matrix<double,cov_dim_,1> DerivedOMinus(const ModelRN& model1, const ModelRN& model2 ) {
     return model1.state_.Ominus(model2);
 }

};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                            Definitions
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename tState, typename tTransformation>
typename ModelRN<tState,tTransformation>::Mat ModelRN<tState,tTransformation>::DerivedGetLinTransFuncMatState(const State& state, const double dt) {    
    this->F_.block(0,g_dim_,g_dim_,g_dim_) = Eigen::Matrix<double,g_dim_,g_dim_>::Identity()*dt;
    return this->F_;
}

//--------------------------------------------------------------------------------------------------------------------------

template <typename State, typename Transformation>
typename ModelRN<State,Transformation>::Mat ModelRN<State,Transformation>::DerivedGetLinTransFuncMatNoise(const State& state, const double dt){
    this->G_.block(0,0,g_dim_, g_dim_) = Eigen::Matrix<double,g_dim_,g_dim_>::Identity()*dt;
    this->G_.block(0,g_dim_,g_dim_,g_dim_) = Eigen::Matrix<double,g_dim_,g_dim_>::Identity()*dt*dt/2;
    this->G_.block(g_dim_,g_dim_,g_dim_,g_dim_)= Eigen::Matrix<double,g_dim_,g_dim_>::Identity()*dt;
    return this->G_;

}

template <typename State, typename Transformation>
void ModelRN<State,Transformation>::DerivedUpdateState(const Eigen::Matrix<double,2*g_dim_,1>& state_update) {
    this->state_.g_.OPlusEq(state_update.block(0,0,g_dim_,1));
    this->state_.u_.data_ += state_update.block(g_dim_,0,g_dim_,1);
}

} // rransac


#endif // RRANSAC_COMMON_MODELS_RN_H_