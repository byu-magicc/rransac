#ifndef RRANSAC_COMMON_MODELS_RN_H_
#define RRANSAC_COMMON_MODELS_RN_H_
#pragma once



#include "common/models/model_base.h"
#include "common/utilities.h"
#include "state.h"
#include "utilities.h"


namespace rransac {

template <typename tState, template <typename > typename tSource, template <typename > typename tTransformation>
class ModelRN : public ModelBase<tSource<tState>, tTransformation<tState>, tState::Group::dim_*2, ModelRN<tState, tSource, tTransformation>> {

public:

typedef tState State;
typedef tSource<tState> Source;
typedef typename State::DataType DataType;
typedef tTransformation<tState> Transformation;

template <typename tScalar, template<typename> typename tStateTemplate>
using ModelTemplate = ModelRN<tStateTemplate<tScalar>,tSource,tTransformation>;

static constexpr unsigned int cov_dim_ = tState::Group::dim_*2;
static constexpr unsigned int g_dim_ = State::Group::dim_;
typedef Eigen::Matrix<DataType,2*g_dim_,2*g_dim_> Mat;

typedef utilities::CompatibleWithModelRN ModelCompatibility;
static_assert(std::is_same<typename tSource<tState>::ModelCompatibility, ModelCompatibility>::value, "The source is not compatible with the model");
static_assert(lie_groups::utilities::StateIsRN_rN<tState>::value, "The state is not compatible with the model");

/**
 * Computes the Jacobian of the state transition function with respect to the state evaluated at the current state estimate.
 * @param[in] state The state to linearize about
 * @param[in] dt A time interval
 * @return The Jacobian \f$ F_k\f$. 
 */ 
static Mat DerivedGetLinTransFuncMatState(const State& state, const DataType dt);

/**
 * Computes the Jacobian of the state transition function with respect to the noise evaluated at the current state estimate.
 * @param[in] state The state to linearize about
 * @param[in] dt  A time interval
 * @return Returns the Jacobian \f$ G_k \f$
 */
static Mat DerivedGetLinTransFuncMatNoise(const State& state, const DataType dt);

/**
* Update the state of the model using the provided state_update
* @param state_update An element of the lie algebra of the state used to update the state. 
*/
void DerivedOPlusEq(const Eigen::Matrix<DataType,2*g_dim_,1>& state_update);

/**
 * Returns a Random State
 */ 
static State DerivedGetRandomState(){ return State::Random();}

/**
 * 
 */
 static Eigen::Matrix<DataType,cov_dim_,1> DerivedOMinus(const ModelRN & model1, const ModelRN& model2 ) {
     return model1.state_.OMinus(model1.state_, model2.state_);
 }

};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                            Definitions
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename tState, template <typename > typename tSource, template <typename > typename tTransformation>
typename ModelRN<tState,tSource,tTransformation>::Mat ModelRN<tState,tSource,tTransformation>::DerivedGetLinTransFuncMatState(const State& state, const DataType dt) {   
    
    // static constexpr unsigned int g_dim_ = tState::Group::dim_;
    // typedef Eigen::Matrix<double,2*g_dim_,2*g_dim_> Mat;
    
    Mat F = Mat::Identity();
    F.block(0,g_dim_,g_dim_,g_dim_) = Eigen::Matrix<DataType,g_dim_,g_dim_>::Identity()*dt;
    return F;
}

//--------------------------------------------------------------------------------------------------------------------------

template <typename tState, template <typename > typename tSource, template <typename > typename tTransformation>
typename ModelRN<tState,tSource,tTransformation>::Mat ModelRN<tState,tSource,tTransformation>::DerivedGetLinTransFuncMatNoise(const State& state, const DataType dt) {
    
    Mat G;
    G.block(g_dim_,0,g_dim_,g_dim_).setZero();
    G.block(0,0,g_dim_, g_dim_) = Eigen::Matrix<DataType,g_dim_,g_dim_>::Identity()*dt;
    G.block(0,g_dim_,g_dim_,g_dim_) = Eigen::Matrix<DataType,g_dim_,g_dim_>::Identity()*dt*dt/2;
    G.block(g_dim_,g_dim_,g_dim_,g_dim_)= Eigen::Matrix<DataType,g_dim_,g_dim_>::Identity()*dt;
    return G;

}

//--------------------------------------------------------------------------------------------------------------------------

template <typename tState, template <typename > typename tSource, template <typename > typename tTransformation>
void ModelRN<tState,tSource,tTransformation>::DerivedOPlusEq(const Eigen::Matrix<DataType,2*g_dim_,1>& state_update) {
    this->state_.g_.OPlusEq(state_update.block(0,0,g_dim_,1));
    this->state_.u_.data_ += state_update.block(g_dim_,0,g_dim_,1);
}

} // rransac


#endif // RRANSAC_COMMON_MODELS_RN_H_