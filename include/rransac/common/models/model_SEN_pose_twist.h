#ifndef RRANSAC_COMMON_MODELS_SEN_POSE_TWIST_H_
#define RRANSAC_COMMON_MODELS_SEN_POSE_TWIST_H_

#include "common/models/model_base.h"


namespace rransac {

template <typename tState, template <typename > typename tTransformation>
class ModelSENPoseTwist : public ModelBase< SourceSENPoseTwist<tState>, tTransformation<tState>, tState::Group::dim_*2, ModelSENPoseTwist<tState, tTransformation>> {

public:

typedef tState State;
typedef SourceSENPoseTwist<tState> Source;
typedef typename State::DataType DataType;
typedef tTransformation<tState> Transformation;

template <typename tScalar, template<typename> typename tStateTemplate>
using ModelTemplate = ModelSENPoseTwist<tStateTemplate<tScalar>,tTransformation>;

static constexpr unsigned int cov_dim_ = State::Group::dim_*2;
static constexpr unsigned int g_dim_ = State::Group::dim_;
typedef Eigen::Matrix<DataType,2*g_dim_,2*g_dim_> Mat;
typedef Eigen::Matrix<DataType,Eigen::Dynamic,Eigen::Dynamic> MatXd;

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
static State DerivedGetRandomState(){ return State::Random(); }

static Eigen::Matrix<DataType,cov_dim_,1> DerivedOMinus(const ModelSENPoseTwist& model1, const ModelSENPoseTwist& model2 ) {
     return model1.state_.OMinus(model2.state_);
}

};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                            Definitions
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename tState, template <typename > typename tTransformation>
typename ModelSENPoseTwist<tState,tTransformation>::Mat ModelSENPoseTwist<tState,tTransformation>::DerivedGetLinTransFuncMatState(const State& state, const DataType dt) {    
    Mat F;
    F.block(g_dim_,0,g_dim_,g_dim_).setZero();
    F.block(g_dim_,g_dim_,g_dim_,g_dim_).setIdentity();
    F.block(0,0,g_dim_, g_dim_) = typename State::Group(State::Algebra::Exp(state.u_.data_*dt)).Adjoint();
    F.block(0,g_dim_,g_dim_,g_dim_) = (state.u_*dt).Jr()*dt;
    return F;
}

//--------------------------------------------------------------------------------------------------------------------

template <typename tState, template <typename > typename tTransformation>
typename ModelSENPoseTwist<tState,tTransformation>::Mat ModelSENPoseTwist<tState,tTransformation>::DerivedGetLinTransFuncMatNoise(const State& state, const DataType dt){
    Mat G;
    Eigen::Matrix<DataType, g_dim_, g_dim_> tmp = (state.u_*dt).Jr()*dt;
    G.block(g_dim_,0,g_dim_,g_dim_).setZero();
    G.block(0,0,g_dim_, g_dim_) = tmp;
    G.block(0,g_dim_,g_dim_,g_dim_) = tmp*dt/2.0;
    G.block(g_dim_,g_dim_,g_dim_,g_dim_)= Eigen::Matrix<DataType,g_dim_,g_dim_>::Identity()*dt;
    return G;

}

//--------------------------------------------------------------------------------------------------------------------

template <typename tState, template <typename > typename tTransformation>
void ModelSENPoseTwist<tState,tTransformation>::DerivedOPlusEq(const Eigen::Matrix<DataType,2*g_dim_,1>& state_update) {
    this->state_.g_.OPlusEq(state_update.block(0,0,g_dim_,1));
    this->state_.u_.data_ += state_update.block(g_dim_,0,g_dim_,1);
}

} // rransac


#endif // RRANSAC_COMMON_MODELS_SEN_POSE_TWIST_H_