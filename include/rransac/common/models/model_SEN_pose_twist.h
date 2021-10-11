#ifndef RRANSAC_COMMON_MODELS_SEN_POSE_TWIST_H_
#define RRANSAC_COMMON_MODELS_SEN_POSE_TWIST_H_
#pragma once

#include "lie_groups/state.h"
#include "lie_groups/utilities.h"
#include "rransac/common/models/model_base.h"
#include "rransac/common/sources/source_SEN_pose_twist.h"
#include "rransac/common/utilities.h"

namespace rransac {

/**
 * \class ModelSENPoseTwist
 * This model is designed to be used for target's whose configuration manifold is SEN and the measurement space is also SEN. See ModelBase for more detail.
 */ 
template <typename tSourceContainer>
class ModelSENPoseTwist : public ModelBase< tSourceContainer, tSourceContainer::State::Group::dim_*2, ModelSENPoseTwist> {

public:

typedef ModelBase< tSourceContainer, tSourceContainer::State::Group::dim_*2, ModelSENPoseTwist> Base;
    typedef typename Base::State State;                             /**< The state of the target. @see State. */
    typedef typename Base::DataType DataType;                       /**< The scalar object for the data. Ex. float, double, etc. */
    typedef typename Base::SourceContainer SourceContainer;         /**< The object type of the source. @see SourceBase. */
    typedef typename Base::Transformation Transformation;           /**< The object type of the measurement and track transformation. */
    typedef typename Base::TransformDataType TransformDataType;     /**< The data type of the transformation data. */
    typedef typename Base::Measurement Measurement;                 /**< The measurement type. */
    typedef typename Base::MatModelCov MatModelCov;                 /**< The object type of the error covariance, Jacobians, and others. */
    typedef typename Base::VecCov VecCov;                           /**< The object type of state update. */
    typedef typename Base::MatS MatS;                               /**< The object type of the innovation covariance. */
    typedef typename Base::MatH MatH;                               /**< The object type of the Jacobian of the observation function w.r.t. to the state. */
    typedef typename Base::MatV MatV;                               /**< The object type of the Jacobian of the observation function w.r.t. the measurement noise. */
    static constexpr unsigned int cov_dim_ = Base::cov_dim_;        /**< The dimension of the error covariance. */

    static constexpr unsigned int g_dim_ = State::Group::dim_;                      /**< The dimension of the pose of the state, i.e. the dimension of the group portion of the state. */

    static_assert(std::is_same<typename SourceContainer::ModelCompatibility,  utilities::CompatibleWithModelSENPoseTwist>::value, "The source is not compatible with the model");
    static_assert(lie_groups::utilities::StateIsSEN_seN<State>::value, "The state is not compatible with the model");



/**
 * Propagates the state estimate forward or backwards in time according to the time interval dt
 * @param[in] dt  The amount of time the state needs to be propagated. The value of dt can be positive or negative. 
 *                a positive value would indicate forward propagation and a negative value would indicate backward propagation.
 * @return Returns the propagated state.
 */ 
static void DerivedPropagateState(State& state, const DataType dt) {        
    state.g_.OPlusEq(state.u_.data_*dt);
}


/**
 * Computes the Jacobian of the state transition function with respect to the state evaluated at the current state estimate.
 * @param[in] state The state to linearize about
 * @param[in] dt A time interval
 * @return The Jacobian \f$ F_k\f$. 
 */ 
static MatModelCov DerivedGetLinTransFuncMatState(const State& state, const DataType dt);

/**
 * Computes the Jacobian of the state transition function with respect to the noise evaluated at the current state estimate.
 * @param[in] state The state to linearize about
 * @param[in] dt  A time interval
 * @return Returns the Jacobian \f$ G_k \f$
 */
static MatModelCov DerivedGetLinTransFuncMatNoise(const State& state, const DataType dt);

/**
* Update the state of the model using the provided state_update
* @param[in] state_update An element of the lie algebra of the state used to update the state. 
*/
static void DerivedOPlus(State& state, const VecCov& state_update);

/**
 * Returns a Random State
 */ 
static State DerivedGetRandomState(const DataType scalar = static_cast<DataType>(1.0)){ return State::Random(scalar); }

/**
 * Computes the OMinus operation for the state. In other words, it computes the geodesic distance between the two track's state estimate.
 * This OMinus operation can be different than the operation defined by the state. 
 * @param[in] track1 A track of this type.
 * @param[in] track2 A track of this type.
 * @return Returns the geodesic distance between the track's state estimate.
 */ 
static VecCov DerivedOMinus(const ModelSENPoseTwist& model1, const ModelSENPoseTwist& model2 ) {
     return model1.state_.OMinus(model2.state_);
}

/**
 * Computes the right Jacobian of the state's Lie group using an element of the 
 * Cartesian algebraic space. 
 * @param u An element of the Cartesian algebraic space.
 */ 
static MatModelCov DerivedJr(const VecCov& u) {
    return State::Jr(u);
}

/**
 * Computes the left Jacobian of the state's Lie group using an element of the 
 * Cartesian algebraic space. 
 * @param u An element of the Cartesian algebraic space.
 */ 
static MatModelCov DerivedJl(const VecCov& u) {
    return State::Jl(u);
}

/**
 * Computes the inverse of the right Jacobian of the state's Lie group using an element of the 
 * Cartesian algebraic space. 
 * @param u An element of the Cartesian algebraic space.
 */ 
static MatModelCov DerivedJrInv(const VecCov& u) {
    return State::JrInv(u);
}

/**
 * Computes the inverse of the left Jacobian of the state's Lie group using an element of the 
 * Cartesian algebraic space. 
 * @param u An element of the Cartesian algebraic space.
 */ 
static MatModelCov DerivedJlInv(const VecCov& u) {
    return State::JlInv(u);
}

};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                            Definitions
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename _SourceContainer>
typename ModelSENPoseTwist<_SourceContainer>::MatModelCov ModelSENPoseTwist<_SourceContainer>::DerivedGetLinTransFuncMatState(const State& state, const DataType dt) {    
    MatModelCov F;
    F.block(g_dim_,0,g_dim_,g_dim_).setZero();
    F.block(g_dim_,g_dim_,g_dim_,g_dim_).setIdentity();
    F.block(0,0,g_dim_, g_dim_) = typename State::Group(State::Algebra::Exp(-state.u_.data_*dt)).Adjoint();
    F.block(0,g_dim_,g_dim_,g_dim_) = (state.u_*dt).Jr()*dt;
    return F;
}

//--------------------------------------------------------------------------------------------------------------------

template <typename _SourceContainer>
typename ModelSENPoseTwist<_SourceContainer>::MatModelCov ModelSENPoseTwist<_SourceContainer>::DerivedGetLinTransFuncMatNoise(const State& state, const DataType dt){
    MatModelCov G;
    Eigen::Matrix<DataType, g_dim_, g_dim_> tmp = (state.u_*dt).Jr();
    G.block(g_dim_,0,g_dim_,g_dim_).setZero();
    G.block(0,0,g_dim_, g_dim_) = tmp;
    G.block(0,g_dim_,g_dim_,g_dim_) = tmp/2.0;
    G.block(g_dim_,g_dim_,g_dim_,g_dim_)= Eigen::Matrix<DataType,g_dim_,g_dim_>::Identity();
    return G;

}

//--------------------------------------------------------------------------------------------------------------------

template <typename _SourceContainer>
void ModelSENPoseTwist<_SourceContainer>::DerivedOPlus(State& state, const VecCov& state_update) {
    state.g_.OPlusEq(state_update.block(0,0,g_dim_,1));
    state.u_.data_ += state_update.block(g_dim_,0,g_dim_,1);
}

} // rransac


#endif // RRANSAC_COMMON_MODELS_SEN_POSE_TWIST_H_