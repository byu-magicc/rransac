#ifndef RRANSAC_COMMON_MODELS_RN_H_
#define RRANSAC_COMMON_MODELS_RN_H_
#pragma once


#include "lie_groups/state.h"
#include "lie_groups/utilities.h"
#include "rransac/common/models/model_base.h"
#include "rransac/common/sources/source_RN.h"
#include "rransac/common/utilities.h"
#include "rransac/common/models/model_RN_jacobian.h"

namespace rransac {



/**
 * \class ModelRN
 * This model is designed to be used for target's whose configuration manifold is RN and measurement space is RN. See ModelBase for more detail.
 */ 



template <typename _SourceContainer>
class ModelRN : public ModelBase<_SourceContainer, _SourceContainer::State::dim_, ModelRN> {

public:

    typedef ModelBase<_SourceContainer, _SourceContainer::State::dim_, ModelRN> Base;
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


    static constexpr unsigned int num_tangent_spaces_ = State::NumTangentSpaces;    /**< The number of tangent spaces in the group. Ex 1 tangent space is velocity and 2 is velocity and acceleration. */
    static constexpr unsigned int g_dim_ = State::Group::dim_;                      /**< The dimension of the group element of the state. */
    static constexpr unsigned int u_dim_ = g_dim_*num_tangent_spaces_;              /**< The dimension of the total tangent space. */
    typedef ConstructF<DataType,g_dim_,num_tangent_spaces_> FConstructor;           /**< A functor used to construct the Jacobian F. For certain states, this construction is optimized. @see ConstructF. */



    static_assert(std::is_same<typename SourceContainer::ModelCompatibility, typename rransac::utilities::CompatibleWithModelRN>::value, "ModelRN: The source is not compatible with the model");
    static_assert(lie_groups::utilities::StateIsRN_rN<State>::value, "ModelRN: The state is not compatible with the model");



    /**
 * Propagates the state estimate forward or backwards in time according to the time interval dt
 * @param[in] dt  The amount of time the state needs to be propagated. The value of dt can be positive or negative. 
 *                a positive value would indicate forward propagation and a negative value would indicate backward propagation.
 * @return Returns the propagated state.
 */ 
static void DerivedPropagateState(State& state, const DataType dt) { 
    MatModelCov F =  FConstructor::GetF(dt);     
    state.g_.data_ += F.block(0,g_dim_,g_dim_,u_dim_)*state.u_.data_;
    state.u_.data_ = F.block(g_dim_,g_dim_,u_dim_,u_dim_)*state.u_.data_;
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
static State DerivedGetRandomState(const DataType scalar = static_cast<DataType>(1.0)){ return State::Random(scalar);}

/**
 * Computes the OMinus operation for the state. In other words, it computes the geodesic distance between the two track's state estimate.
 * This OMinus operation can be different than the operation defined by the state. 
 * @param[in] track1 A track of this type.
 * @param[in] track2 A track of this type.
 * @return Returns the geodesic distance between the track's state estimate.
 */ 
 static VecCov DerivedOMinus(const ModelRN & model1, const ModelRN& model2 ) {
     return model1.state_.OMinus(model1.state_, model2.state_);
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



 private:

 

};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                            Definitions
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename tSourceContainer>
typename ModelRN<tSourceContainer>::MatModelCov ModelRN<tSourceContainer>::DerivedGetLinTransFuncMatState(const State& state, const DataType dt) {   
    

    return FConstructor::GetF(dt);
}

//--------------------------------------------------------------------------------------------------------------------------

template <typename tSourceContainer>
typename ModelRN<tSourceContainer>::MatModelCov ModelRN<tSourceContainer>::DerivedGetLinTransFuncMatNoise(const State& state, const DataType dt) {
    
    return MatModelCov::Identity();

}

//---------------------------------------------------------------------------------------------------------------------------
template <typename tSourceContainer>
void ModelRN<tSourceContainer>::DerivedOPlus(State& state, const VecCov& state_update) {
    state.g_.data_ += state_update.block(0,0,g_dim_,1);
    state.u_.data_ += state_update.block(g_dim_,0,u_dim_,1);
}












} // rransac


#endif // RRANSAC_COMMON_MODELS_RN_H_