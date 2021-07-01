#ifndef RRANSAC_COMMON_MODELS_RN_H_
#define RRANSAC_COMMON_MODELS_RN_H_
#pragma once


#include "lie_groups/state.h"
#include "lie_groups/utilities.h"
#include "rransac/common/models/model_base.h"
#include "rransac/common/sources/source_RN.h"
#include "rransac/common/utilities.h"






/**
 * \class ModelRN
 * This model is designed to be used for target's whose configuration manifold is RN and measurement space is RN. See ModelBase for more detail.
 */ 

namespace rransac {

template <typename tSourceContainer>
class ModelRN : public ModelBase<tSourceContainer, tSourceContainer::State::Group::dim_*2, ModelRN> {

public:

typedef tSourceContainer SourceContainer;                                   /**< The Source container. */
typedef typename tSourceContainer::State State;                             /**< The state of the target. @see State. */
typedef typename State::DataType DataType;                                  /**< The scalar object for the data. Ex. float, double, etc. */
typedef typename SourceContainer::Transformation Transformation;            /**< The transformation data type. */

static constexpr unsigned int cov_dim_ = State::Group::dim_*2;                  /**< The dimension of the error covariance. */
static constexpr unsigned int g_dim_ = State::Group::dim_;                      /**< The dimension of the pose of the state, i.e. the dimension of the group portion of the state. */
typedef Eigen::Matrix<DataType,2*g_dim_,2*g_dim_> Mat;                          /**< The object type of the error covariance, Jacobians, and others. */




static_assert(std::is_same<typename SourceContainer::ModelCompatibility, typename rransac::utilities::CompatibleWithModelRN>::value, "ModelRN: The source is not compatible with the model");
static_assert(lie_groups::utilities::StateIsRN_rN<State>::value, "ModelRN: The state is not compatible with the model");

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
* @param[in] state_update An element of the lie algebra of the state used to update the state. 
*/
void DerivedOPlusEq(const Eigen::Matrix<DataType,2*g_dim_,1>& state_update);

/**
 * Returns a Random State
 */ 
static State DerivedGetRandomState(){ return State::Random();}

/**
 * Computes the OMinus operation for the state. In other words, it computes the geodesic distance between the two track's state estimate.
 * This OMinus operation can be different than the operation defined by the state. 
 * @param[in] track1 A track of this type.
 * @param[in] track2 A track of this type.
 * @return Returns the geodesic distance between the track's state estimate.
 */ 
 static Eigen::Matrix<DataType,cov_dim_,1> DerivedOMinus(const ModelRN & model1, const ModelRN& model2 ) {
     return model1.state_.OMinus(model1.state_, model2.state_);
 }

};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                            Definitions
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename tSourceContainer>
typename ModelRN<tSourceContainer>::Mat ModelRN<tSourceContainer>::DerivedGetLinTransFuncMatState(const State& state, const DataType dt) {   
    
    Mat F = Mat::Identity();
    F.block(0,g_dim_,g_dim_,g_dim_) = Eigen::Matrix<DataType,g_dim_,g_dim_>::Identity()*dt;
    return F;
}

//--------------------------------------------------------------------------------------------------------------------------

template <typename tSourceContainer>
typename ModelRN<tSourceContainer>::Mat ModelRN<tSourceContainer>::DerivedGetLinTransFuncMatNoise(const State& state, const DataType dt) {
    
    Mat G;
    G.block(g_dim_,0,g_dim_,g_dim_).setZero();
    G.block(0,0,g_dim_, g_dim_) = Eigen::Matrix<DataType,g_dim_,g_dim_>::Identity();
    G.block(0,g_dim_,g_dim_,g_dim_) = Eigen::Matrix<DataType,g_dim_,g_dim_>::Identity()/2.0;
    G.block(g_dim_,g_dim_,g_dim_,g_dim_)= Eigen::Matrix<DataType,g_dim_,g_dim_>::Identity();
    return G;

}

//--------------------------------------------------------------------------------------------------------------------------

template <typename tSourceContainer>
void ModelRN<tSourceContainer>::DerivedOPlusEq(const Eigen::Matrix<DataType,2*g_dim_,1>& state_update) {
    this->state_.g_.OPlusEq(state_update.block(0,0,g_dim_,1));
    this->state_.u_.data_ += state_update.block(g_dim_,0,g_dim_,1);
}

} // rransac


#endif // RRANSAC_COMMON_MODELS_RN_H_