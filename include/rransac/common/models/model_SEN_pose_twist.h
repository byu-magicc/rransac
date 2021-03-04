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
 * \class ModelRN
 * This model is designed to be used for target's whose configuration manifold is SEN and the measurement space is also SEN. See ModelBase for more detail.
 */ 
template <typename tState, template <typename > typename tTransformation, template <typename > typename tSource = SourceSENPoseTwist>
class ModelSENPoseTwist : public ModelBase< SourceSENPoseTwist<tState>, tTransformation<tState>, tState::Group::dim_*2, ModelSENPoseTwist<tState, tTransformation,tSource>> {

public:

typedef tState State;                                                       /**< The state of the target. @see State. */
typedef typename State::DataType DataType;                                  /**< The scalar object for the data. Ex. float, double, etc. */
typedef tSource<tState> Source;                                             /**< The object type of the source. @see SourceBase. */
typedef tTransformation<tState> Transformation;                             /**< The object type of the measurement and track transformation. */

template <typename tScalar, template<typename> typename tStateTemplate>
using ModelTemplate = ModelSENPoseTwist<tStateTemplate<tScalar>,tTransformation,tSource>; /**< Used to create a model of the state, source and transformation, but with a different DataType. This is needed to solve the 
                                                                                     nonlinear log maximum likelihood estimation problem by Ceres. */

static constexpr unsigned int cov_dim_ = State::Group::dim_*2;                 /**< The dimension of the error covariance. */
static constexpr unsigned int g_dim_ = State::Group::dim_;                      /**< The dimension of the pose of the state, i.e. the dimension of the group portion of the state. */
typedef Eigen::Matrix<DataType,2*g_dim_,2*g_dim_> Mat;                          /**< The object type of the error covariance, Jacobians, and others. */

static_assert(std::is_same<typename tSource<tState>::ModelCompatibility,  utilities::CompatibleWithModelSENPoseTwist>::value, "The source is not compatible with the model");
static_assert(lie_groups::utilities::StateIsSEN_seN<tState>::value, "The state is not compatible with the model");


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
static State DerivedGetRandomState(){ return State::Random(); }

/**
 * Computes the OMinus operation for the state. In other words, it computes the geodesic distance between the two track's state estimate.
 * This OMinus operation can be different than the operation defined by the state. 
 * @param[in] track1 A track of this type.
 * @param[in] track2 A track of this type.
 * @return Returns the geodesic distance between the track's state estimate.
 */ 
static Eigen::Matrix<DataType,cov_dim_,1> DerivedOMinus(const ModelSENPoseTwist& model1, const ModelSENPoseTwist& model2 ) {
     return model1.state_.OMinus(model2.state_);
}

};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                            Definitions
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename tState, template <typename > typename tTransformation, template <typename > typename tSource>
typename ModelSENPoseTwist<tState,tTransformation,tSource>::Mat ModelSENPoseTwist<tState,tTransformation,tSource>::DerivedGetLinTransFuncMatState(const State& state, const DataType dt) {    
    Mat F;
    F.block(g_dim_,0,g_dim_,g_dim_).setZero();
    F.block(g_dim_,g_dim_,g_dim_,g_dim_).setIdentity();
    F.block(0,0,g_dim_, g_dim_) = typename State::Group(State::Algebra::Exp(-state.u_.data_*dt)).Adjoint();
    F.block(0,g_dim_,g_dim_,g_dim_) = (state.u_*dt).Jr()*dt;
    return F;
}

//--------------------------------------------------------------------------------------------------------------------

template <typename tState, template <typename > typename tTransformation, template <typename > typename tSource>
typename ModelSENPoseTwist<tState,tTransformation,tSource>::Mat ModelSENPoseTwist<tState,tTransformation,tSource>::DerivedGetLinTransFuncMatNoise(const State& state, const DataType dt){
    Mat G;
    Eigen::Matrix<DataType, g_dim_, g_dim_> tmp = (state.u_*dt).Jr();
    G.block(g_dim_,0,g_dim_,g_dim_).setZero();
    G.block(0,0,g_dim_, g_dim_) = tmp;
    G.block(0,g_dim_,g_dim_,g_dim_) = tmp/2.0;
    G.block(g_dim_,g_dim_,g_dim_,g_dim_)= Eigen::Matrix<DataType,g_dim_,g_dim_>::Identity();
    return G;

}

//--------------------------------------------------------------------------------------------------------------------

template <typename tState, template <typename > typename tTransformation, template <typename > typename tSource>
void ModelSENPoseTwist<tState,tTransformation,tSource>::DerivedOPlusEq(const Eigen::Matrix<DataType,2*g_dim_,1>& state_update) {
    this->state_.g_.OPlusEq(state_update.block(0,0,g_dim_,1));
    this->state_.u_.data_ += state_update.block(g_dim_,0,g_dim_,1);
}

} // rransac


#endif // RRANSAC_COMMON_MODELS_SEN_POSE_TWIST_H_