#ifndef RRANSAC_COMMON_MODELS_SEN_POS_VEL_H_
#define RRANSAC_COMMON_MODELS_SEN_POS_VEL_H_
#pragma once


#include "lie_groups/state.h"
#include "lie_groups/utilities.h"
#include "rransac/common/models/model_base.h"
#include "rransac/common/sources/source_SEN_pos_vel.h"
#include "rransac/common/utilities.h"

namespace rransac {

/**
 * \class ModelSENPosVel
 * This model is designed to be used for target's whose configuration manifold is SEN and whose measurement space 
 * is RN. It is assumed that the sources used with this model can only observe the position and translational velocity of the target and that the orientation and
 * angular velocity is not measured. The model is only observable under the assumption that the velocity of the target is oriented with the heading of the target and 
 * that the translation velocity in the body frame is positive. These assumptions enforce that the translational velocities perpindicular to the heading are zero. Because
 * of this, the dimension of the error covariance is less than the dimension of the state. So, in order to properly model the target, we had to make adjustments. 
 */ 
template <typename tState, template <typename > typename tTransformation, template <typename > typename tSource=SourceSENPosVel>
class ModelSENPosVel : public ModelBase<tSource<tState>, tTransformation<tState>,  tState::Group::dim_ + tState::Algebra::dim_ - tState::Algebra::dim_t_vel_ + 1, ModelSENPosVel<tState, tTransformation,tSource>> {

public:

typedef tState State;                                                       /**< The state of the target. @see State. */
typedef typename State::DataType DataType;                                  /**< The scalar object for the data. Ex. float, double, etc. */
typedef tSource<tState> Source;                                             /**< The object type of the source. @see SourceBase. */
typedef tTransformation<tState> Transformation;                             /**< The object type of the measurement and track transformation. */

template <typename tScalar, template<typename> typename tStateTemplate>
using ModelTemplate = ModelSENPosVel<tStateTemplate<tScalar>,tTransformation,tSource>; /**< Used to create a model of the state, source and transformation, but with a different DataType. This is needed to solve the 
                                                                                     nonlinear log maximum likelihood estimation problem by Ceres. */

static constexpr int cov_dim_ = State::Group::dim_ + State::Algebra::dim_ - State::Algebra::dim_t_vel_ + 1;   /**< The dimension of the error covariance. */
static constexpr unsigned int g_dim_ = State::Group::dim_;                                                    /**< The dimension of the pose of the state, i.e. the dimension of the group portion of the state. */                                          
static constexpr unsigned int l_dim_ =  State::Algebra::dim_a_vel_ + 1;                                       /**< The dimension of the angular velocity plus 1. */
typedef Eigen::Matrix<DataType,cov_dim_,cov_dim_> Mat;                                                        /**< The object type of the error covariance, Jacobians, and others. */
typedef Eigen::Matrix<DataType,cov_dim_,1> VecCov;                                                            /**< The object type for the OMinus operation. */


static_assert(std::is_same<typename tSource<tState>::ModelCompatibility, utilities::CompatibleWithModelSENPosVel>::value, "The source is not compatible with the model");
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
void DerivedOPlusEq(const VecCov& state_update);

/**
 * Returns a Random State
 */ 
static State DerivedGetRandomState();

/**
 * Computes the OMinus operation for the state. In other words, it computes the geodesic distance between the two track's state estimate.
 * This OMinus operation can be different than the operation defined by the state. 
 * @param[in] track1 A track of this type.
 * @param[in] track2 A track of this type.
 * @return Returns the geodesic distance between the track's state estimate.
 */ 
static VecCov DerivedOMinus(const ModelSENPosVel& model1, const ModelSENPosVel& model2 ) {

    VecCov tmp;
    Eigen::Matrix<DataType,tState::dim_,1> err = model1.state_.OMinus(model2.state_);
    tmp.block(0,0,tState::Group::dim_+1,1) = err.block(0,0,tState::Group::dim_+1,1);
    tmp.block(tState::Group::dim_+1,0,tState::Group::dim_rot_,1) = err.block(tState::Group::dim_+tState::Group::dim_pos_,0,tState::Group::dim_rot_,1);

    return tmp;
}

};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                            Definitions
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename tState, template <typename > typename tTransformation, template <typename > typename tSource>
typename ModelSENPosVel<tState,tTransformation,tSource>::Mat  ModelSENPosVel<tState,tTransformation,tSource>::DerivedGetLinTransFuncMatState(const State& state, const DataType dt) {  
    Mat F;
    F.block(g_dim_,0,l_dim_,g_dim_).setZero();
    F.block(g_dim_, g_dim_, l_dim_, l_dim_).setIdentity();
    Eigen::Matrix<DataType, g_dim_,g_dim_> tmp = (state.u_*dt).Jr()*dt;  
    F.block(0,0,g_dim_,g_dim_) = typename State::Group(State::Algebra::Exp(state.u_.data_*dt)).Adjoint();
    F.block(0,g_dim_, g_dim_, 1) = tmp.block(0,0,g_dim_,1); // Jacobian w.r.t. rho x
    F.block(0,g_dim_+1, g_dim_, State::Algebra::dim_a_vel_) = tmp.block(0, State::Algebra::dim_t_vel_, g_dim_, State::Algebra::dim_a_vel_); // Jacobian w.r.t. angular velocities
    return F;
}

//--------------------------------------------------------------------------------------------------------------------

template <typename tState, template <typename > typename tTransformation, template <typename > typename tSource>
typename ModelSENPosVel<tState,tTransformation,tSource>::Mat ModelSENPosVel<tState,tTransformation,tSource>::DerivedGetLinTransFuncMatNoise(const State& state, const DataType dt){
    Mat G;
    Eigen::Matrix<DataType, g_dim_,g_dim_> tmp = (state.u_*dt).Jr()*dt; 
    G.block(g_dim_,0,l_dim_,g_dim_).setZero();
    G.block(0,0,g_dim_, g_dim_) = tmp;
    G.block(0,g_dim_, g_dim_, 1) = tmp.block(0,0,g_dim_,1)*dt/2.0;
    G.block(0,g_dim_+1, g_dim_, State::Algebra::dim_a_vel_) = tmp.block(0, State::Algebra::dim_t_vel_, g_dim_, State::Algebra::dim_a_vel_)*dt/2.0; // Jacobian w.r.t. angular velocities

    G.block(g_dim_,g_dim_,l_dim_,l_dim_)= Eigen::Matrix<DataType,l_dim_,l_dim_>::Identity()*dt;
    return G;

}

//--------------------------------------------------------------------------------------------------------------------

template <typename tState, template <typename > typename tTransformation, template <typename > typename tSource>
void ModelSENPosVel<tState,tTransformation,tSource>::DerivedOPlusEq(const VecCov& state_update) {
    
    Eigen::Matrix<DataType,g_dim_,1> twist_update;
    twist_update.setZero();
    twist_update(0,0) = state_update(g_dim_,0);          // get rho_x
    twist_update.block(State::Algebra::dim_t_vel_,0, State::Algebra::dim_a_vel_,1) = state_update.block(g_dim_+1,0,State::Algebra::dim_a_vel_,1);
    this->state_.g_.OPlusEq(state_update.block(0,0,g_dim_,1));
    this->state_.u_.data_ += twist_update;
}

//--------------------------------------------------------------------------------------------------------------------

template <typename tState, template <typename > typename tTransformation, template <typename > typename tSource>
tState ModelSENPosVel<tState,tTransformation,tSource>::DerivedGetRandomState(){
    State state = State::Random();

    state.g_.R_.block(0,0,state.u_.p_.rows(),1) = state.u_.p_.normalized(); 

    if(state.g_.dim_ == 3) {
        state.g_.R_.block(0,1,state.u_.p_.rows(),1) << - state.g_.data_(1,0), state.g_.data_(0,0);
    } else {

        state.g_.R_.block(0,1,state.u_.p_.rows(),1) << -state.g_.data_(1,0), state.g_.data_(0,0), state.g_.data_(2,0);
        state.g_.R_.block(0,2,state.u_.p_.rows(),1) = State::Group::RotAlgebra::Wedge(state.g_.data_.block(0,0,state.u_.p_.rows(),1))*state.g_.data_.block(0,1,state.u_.p_.rows(),1);

    }

    DataType px = state.u_.p_.norm();
    state.u_.p_.setZero();
    state.u_.p_(0,0) = px;
    return  state;
}

} // rransac


#endif // RRANSAC_COMMON_MODELS_SEN_POS_VEL_H_