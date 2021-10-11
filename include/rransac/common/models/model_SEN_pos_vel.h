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
template <typename tSourceContainer>
class ModelSENPosVel : public ModelBase<tSourceContainer, tSourceContainer::State::Group::dim_ + tSourceContainer::State::Algebra::dim_ - tSourceContainer::State::Algebra::dim_t_vel_ + 1, ModelSENPosVel> {

public:

    typedef ModelBase<tSourceContainer, tSourceContainer::State::Group::dim_ + tSourceContainer::State::Algebra::dim_ - tSourceContainer::State::Algebra::dim_t_vel_ + 1, ModelSENPosVel> Base;
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


    static constexpr unsigned int g_dim_ = State::Group::dim_;                                                    /**< The dimension of the pose of the state, i.e. the dimension of the group portion of the state. */                                          
    static constexpr unsigned int l_dim_ =  State::Algebra::dim_a_vel_ + 1;                                       /**< The dimension of the angular velocity plus 1. */


    static_assert(std::is_same<typename SourceContainer::ModelCompatibility, utilities::CompatibleWithModelSENPosVel>::value, "The source is not compatible with the model");
    static_assert(lie_groups::utilities::StateIsSEN_seN<typename SourceContainer::State>::value, "The state is not compatible with the model");



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
static State DerivedGetRandomState(const DataType scalar = static_cast<DataType>(1.0));

/**
 * Computes the OMinus operation for the state. In other words, it computes the geodesic distance between the two track's state estimate.
 * This OMinus operation can be different than the operation defined by the state. 
 * @param[in] track1 A track of this type.
 * @param[in] track2 A track of this type.
 * @return Returns the geodesic distance between the track's state estimate.
 */ 
static VecCov DerivedOMinus(const ModelSENPosVel& model1, const ModelSENPosVel& model2 ) {

    VecCov tmp;
    Eigen::Matrix<DataType,State::dim_,1> err = model1.state_.OMinus(model2.state_);
    tmp.block(0,0,g_dim_+1,1) = err.block(0,0,g_dim_+1,1);
    tmp.block(g_dim_+1,0,State::Group::dim_rot_,1) = err.block(g_dim_+State::Group::dim_pos_,0,State::Group::dim_rot_,1);

    return tmp;
}

/**
 * Computes the right Jacobian of the state's Lie group using an element of the 
 * Cartesian algebraic space. 
 * @param u An element of the Cartesian algebraic space.
 */ 
static MatModelCov DerivedJr(const VecCov& u) {
    // typename State::Vec_SC c;
    // c.setZero();
    // c.block(0,0,g_dim_+1,1) = u.block(0,0,g_dim_+1,1);
    // c.block(g_dim_ + State::Group::dim_pos_,0, State::Group::dim_rot_,1) = u.block(g_dim_+1,0, State::Group::dim_rot_,1);

    typename State::Algebra v(u.block(0,0,g_dim_,1));

    MatModelCov jacobian;
    jacobian.setIdentity();
    jacobian.block(0,0,g_dim_,g_dim_) = v.Jr();


    // typename State::Mat_SC jacobian_c = State::Jr(c);

    // jacobian.block(0,0,g_dim_+1,g_dim_+1) = jacobian_c.block(0,0,g_dim_+1,g_dim_+1);
    // jacobian.block(g_dim_+1,g_dim_+1,State::Group::dim_rot_,State::Group::dim_rot_) =


    return jacobian;
}

/**
 * Computes the left Jacobian of the state's Lie group using an element of the 
 * Cartesian algebraic space. 
 * @param u An element of the Cartesian algebraic space.
 */ 
static MatModelCov DerivedJl(const VecCov& u) {
    
    typename State::Algebra v(u.block(0,0,g_dim_,1));
    MatModelCov jacobian;
    jacobian.setIdentity();
    jacobian.block(0,0,g_dim_,g_dim_) = v.Jl();
    return jacobian;

}

/**
 * Computes the inverse of the right Jacobian of the state's Lie group using an element of the 
 * Cartesian algebraic space. 
 * @param u An element of the Cartesian algebraic space.
 */ 
static MatModelCov DerivedJrInv(const VecCov& u) {

    typename State::Algebra v(u.block(0,0,g_dim_,1));
    MatModelCov jacobian;
    jacobian.setIdentity();
    jacobian.block(0,0,g_dim_,g_dim_) = v.JrInv();
    return jacobian;
}

/**
 * Computes the inverse of the left Jacobian of the state's Lie group using an element of the 
 * Cartesian algebraic space. 
 * @param u An element of the Cartesian algebraic space.
 */ 
static MatModelCov DerivedJlInv(const VecCov& u) {
    
    typename State::Algebra v(u.block(0,0,g_dim_,1));
    MatModelCov jacobian;
    jacobian.setIdentity();
    jacobian.block(0,0,g_dim_,g_dim_) = v.JlInv();
    return jacobian;
}

};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                            Definitions
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename _SourceContainer>
typename ModelSENPosVel<_SourceContainer>::MatModelCov  ModelSENPosVel<_SourceContainer>::DerivedGetLinTransFuncMatState(const State& state, const DataType dt) {  
    MatModelCov F;
    F.block(g_dim_,0,l_dim_,g_dim_).setZero();
    F.block(g_dim_, g_dim_, l_dim_, l_dim_).setIdentity();
    Eigen::Matrix<DataType, g_dim_,g_dim_> tmp = (state.u_*dt).Jr()*dt;  
    F.block(0,0,g_dim_,g_dim_) = typename State::Group(State::Algebra::Exp(-state.u_.data_*dt)).Adjoint();
    F.block(0,g_dim_, g_dim_, 1) = tmp.block(0,0,g_dim_,1); // Jacobian w.r.t. rho x
    F.block(0,g_dim_+1, g_dim_, State::Algebra::dim_a_vel_) = tmp.block(0, State::Algebra::dim_t_vel_, g_dim_, State::Algebra::dim_a_vel_); // Jacobian w.r.t. angular velocities
      
    return F;
}

//--------------------------------------------------------------------------------------------------------------------

template <typename _SourceContainer>
typename ModelSENPosVel<_SourceContainer>::MatModelCov ModelSENPosVel<_SourceContainer>::DerivedGetLinTransFuncMatNoise(const State& state, const DataType dt){
    MatModelCov G;
    Eigen::Matrix<DataType, g_dim_,g_dim_> tmp = (state.u_*dt).Jr(); 
    G.block(g_dim_,0,l_dim_,g_dim_).setZero();
    G.block(0,0,g_dim_, g_dim_) = tmp;
    G.block(0,g_dim_, g_dim_, 1) = tmp.block(0,0,g_dim_,1)*0.5;
    G.block(0,g_dim_+1, g_dim_, State::Algebra::dim_a_vel_) = tmp.block(0, State::Algebra::dim_t_vel_, g_dim_, State::Algebra::dim_a_vel_)*0.5; // Jacobian w.r.t. angular velocities
    G.block(g_dim_,g_dim_,l_dim_,l_dim_)= Eigen::Matrix<DataType,l_dim_,l_dim_>::Identity();
    return G;

}


//-------------------------------------------------------------------------------------------------------------------

template <typename _SourceContainer>
void ModelSENPosVel<_SourceContainer>::DerivedOPlus(State& state, const VecCov& state_update) {
    Eigen::Matrix<DataType,g_dim_,1> twist_update;
    twist_update.setZero();
    twist_update(0,0) = state_update(g_dim_,0);          // get rho_x
    twist_update.block(State::Algebra::dim_t_vel_,0, State::Algebra::dim_a_vel_,1) = state_update.block(g_dim_+1,0,State::Algebra::dim_a_vel_,1);
    state.g_.OPlusEq(state_update.block(0,0,g_dim_,1));
    state.u_.data_ += twist_update;
}

//--------------------------------------------------------------------------------------------------------------------

template <typename _SourceContainer>
typename ModelSENPosVel<_SourceContainer>::State ModelSENPosVel<_SourceContainer>::DerivedGetRandomState(const DataType scalar){
    State state = State::Random(scalar);

    // state.g_.R_.block(0,0,state.u_.p_.rows(),1) = state.u_.p_.normalized(); 

    // if(state.g_.dim_ == 3) {
    //     state.g_.R_.block(0,1,state.u_.p_.rows(),1) << - state.g_.data_(1,0), state.g_.data_(0,0);
    // } else {

    //     state.g_.R_.block(0,1,state.u_.p_.rows(),1) << -state.g_.data_(1,0), state.g_.data_(0,0), state.g_.data_(2,0);
    //     state.g_.R_.block(0,2,state.u_.p_.rows(),1) = State::Group::RotAlgebra::Wedge(state.g_.data_.block(0,0,state.u_.p_.rows(),1))*state.g_.data_.block(0,1,state.u_.p_.rows(),1);

    // }
    

    DataType px = state.u_.p_.norm();
    state.u_.p_.setZero();
    state.u_.p_(0,0) = px;
    return  state;
}

} // rransac


#endif // RRANSAC_COMMON_MODELS_SEN_POS_VEL_H_