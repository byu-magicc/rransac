#ifndef RRANSAC_COMMON_MODELS_SEN_POS_VEL_H_
#define RRANSAC_COMMON_MODELS_SEN_POS_VEL_H_

#include "common/models/model_base.h"


namespace rransac {

template <typename tState, template <typename > typename tTransformation>
class ModelSENPosVel : public ModelBase<SourceSENPosVel<tState>, tTransformation<tState>,  tState::Group::dim_ + tState::Algebra::dim_ - tState::Algebra::dim_t_vel_ + 1, ModelSENPosVel<tState, tTransformation>> {

public:

typedef tState State;
typedef typename State::DataType DataType;
typedef tTransformation<tState> Transformation;

template <typename tScalar, template<typename> typename tStateTemplate>
using ModelTemplate = ModelSENPosVel<tStateTemplate<tScalar>,tTransformation>;

static constexpr unsigned int g_dim_ = State::Group::dim_;
static constexpr unsigned int cov_dim_ = State::Group::dim_ + State::Algebra::dim_ - State::Algebra::dim_t_vel_ + 1;
static constexpr unsigned int l_dim_ =  State::Algebra::dim_a_vel_ + 1;
typedef Eigen::Matrix<DataType,cov_dim_,cov_dim_> Mat;

/**
 * Computes the Jacobian of the state transition function with respect to the state evaluated at the current state estimate.
 * @param[in] state The state to linearize about
 * @param[in] dt A time interval
 * @return The Jacobian \f$ F_k\f$. 
 */ 
static Mat DerivedGetLinTransFuncMatState(const State& state, const double dt);

/**
 * Computes the Jacobian of the state transition function with respect to the noise evaluated at the current state estimate.
 * @param[in] state The state to linearize about
 * @param[in] dt  A time interval
 * @return Returns the Jacobian \f$ G_k \f$
 */
static Mat DerivedGetLinTransFuncMatNoise(const State& state, const double dt);

/**
* Update the state of the model using the provided state_update. The state_update provided is augmented to account
* for the additional translationals velocities (which are zero) before being added to state.
* @param state_update An element of the lie algebra of the state used to update the state. 
*/
void DerivedOPlusEq(const Eigen::Matrix<DataType,cov_dim_,1>& state_update);

/**
 * Returns a Random State
 */ 
static State DerivedGetRandomState();

static Eigen::Matrix<DataType,cov_dim_,1> DerivedOMinus(const ModelSENPosVel& model1, const ModelSENPosVel& model2 ) {

    Eigen::Matrix<DataType,cov_dim_,1> tmp;
    Eigen::Matrix<DataType,tState::dim_,1> err = model1.state_.OMinus(model2.state_);
    tmp.block(0,0,tState::Group::dim_+1,1) = err.block(0,0,tState::Group::dim_+1,1);
    tmp.block(tState::Group::dim_+1,0,tState::Group::dim_rot_,1) = err.block(tState::Group::dim_+tState::Group::dim_pos_,0,tState::Group::dim_rot_,1);

    return tmp;
}

};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                            Definitions
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename tState, template <typename > typename tTransformation>
typename ModelSENPosVel<tState, tTransformation>::Mat  ModelSENPosVel<tState, tTransformation>::DerivedGetLinTransFuncMatState(const State& state, const double dt) {  
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

template <typename tState, template <typename > typename tTransformation>
typename ModelSENPosVel<tState, tTransformation>::Mat ModelSENPosVel<tState, tTransformation>::DerivedGetLinTransFuncMatNoise(const State& state, const double dt){
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

template <typename tState, template <typename > typename tTransformation>
void ModelSENPosVel<tState, tTransformation>::DerivedOPlusEq(const Eigen::Matrix<DataType,cov_dim_,1>& state_update) {
    Eigen::Matrix<DataType,g_dim_,1> twist_update;
    twist_update.setZero();
    twist_update(0,0) = state_update(g_dim_,0);          // get rho_x
    twist_update.block(State::Algebra::dim_t_vel_,0, State::Algebra::dim_a_vel_,1) = state_update.block(g_dim_+1,0,State::Algebra::dim_a_vel_,1);
    this->state_.g_.OPlusEq(state_update.block(0,0,g_dim_,1));
    this->state_.u_.data_ += twist_update;
}

//--------------------------------------------------------------------------------------------------------------------

template <typename tState, template <typename > typename tTransformation>
tState ModelSENPosVel<tState, tTransformation>::DerivedGetRandomState(){
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