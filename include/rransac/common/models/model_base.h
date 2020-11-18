#ifndef RRANSAC_COMMON_MODELS_BASE_H_
#define RRANSAC_COMMON_MODELS_BASE_H_



#include <Eigen/Dense>
#include <vector>

#include "state.h"
#include "common/measurement/measurement_base.h"
#include "data_structures/consensus_set.h"
#include "parameters.h"

namespace rransac {


/**
 * \class ModelBase
 * R-RANSAC is designed to be modular and work with a variety of models. However, this 
 * version of R-RANSAC is extended to work with any Lie group. See the paper ___ for detailed
 * information regarding the system model. 
 * 
 * In order to work with any Lie group, the model base is a template class that requires 
 * The state S
 * 
 * Note: You can use any model. 
 */ 
template <typename State, typename StateType, typename Source, typename Transformation> 
class ModelBase
{

typedef State State;
typedef StateType StateType;
typedef Source Source;
typedef Transformation Transformation;
static constexpr unsigned int g_dim_ = State::g_type_::dim_;
typedef Eigen::Matrix<double,2.0*g_dim_,2.0*g_dim_> Mat;



public:
    State state_;                     /** < The estimated state of the phenomenon or target.*/
    Mat err_cov_;                 /** < The error covariance. */
    double model_likelihood_;         /** < The likelihood that the model represents an actual phenomenon. This value will be between 0 and 1. */
    std::vector<std::vector<Meas>> new_assoc_meas_;   /** < Measurements recently associated with the model. These measurements have not been used to update 
                                                     the model. Once an associated measurement has been used to update the model, they are added to the 
                                                     consensus set and removed from this vector. These measurements should have weights assigned to them. 
                                                     Each vector of measurements corresponds to a unique source ID. */
    
    ConsensusSet<Meas> cs_;           /** < The consensus set. */

    Mat Q_;                           /** < Process noise covariance */

    std::vector<Source>& sources_;    /** < Reference to the sources contained in system */

    // Useful matrices for propagating
    Mat F_;                           /** < The Jacobian of the state transition function w.r.t. the states */
    Mat G_;                           /** < The Jacobian of the state transition function w.r.t. the noise  */


     /**
     * Initializes the model.
     * @param[in] sources A reference to the sources contained in system
     * @param[in] params  The system parameters specified by the user
     */ 
    void Init(std::vector<Source>& sources, const Parameters& params) {
        sources_ = sources;
        // static_cast<Derived*>(this)->Init(params);
    }
    
    /**
     * Sets the user defined parameters
     * @param[in] params  The system parameters specified by the user
     */ 
    void SetParameters(const Parameters& params) {
        Q_ = params.process_noise_covariance_;
    }

     /**
     * Propagates the state estimate forward or backwards in time according to the time interval dt
     * @param[in] dt  The amount of time the state needs to be propagated. The value of dt can be positive or negative. 
     *                a positive value would indicate forward propagation and a negative value would indicate backward propagation.
     * @return Returns the propagated state.
     */ 
    static State PropagateState(const State& state, const double dt) {
        State tmp = state;
        // lie_groups::SE2_se2 tmp;
        tmp.g_.OPlusEq(tmp.u_.data_*dt);
        return tmp;
    }
    
    
    /**
     * Computes the Jacobian of the state transition function with respect to the state evaluated at the current state estimate.
     * @param[in] state The state to linearize about
     * @param[in] dt A time interval
     * @return The Jacobian \f$ F_k\f$. 
     */ 
    static Mat GetLinTransFuncMatState(const State& state, const double dt) {
        Mat F;
        F.block(0,0,g_dim_, g_dim_) = State::g_type_::Adjoint(State::u_type_::Exp(state.u_.data_*dt));
        F.block(0,g_dim_,g_dim_,g_dim_) = State::u_type_::Jr(state.u_.data_*dt)*dt;
        F.block(g_dim_,0,g_dim_,g_dim_).setZero();
        F.block(g_dim_,g_dim_,g_dim_,g_dim_).setIdentity(); 
        return F;
    }

    /**
     * Computes the Jacobian of the state transition function with respect to the noise evaluated at the current state estimate.
     * @param[in] state The state to linearize about
     * @param[in] dt  A time interval
     * @return Returns the Jacobian \f$ G_k \f$
     */
    static Mat GetLinTransFuncMatNoise(const State& state, const double dt){
        Mat G;
        G.block(0,0,g_dim_, g_dim_) = State::u_type_::Jr(state.u_.data_*dt)*dt;
        G.block(0,g_dim_,g_dim_,g_dim_) = State::u_type_::Jr(state.u_.data_*dt)*dt*dt/2;
        G.block(g_dim_,0,g_dim_,g_dim_).setZero(); 
        G.block(g_dim_,g_dim_,g_dim_,g_dim_)= Eigen::Matrix<double,g_dim_,g_dim_>::Identity()*dt;
        return G;

    }

    /**
     * Propagates the state estimate and error covariance to the current time.
     * @param[in] dt  The amount of time the model needs to be propagated.
     */ 
    void PropagateModel(const double dt) {

        // Construct matrices to transform covariance.
        F_ = GetLinTransFuncMatState(state_,dt);
        G_ = GetLinTransFuncMatNoise(state_,dt);

        // Transform covariance
        err_cov_ = F_*err_cov_*F_.transpose() + G_*err_cov_*G_.transpose();

        // Propagate state
        state_.g_.OPlusEq(state_.u_.data_*dt);

    }

    /**
     * Uses the newly associated measurements to update the state estimate, error covariance, and consensus set using a 
     * centralized measurement fusion.
     * @param[in] param Contains all of the user defined parameters.
     */ 
    void UpdateModel(const Parameters& params);


    /**
     * Calculates the Jacobian of the observation matrix with respect to the state estimate
     * evaluated at the current state estimate.
     * The Jacobian  is dependent on the measurement source.
     * @param source_ID A unique identifier to identify the source. 
     * @return Returns the Jacobian \f$H_k\f$
     */ 
    static Eigen::MatrixXd GetLinObsMatState(const State& state, const unsigned int source_ID){
        return sources_[source_ID].GetLinObsMatState(state);
    }

    /**
     * Calculates the Jacobian of the observation matrix with respect to the measurement noise
     * evaluated at the current state estimate.
     * The Jacobian  is dependent on the measurement source.
     * @param source_ID A unique identifier to identify the source. 
     * @return Returns the Jacobian \f$V_k\f$
     */ 
    static Eigen::MatrixXd GetLinObsMatMeasNoise(const State& state, const unsigned int source_ID){
        return sources_[source_ID].GetLinObsMatMeasNoise(state);
    }

    /**
     * Calculates an estimated measurement given a state of the model and the source ID.
     * Since there can be multiple sources, the measurement can be different based on the source. 
     * @param source_ID A unique identifier to identify the source. 
     * @return Returns a measurement \f$ y \f$ based on the state and source.
     */ 
    Eigen::MatrixXd GetEstMeas(const S& state, const unsigned int source_ID){
        return sources_[source_ID].GetEstMeas(state);
    }

    /**
     * Using the transformation data provided by the user, this function transforms the state estimate and error covariance
     * from the previous global frame to the current global frame.
     * @param[in] T The transformation object provided by the user. The object should already have the data it needs to transform the model.
     * @param[in] dt The time interval between the previous global frame and the current global frame. 
     */ 
    virtual void TransformModel(const Transformation& T){
        T.TransformTrack(state_,err_cov_);
    }

private:

void GetInnovationAndCovarianceFixedMeasurementCov(const std::vector<Meas>& meas, Eigen::Matrix<double,g_dim_*2,1>& state_update, Mat& cov_update);

void GetInnovationAndCovarianceNonFixedMeasurementCov(const std::vector<Meas>& meas, Eigen::Matrix<double,g_dim_*2,1>& state_update, Mat& cov_update);


};

 //---------------------------------------------------------------------------------------------------------



template <typename State, typename SourceType, typename Source, typename Transformation> 
void ModelBase<State, SourceType, Source, Transformation>::UpdateModel(const Parameters& params) {

Eigen::Matrix<double,g_dim_*2,1> state_update_sum;
Eigen::Matrix<double,g_dim_*2,1> state_update;
Eigen::Matrix<double,g_dim_*2,1> update;
Mat cov;
Mat cov_sum;
Mat error_cov_inverse = err_cov_.inverse();
state_update_sum.setZero();
cov_sum.setZero();

// loop through the measurements per source
for (std::vector<Meas> meas : new_assoc_meas_) {

    if(sources_[meas.front()].params_.meas_cov_fixed_) {
        GetInnovationAndCovarianceFixedMeasurementCov(meas, state_update, cov);
    } else {
        GetInnovationAndCovarianceNonFixedMeasurementCov(meas, state_update, cov);
    }

    state_update_sum+= state_update;
    cov_sum += (cov.inverse() - error_cov_inverse);
}

// Update the error covariance
error_cov_inverse += cov_sum;
err_cov_ = error_cov_inverse.inverse();

// Update the state
update = err_cov_ * state_update_sum;
state_.g_.OPlusEq(update.block(0,0,g_dim_,1));
state_.u_.data += update.block(g_dim_,0,g_dim_,1)

}

//---------------------------------------------------------------------------------------------------------

template <typename State, typename SourceType, typename Source, typename Transformation> 
void ModelBase<State, SourceType, Source, Transformation>::GetInnovationAndCovarianceFixedMeasurementCov(const std::vector<Meas>& meas, Eigen::Matrix<double,g_dim_*2,1>& state_update, Mat& cov_update) {

state_update.setZero();
cov_update = err_cov_;

Eigen::Matrix<double,g_dim_*2,1> nu = Eigen::Matrix<double,g_dim_*2,1>::Zero(); // Total innovation term
Eigen::Matrix<double,g_dim_*2,1> nu_i;
Eigen::MatrixXd H = GetLinObsMatState(state_);                                         // Jacobian of observation function w.r.t. state
Eigen::MatrixXd V = GetLinObsMatMeasNoise(state_);                                     // Jacobian of observation function w.r.t. noise
Eigen::MatrixXd K;                                                                     // Kalman Gain
Mat S;                                                                                 // Innovation covariance
Mat cov_tilde=Mat::Zero();
Meas estimated_meas = sources_[meas.front().source_index].GetEstMeas(state_);
Eigen::MatrixXd tmp = V*sources_[meas.front().source_index].params_.meas_cov_*V.transpose();


S = (H*err_cov_*H.transpose() + tmp);
K = err_cov_*H.transpose()*S.inverse();

double B0 = 1;

// Get total weighted innovation and part of the cov_tilde
for (Meas m : meas) {
    nu_i = sources_[m.source_index].OMinus(m, estimated_meas);
    nu += m.weight*nu_i;
    cov_tilde+= m.weight*nu_i*nu_i.transpose();
    B0 -= m.weight;
}

// Finish constructing cov_tilde
cov_tilde -= nu*nu.transpose()*meas.size();
cov_tilde = K*cov_tilde*K.transpose();

// construct covariance
cov_update +=  cov_tilde - (1-B0)*K*S*K.transpose();
state_update = H.transpose()*tmp.inverse()*nu;


}


//---------------------------------------------------------------------------------------------------------
template <typename State, typename SourceType, typename Source, typename Transformation> 
void ModelBase<State, SourceType, Source, Transformation>::GetInnovationAndCovarianceNonFixedMeasurementCov(const std::vector<Meas>& meas, Eigen::Matrix<double,g_dim_*2,1>& state_update, Mat& cov_update) {

state_update.setZero();
cov_update = err_cov_;

Eigen::Matrix<double,g_dim_*2,1> nu = Eigen::Matrix<double,g_dim_*2,1>::Zero();        // Total innovation term
std::vector<Eigen::Matrix<double,g_dim_*2,1>> nu_i(meas.size());                       // Innovation term for each measurement
Eigen::MatrixXd H = GetLinObsMatState(state_);                                         // Jacobian of observation function w.r.t. state
Eigen::MatrixXd V = GetLinObsMatMeasNoise(state_);                                     // Jacobian of observation function w.r.t. noise
Eigen::MatrixXd K;                                                                     // Kalman Gain
Mat S;                                                                                 // Innovation covariance
Meas estimated_meas = sources_[meas.front().source_index].GetEstMeas(state_);


// Get individual innovations and weighted total innovations
for (unsigned long int ii=0; ii< meas.size(); ++ii) {
    nu_i[ii] = sources_[meas[ii].source_index].OMinus(meas[ii], estimated_meas);
    nu += meas[ii].weight*nu_i[ii];
}

// Construct covariance update and state update
for (unsigned long int ii=0; ii< meas.size(); ++ii) {

    S = (H*err_cov_*H.transpose() + V*meas[ii].meas_cov*V.transpose());
    K = err_cov_*H.transpose()*S.inverse();
    cov_update += meas[ii].weight*K*nu_i[ii]*nu_i[ii].transpose()-nu*nu.transpose()*K.transpose() - meas[ii].weight*K*S*K.transpose();

    state_update += meas[ii].weight*H.transpose()*(V*meas[ii].meas_cov*V.transpose()).inverse()*nu_i[ii];
}



}


} // namespace rransac

#endif // RRANSAC_COMMON_MODEL_BASE_H_