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
typedef Eigen::Matrix<double,2.0*S::G::dim_,2.0*S::G::dim_> Mat;


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
        F.block(0,0,State::G::dim_, State::G::dim_) = State::G::Adjoint(State::U::Exp(state.u_.data_*dt));
        F.block(0,State::G::dim_,State::G::dim_,State::G::dim_) = State::U::Jr(state.u_.data_*dt)*dt;
        F.block(State::G::dim_,0,State::G::dim_,State::G::dim_).setZero();
        F.block(State::G::dim_,State::G::dim_,State::G::dim_,State::G::dim_).setIdentity(); 
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
        G.block(0,0,State::G::dim_, State::G::dim_) = State::U::Jr(state.u_.data_*dt)*dt;
        G.block(0,State::G::dim_,State::G::dim_,State::G::dim_) = State::U::Jr(state.u_.data_*dt)*dt*dt/2;
        G.block(State::G::dim_,0,State::G::dim_,State::G::dim_).setZero(); 
        G.block(State::G::dim_,State::G::dim_,State::G::dim_,State::G::dim_)= Eigen::Matrix<double,State::G::dim_,State::G::dim_>::Identity()*dt;
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

void GetInnovationAndCovarianceFixedMeasurementCov(const std::vector<Meas>& meas, Eigen::Matrix<double,State::g_type_::dim_*2,1>& innovation, Eigen::Matrix<double,State::g_type_::dim_*2,State::g_type_::dim_*2>& cov);

void GetInnovationAndCovarianceNonFixedMeasurementCov(const std::vector<Meas>& meas, Eigen::Matrix<double,State::g_type_::dim_*2,1>& innovation, Eigen::Matrix<double,State::g_type_::dim_*2,State::g_type_::dim_*2>& cov);


};


// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// //                                               Abelian Models
// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// template <typename State, typename Source, typename Transformation> 
// void ModelBase<State, lie_groups::Abelian, Source, Transformation>::Init(std::vector<Source>& sources, const Parameters& params) {
//     sources_ = sources;
//     F_.setIdentity();
//     G_.setIdentity(); 
// }

// //---------------------------------------------------------------------------------------------------------



template <typename State, typename SourceType, typename Source, typename Transformation> 
void ModelBase<State, SourceType, Source, Transformation>::UpdateModel(const Parameters& params) {

Eigen::Matrix<double,State::g_type_::dim_*2,1> innovation_sum;
Eigen::Matrix<double,State::g_type_::dim_*2,1> innovation;
Eigen::Matrix<double,State::g_type_::dim_*2,1> update;
Eigen::Matrix<double,State::g_type_::dim_*2,State::g_type_::dim_*2> cov;
Eigen::Matrix<double,State::g_type_::dim_*2,State::g_type_::dim_*2> cov_sum;
Eigen::Matrix<double,State::g_type_::dim_*2,State::g_type_::dim_*2> error_cov_inverse = err_cov_.inverse();
innovation_sum.setZero();
cov_sum.setZero();

// loop through the measurements per source
for (std::vector<Meas> meas : new_assoc_meas_) {

    if(sources_[meas.front()].params_.meas_cov_fixed_) {
        GetInnovationAndCovarianceFixedMeasurementCov(meas, innovation, cov);
    } else {
        GetInnovationAndCovarianceNonFixedMeasurementCov(meas, innovation, cov);
    }

    innovation_sum+= innovation;
    cov_sum += (cov.inverse() - error_cov_inverse);
}

error_cov_inverse = error_cov_inverse + cov_sum;
err_cov_ = error_cov_inverse.inverse();

update = err_cov_ * innovation_sum;
state_.g_.OPlusEq(update.block(0,0,State::g_type_::dim_,1));
state_.u_.data += update.block(State::g_type_::dim_,0,State::g_type_::dim_,1)

}

//---------------------------------------------------------------------------------------------------------

template <typename State, typename SourceType, typename Source, typename Transformation> 
void ModelBase<State, SourceType, Source, Transformation>::GetInnovationAndCovarianceFixedMeasurementCov(const std::vector<Meas>& meas, Eigen::Matrix<double,State::g_type_::dim_*2,1>& innovation, Eigen::Matrix<double,State::g_type_::dim_*2,State::g_type_::dim_*2>& cov) {

innovation.setZero();
covariance.setZero();

// std::vector<Eigen::Matrix<double,State::g_type_::dim_*2,1>> innovations(meas.size());
Eigen::Matrix<double,State::g_type_::dim_*2,1> error
Eigen::MatrixXd H = GetLinObsMatState(state_);
Eigen::MatrixXd V = GetLinObsMatMeasNoise(state_);
Eigen::MatrixXd K;
Eigen::MatrixXd tmp = V*sources_[meas.front().source_index].params_.meas_cov_*V.transpose();
Eigen::Matrix<double,State::g_type_::dim_*2,State::g_type_::dim_*2> cov_tilde=Eigen::Matrix<double,State::g_type_::dim_*2,State::g_type_::dim_*2>::Zero();
Eigen::Matrix<double,State::g_type_::dim_*2,State::g_type_::dim_*2> innovation_cov;
Meas estimated_meas = sources_[meas.front().source_index].GetEstMeas(state_);

innovation_cov = (H*err_cov_*H.transpose() + tmp);
K = err_cov_*H.transpose()*innovation_cov.inverse();

double B0 = 1;

// Get total innovation
for (Meas m : meas) {
    error = sources_[m.source_index].OMinus(m, estimated_meas);
    innovation += m.weight*error;
    cov_tilde+= m.weight*error*error.transpose();
    B0 -= m.weight;
}

// Construct cov tilde
cov_tilde -= innovation*innovation.transpose()*meas.size();
cov_tilde = K*cov_tilde*K.transpose();

// construct covariance
cov = err_cov_ - (1-B0)*K*innovation_cov*K.transpose();


}


} // namespace rransac

#endif // RRANSAC_COMMON_MODEL_BASE_H_