#ifndef RRANSAC_COMMON_MODELS_BASE_H_
#define RRANSAC_COMMON_MODELS_BASE_H_



#include <Eigen/Dense>
#include <vector>

#include "state.h"
#include "common/measurement/measurement_base.h"
#include "data_structures/consensus_set.h"
#include "parameters.h"
#include "common/sources/source_base.h"
#include "common/sources/source_RN.h"
#include "common/sources/source_SEN_pos_vel.h"
#include "common/sources/source_SEN_pose_twist.h"

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
template <typename State, typename Source, typename Transformation, int Cov_DIM, typename Derived> 
class ModelBase
{
public:
typedef State MB_State;
// typedef StateType MB_StateType;
typedef Source MB_Source;
typedef Transformation MB_Transformation;
static constexpr unsigned int cov_dim_ = Cov_DIM;
typedef Eigen::Matrix<double,cov_dim_,cov_dim_> Mat;



public:
    State state_;                     /** < The estimated state of the phenomenon or target.*/
    Mat err_cov_;                     /** < The error covariance. */
    double model_likelihood_;         /** < The likelihood that the model represents an actual phenomenon. This value will be between 0 and 1. */
    std::vector<std::vector<Meas>> new_assoc_meas_;   /** < Measurements recently associated with the model. These measurements have not been used to update 
                                                     the model. Once an associated measurement has been used to update the model, they are added to the 
                                                     consensus set and removed from this vector. These measurements should have weights assigned to them. 
                                                     Each vector of measurements corresponds to a unique source ID. */
    
    ConsensusSet<Meas> cs_;           /** < The consensus set. */

    Mat Q_;                           /** < Process noise covariance */

    std::vector<Source>* sources_;    /** < Reference to the sources contained in system */

    // Useful matrices for propagating
    Mat F_;                           /** < The Jacobian of the state transition function w.r.t. the states */
    Mat G_;                           /** < The Jacobian of the state transition function w.r.t. the noise  */


    ModelBase()=default;
    ~ModelBase()=default;

     /**
     * Initializes the model.
     * @param[in] sources A reference to the sources contained in system
     * @param[in] params  The system parameters specified by the user
     */ 
    void Init(std::vector<Source>& sources, const Parameters& params);
    
    /**
     * Sets the user defined parameters
     * @param[in] params  The system parameters specified by the user
     */ 
    void SetParameters(const Parameters& params) { Q_ = params.process_noise_covariance_; }

     /**
     * Propagates the state estimate forward or backwards in time according to the time interval dt
     * @param[in] dt  The amount of time the state needs to be propagated. The value of dt can be positive or negative. 
     *                a positive value would indicate forward propagation and a negative value would indicate backward propagation.
     * @return Returns the propagated state.
     */ 
    static State PropagateState(const State& state, const double dt) {
        State tmp = state;
        tmp.g_.OPlusEq(tmp.u_.data_*dt);
        return tmp;
    }
    
    
    /**
     * Computes the Jacobian of the state transition function with respect to the state evaluated at the current state estimate.
     * @param[in] state The state to linearize about
     * @param[in] dt A time interval
     * @return The Jacobian \f$ F_k\f$. 
     */ 
    Mat GetLinTransFuncMatState(const State& state, const double dt) {
        return static_cast<Derived*>(this)->GetLinTransFuncMatState(state, dt);
    }

    /**
     * Computes the Jacobian of the state transition function with respect to the noise evaluated at the current state estimate.
     * @param[in] state The state to linearize about
     * @param[in] dt  A time interval
     * @return Returns the Jacobian \f$ G_k \f$
     */
    Mat GetLinTransFuncMatNoise(const State& state, const double dt){
        return static_cast<Derived*>(this)->GetLinTransFuncMatNoise(state, dt);

    }

    /**
     * Propagates the state estimate and error covariance to the current time.
     * @param[in] dt  The amount of time the model needs to be propagated.
     */ 
    void PropagateModel(const double dt);

    /**
     * Uses the newly associated measurements to update the statMeas_dime estimate, error covariance, and consensus set using a 
     * centralized measurement fusion.
     * @param[in] param Contains all of the user defined parameters.
     */ 
    void UpdateModel(const Parameters& params) {
        static_cast<Derived*>(this)->UpdateState(GetStateUpdate( params));
    }


    /**
     * Calculates the Jacobian of the observation matrix with respect to the state estimate
     * evaluated at the current state estimate.
     * The Jacobian  is dependent on the measurement source.
     * @param source_ID A unique identifier to identify the source. 
     * @return Returns the Jacobian \f$H_k\f$
     */ 
    Eigen::MatrixXd GetLinObsMatState(const State& state, const unsigned int source_ID){
        return (*sources_)[source_ID].GetLinObsMatState(state);
    }

    /**
     * Calculates the Jacobian of the observation matrix with respect to the measurement noise
     * evaluated at the current state estimate.
     * The Jacobian  is dependent on the measurement source.
     * @param source_ID A unique identifier to identify the source. 
     * @return Returns the Jacobian \f$V_k\f$
     */ 
    Eigen::MatrixXd GetLinObsMatSensorNoise(const State& state, const unsigned int source_ID){
        return (*sources_)[source_ID].GetLinObsMatSensorNoise(state);
    }

    /**
     * Calculates an estimated measurement given a state of the model and the source ID.
     * Since there can be multiple sources, the measurement can be different based on the source. 
     * @param source_ID A unique identifier to identify the source. 
     * @return Returns a measurement \f$ y \f$ based on the state and source.
     */ 
    Eigen::MatrixXd GetEstMeas(const State& state, const unsigned int source_ID){
        return (*sources_)[source_ID].GetEstMeas(state);
    }

    /**
     * Using the transformation data provided by the user, this function transforms the state estimate and error covariance
     * from the previous global frame to the current global frame.
     * @param[in] T The transformation object provided by the user. The object should already have the data it needs to transform the model.
     * @param[in] dt The time interval between the previous global frame and the current global frame. 
     */ 
    void TransformModel(Transformation& T){
        T.TransformTrack(state_,this->err_cov_);
    }

    /**
     * Returns a Random State
     */ 
    static State GetRandomState(){
        return Derived::GetRandomState();
    }

private:

void GetInnovationAndCovarianceFixedMeasurementCov(const std::vector<Meas>& meas, Eigen::Matrix<double,cov_dim_,1>& state_update, Mat& cov_update);

void GetInnovationAndCovarianceNonFixedMeasurementCov(const std::vector<Meas>& meas, Eigen::Matrix<double,cov_dim_,1>& state_update, Mat& cov_update);

Eigen::Matrix<double,Cov_DIM,1> GetStateUpdate(const Parameters& params);


};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                            Definitions
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename State, typename Source, typename Transformation, int Cov_DIM, typename Derived>
void ModelBase<State, Source, Transformation,Cov_DIM, Derived>::Init(std::vector<Source>& sources, const Parameters& params) {
    sources_ = &sources;
    err_cov_.setIdentity();
    F_.setIdentity();
    G_.setIdentity();
    SetParameters(params);
}

//-------------------------------------------------------------------------------------------------------------------

template <typename State, typename Source, typename Transformation, int Cov_DIM, typename Derived> 
void ModelBase<State, Source, Transformation,Cov_DIM, Derived>::PropagateModel(const double dt) {

    // Construct matrices to transform covariance.
    F_ = GetLinTransFuncMatState(state_,dt);
    G_ = GetLinTransFuncMatNoise(state_,dt);


    // Transform covariance
    err_cov_ = F_*err_cov_*F_.transpose() + G_*Q_*G_.transpose();

    // Propagate state
    state_.g_.OPlusEq(state_.u_.data_*dt);

}

//-------------------------------------------------------------------------------------------------------------------

template <typename State, typename Source, typename Transformation, int Cov_DIM, typename Derived> 
Eigen::Matrix<double,Cov_DIM,1> ModelBase<State, Source, Transformation,Cov_DIM, Derived>::GetStateUpdate(const Parameters& params) {

Eigen::Matrix<double,cov_dim_,1> state_update_sum;
Eigen::Matrix<double,cov_dim_,1> state_update;
Eigen::Matrix<double,cov_dim_,1> update;
Mat cov;
Mat cov_sum;
Mat error_cov_inverse = err_cov_.inverse();
state_update_sum.setZero();
cov_sum.setZero();

// std::cout << "state: " << std::endl << state_.g_.data_ << std::endl;

// loop through the measurements per source
for (std::vector<Meas> meas : new_assoc_meas_) {

    if((*sources_)[meas.front().source_index].params_.meas_cov_fixed_) {
        GetInnovationAndCovarianceFixedMeasurementCov(meas, state_update, cov);
    } else {
        GetInnovationAndCovarianceNonFixedMeasurementCov(meas, state_update, cov);
    }

    // std::cout << "state_update: " << std::endl << state_update << std::endl;
    // std::cout << "cov: " << std::endl << cov << std::endl;


    state_update_sum+= state_update;
    cov_sum += (cov.inverse() - error_cov_inverse);
    // cov_sum += cov;
}

// std::cout << "cov_sum: " << std::endl << cov_sum << std::endl;



// Update the state
update = err_cov_ * state_update_sum;
// std::cout << "state_update_sum: " << std::endl << state_update_sum << std::endl;
// std::cout << "update: " << std::endl << update << std::endl;


// Update the error covariance
error_cov_inverse += cov_sum;
err_cov_ = error_cov_inverse.inverse();

return update;

}

//---------------------------------------------------------------------------------------------------------

template <typename State, typename Source, typename Transformation, int Cov_DIM, typename Derived> 
void ModelBase<State, Source, Transformation,Cov_DIM, Derived>::GetInnovationAndCovarianceFixedMeasurementCov(const std::vector<Meas>& meas, Eigen::Matrix<double,cov_dim_,1>& state_update, Mat& cov_update) {

state_update.setZero();
cov_update = err_cov_;


Eigen::MatrixXd H = GetLinObsMatState(state_,meas.front().source_index);              // Jacobian of observation function w.r.t. state
Eigen::MatrixXd V = GetLinObsMatSensorNoise(state_,meas.front().source_index);        // Jacobian of observation function w.r.t. noise
Eigen::MatrixXd K;                                                                    // Kalman Gain
Eigen::MatrixXd S_inverse;                                                            // Innovation covariance inverse
Eigen::MatrixXd nu_i;                                                                 
Eigen::MatrixXd nu(V.rows(),1);                                                       // Total innovation term
nu.setZero();
Eigen::MatrixXd covSum(V.rows(),V.rows());
covSum.setZero();
// std::cerr << "H: " << std::endl << H << std::endl << std::endl;
// std::cerr << "V: " << std::endl << V << std::endl << std::endl;
// std::cerr << "R: " << std::endl << (*sources_)[meas.front().source_index].params_.meas_cov_ << std::endl << std::endl;
Meas estimated_meas = (*sources_)[meas.front().source_index].GetEstMeas(state_);
// Eigen::MatrixXd tmp = (V*(*sources_)[meas.front().source_index].params_.meas_cov_*V.transpose()).inverse();


// std::cerr << "tmp: " << std::endl << tmp << std::endl << std::endl;
// std::cerr << "H*err_cov_*H.transpose() : " << std::endl << H*err_cov_*H.transpose()  << std::endl << std::endl;


S_inverse = (H*err_cov_*H.transpose() + V*(*sources_)[meas.front().source_index].params_.meas_cov_*V.transpose()).inverse();
K = err_cov_*H.transpose()*S_inverse;

// std::cerr << "K: " << std::endl << K << std::endl << std::endl;


double B0 = 1;

// Get total weighted innovation and part of the cov_tilde
for (Meas m : meas) {
    nu_i = (*sources_)[m.source_index].OMinus(m, estimated_meas);
    nu += m.weight*nu_i;
    covSum+= m.weight*nu_i*nu_i.transpose();
    B0 -= m.weight;

// std::cerr << "meas: " << std::endl << m.pose << std::endl << std::endl;
// std::cerr << "nu_i: " << std::endl << nu_i << std::endl << std::endl;

}

// Finish constructing cov_sum
covSum -= nu*nu.transpose();

// std::cerr << "covS: " << std::endl << covSum << std::endl << std::endl;
// std::cerr << "tmp: " << std::endl << tmp << std::endl << std::endl;
// std::cerr << "nu: " << std::endl << nu << std::endl << std::endl;
// std::cerr << "Q: " << std::endl << Q_ << std::endl << std::endl;

// construct covariance
cov_update +=  K*(covSum*K.transpose() -(1-B0)*H*err_cov_); 
// cov_update = H.transpose()*tmp*H;
state_update = H.transpose()*S_inverse*nu;


}


//---------------------------------------------------------------------------------------------------------
template <typename State, typename Source, typename Transformation, int Cov_DIM, typename Derived> 
void ModelBase<State, Source, Transformation,Cov_DIM, Derived>::GetInnovationAndCovarianceNonFixedMeasurementCov(const std::vector<Meas>& meas, Eigen::Matrix<double,cov_dim_,1>& state_update, Mat& cov_update) {

state_update.setZero();
cov_update = err_cov_;


std::vector<Eigen::MatrixXd> nu_i(meas.size());                                  // Innovation term for each measurement
Eigen::MatrixXd H = GetLinObsMatState(state_,meas.front().source_index);         // Jacobian of observation function w.r.t. state
Eigen::MatrixXd V = GetLinObsMatSensorNoise(state_,meas.front().source_index);   // Jacobian of observation function w.r.t. noise
Eigen::MatrixXd K;                                                               // Kalman Gain
Eigen::MatrixXd S_inverse;                                                       // Innovation covariance inverse
Eigen::MatrixXd nu(err_cov_.rows(),1);                                         // Total innovation term
nu.setZero();
Meas estimated_meas = (*sources_)[meas.front().source_index].GetEstMeas(state_);

// std::cerr << "err_cov_: " << std::endl << err_cov_ << std::endl << std::endl;
// std::cerr << "H: " << std::endl << H << std::endl << std::endl;


// Helper calculations
Eigen::MatrixXd H_cov_HT = H*err_cov_*H.transpose();
Eigen::MatrixXd cov_HT = err_cov_*H.transpose();

// Get individual innovations and weighted total innovations
for (unsigned long int ii=0; ii< meas.size(); ++ii) {
    nu_i[ii] = (*sources_)[meas[ii].source_index].OMinus(meas[ii], estimated_meas);
    

}

// Construct covariance update and state update
for (unsigned long int ii=0; ii< meas.size(); ++ii) {
    
    S_inverse = (H_cov_HT + V*meas[ii].meas_cov*V.transpose()).inverse();
    K = cov_HT*S_inverse;
    nu += meas[ii].weight*K*nu_i[ii];
    cov_update += meas[ii].weight*(K*nu_i[ii]*nu_i[ii].transpose() - cov_HT)*K.transpose();

    state_update += meas[ii].weight*H.transpose()*S_inverse*nu_i[ii];
}

cov_update -= nu*nu.transpose();



}


} // namespace rransac

#endif // RRANSAC_COMMON_MODEL_BASE_H_