#ifndef RRANSAC_COMMON_MODELS_BASE_H_
#define RRANSAC_COMMON_MODELS_BASE_H_



#include <Eigen/Dense>
#include <vector>

#include "state.h"
#include "common/measurement/measurement_base.h"
#include "data_structures/consensus_set.h"
#include "common/parameters.h"

namespace rransac {


/**
 * \class ModelBase
 * R-RANSAC is designed to be modular and work with a variety of models. However, this 
 * version of R-RANSAC is extended to work with any Lie group. See the paper ___ for detailed
 * information regarding the system model. 
 * 
 * In order to work with any Lie group, the model base is a template class that requires 
 * The lie group G, lie algebra U and the measurement M
 * 
 * Note: You can use any model. 
 */ 
template <class S, class M> 
class ModelBase
{

typedef Eigen::Matrix<double, S::G::size1_ + S::U::size1_, S.g_.size1_ + S.u_.size1_> Mat;

public:
    lie_groups::State<G,U> state_;    /** < The estimated state of the phenomenon or target.*/
    Mat err_cov_;                     /** < The error covariance. */
    double model_likelihood_;         /** < The likelihood that the model represents an actual phenomenon. This value will be between 0 and 1. */
    std::vector<std::vector<M>> new_assoc_meas_;   /** < Measurements recently associated with the model. These measurements have not been used to update 
                                                     the model. Once an associated measurement has been used to update the model, they are added to the 
                                                     consensus set and removed from this vector. These measurements should have weights assigned to them. 
                                                     Each vector of measurements corresponds to a unique source ID. */
    
    ConsensusSet<M> cs_;              /** < The consensus set. */

    Mat Q_;  /** < Process noise covariance */




     /**
     * Initializes the model
     * @param[in] params  The system parameters specified by the user
     */ 
    virtual void Init(const Parameters& params)=0;
    

     /**
     * Propagates the state estimate forward or backwards in time according to the time interval dt
     * @param[in] dt  The amount of time the state needs to be propagated. The value of dt can be positive or negative. 
     *                a positive value would indicate forward propagation and a negative value would indicate backward propagation.
     * @return Returns the propagated state.
     */ 
    virtual lie_groups::State<G,U> PropagateState(const lie_groups::State<G,U>& state, const double dt)=0;
    
    
    /**
     * Computes the Jacobian of the state transition function with respect to the state evaluated at the current state estimate.
     * @param[in] state The state to linearize about
     * @param[in] dt A time interval
     * @return The Jacobian \f$ F_k\f$. 
     */ 
    virtual Mat GetLinTransFuncMatState(const lie_groups::State<G,U>& state, const double dt)=0;

    /**
     * Computes the Jacobian of the state transition function with respect to the noise evaluated at the current state estimate.
     * @param[in] state The state to linearize about
     * @param[in] dt  A time interval
     * @return Returns the Jacobian \f$ G_k \f$
     */
    virtual Mat GetLinTransFuncMatNoise(const lie_groups::State<G,U>& state, const double dt)=0;

    /**
     * Propagates the state estimate and error covariance to the current time.
     * @param[in] dt  The amount of time the model needs to be propagated.
     */ 
    virtual void PropagateModel(const double dt)=0; 


    /**
     * Uses the newly associated measurements to update the state estimate, error covariance, and consensus set. 
     * @param[in] param Contains all of the user defined parameters including the sources.
     */ 
    virtual void UpdateModel(const Parameters& param)=0;


    /**
     * Calculates the Jacobian of the observation matrix with respect to the state estimate
     * evaluated at the current state estimate.
     * The Jacobian  is dependent on the measurement source.
     * @param source_ID A unique identifier to identify the source. 
     * @return Returns the Jacobian \f$H_k\f$
     */ 
    virtual Eigen::MatrixXd GetLinObsMatState(const unsigned int source_ID)=0;

    /**
     * Calculates the Jacobian of the observation matrix with respect to the measurement noise
     * evaluated at the current state estimate.
     * The Jacobian  is dependent on the measurement source.
     * @param source_ID A unique identifier to identify the source. 
     * @return Returns the Jacobian \f$V_k\f$
     */ 
    virtual Eigen::MatrixXd GetLinObsMatMeasNoise(const unsigned int source_ID)=0;

    /**
     * Calculates the Jacobian of the observation matrix with respect to the sensor noise
     * evaluated at the current state estimate.
     * The Jacobian  is dependent on the measurement source.
     * @param 
     * @param source_ID A unique identifier to identify the source. 
     * @return Returns the Jacobian \f$V_k\f$
     */ 
    virtual Eigen::MatrixXd GetLinObsMatSensorNoise(const double dt, const unsigned int source_ID)=0;


    /**
     * Calculates an estimated measurement given a state of the model and the source ID.
     * Since there can be multiple sources, the measurement can be different based on the source. 
     * @param source_ID A unique identifier to identify the source. 
     * @return Returns a measurement \f$ y \f$ based on the state and source.
     */ 
    virtual Eigen::MatrixXd GetEstMeas(const Eigen::MatrixXd& state, const unsigned int source_ID)=0;

    /**
     * Using the transformation data provided by the user, this function transforms the state estimate and error covariance
     * from the previous global frame to the current global frame.
     * @param[in] T The transformation data provided by the user. The data should contain the information needed to transform the model.
     * @param[in] dt The time interval between the previous global frame and the current global frame. 
     */ 
    virtual void TransformModel(const Transformation& T, const double dt)=0;

};
}

#endif // RRANSAC_COMMON_MODEL_BASE_H_