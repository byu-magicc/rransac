#ifndef RRANSAC_COMMON_MODELS_BASE_H_
#define RRANSAC_COMMON_MODELS_BASE_H_


/**
 * \class ModelBase
 * R-RANSAC is designed to be modular and work with a variety of models. However, this 
 * version of R-RANSAC is extended to work with any Lie group. See the paper ___ for detailed
 * information regarding the system model. 
 * 
 * Note: You can use any model. 
 */ 

class ModelBase: ConsensusSet
{

 

public:
    Eigen::MatrixXd state_;    /** < The estimated state of the phenomenon or target.*/
    Eigen::MatrixXd err_cov_;  /** < The error covariance. */
    double model_likelihood_;  /** < The likelihood that the model represents an actual phenomenon. This value will be between 0 and 1. */
    std::vector<Meas> new_assoc_meas_; /** < Measurements recently associated with the model. These measurements have not been used to update 
                                                     the model. Once an associated measurement has been used to update the model, they are added to the 
                                                     consensus set and removed from this vector.*/
    bool linear_model_;        /** < If the variable is set to true, it is assumed that the system is linear; otherwise, non-linear.*/
    
    
    

     /**
     * Propagates the state estimate forward or backwards in time according to the time interval dt
     * @param[in] dt  The amount of time the state needs to be propagated. The value of dt can be positive or negative. 
     *                a positive value would indicate forward propagation and a negative value would indicate backward propagation.
     * @return Returns the propagated state.
     */ 
    virtual Eigen::MatrixXd PropagateState(const Eigen::MatrixXd& state, const double dt)=0;


    /**
     * Computes the Jacobian of the state transition function with respect to the state evaluated at the current state estimate.
     * @param[in] dt A time interval
     * @return The Jacobian \f$ F_k\f$. 
     */ 
    virtual Eigen::MatrixXd GetLinTransFuncMatState(const double dt)=0;

    /**
     * * Computes the Jacobian of the state transition function with respect to the noise evaluated at the current state estimate.
     * @param[in] dt  A time interval
     * @return Returns the Jacobian \f$ G_k \f$
     */
    virtual Eigen::MatrixXd GetLinTransFuncMatNoise(const double dt)=0;

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
     * Calculates the Jacobian of the observation matrix with respect to the process noise
     * evaluated at the current state estimate.
     * The Jacobian  is dependent on the measurement source.
     * @param 
     * @param source_ID A unique identifier to identify the source. 
     * @return Returns the Jacobian \f$V_k\f$
     */ 
    virtual Eigen::MatrixXd GetLinObsMatProcessNoise(const double dt, const unsigned int source_ID)=0;


    /**
     * Calculates an estimated measurement given a state of the model and the source ID.
     * Since there can be multiple sources, the measurement can be different based on the source. 
     * @param source_ID A unique identifier to identify the source. 
     * @return Returns a measurement \f$ y \f$ based on the state and source.
     */ 
    virtual Eigen::MatrixXd GetEstMeas(const Eigen::MatrixXd& state, const unsigned int source_ID)=0

    /**
     * Using the transformation data provided by the user, this function transforms the state estimate and error covariance
     * from the previous global frame to the current global frame.
     * @param[in] T The transformation data provided by the user. The data should contain the information needed to transform the model.
     * @param[in] dt The time interval between the previous global frame and the current global frame. 
     */ 
    virtual void TransformModel(const Transformation& T, const double dt)=0;

};


#endif // RRANSAC_COMMON_MODEL_BASE_H_