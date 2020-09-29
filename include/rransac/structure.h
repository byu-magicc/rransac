

//////////////////////////////////////////////////
//          transformation.h
//////////////////////////////////////////////////

#include <Eigen/Core>

/** \struct Transformation
 * When the global frame changes, the measurements and models need to be transformed into the current global frame.
 * This struct provides the necessary data in order to transform the measurements and models. It is provided by the
 * user and used by R-RANSAC.
*/


struct Transformation
{
    Eigen::MatrixXd T;
};



//////////////////////////////////////////////////
//          measurement_base.h
//////////////////////////////////////////////////
#include <Eigen/Core>

#include "rransac/common/transformation.h"
#include "rransac/common/parameters.h"


/** \struct Meas
 * This struct contains information regarding a measurement. The user uses this object to supply R-RANSAC with measurements. 
 * The user is responsible to provide Meas::data, Meas::time_stamp and Meas::source_id. If the measurement covariance is not provided 
 * to the Source object, the user must also provide Meas::meas_cov. The other member variables are set by R-RANSAC.
*/

struct Meas
{
    Eigen::MatrixXd data;       /**< The data that represents the measurement. */
    double time_stamp;          /**< The time the measurement was taken. */
    unsigned int source_id;     /**< A unique identifier that indicates which source the measurement came from. */
    Eigen::MatrixXf meas_cov;   /**< The measurement covariance. Only used if the measurement covariance changes with different measurements; otherwise, the 
                                     measurement covariance given to the Source class is used for every measurement. */ 


    double likelihood;          /**< The likelihood that the measurement came from the phenomenon it was associated with. This value is set during the data 
                                      association process.*/
    double weight;              /**< The weight of the measurement when updating the model is was associated with. This value is set during the data association 
                                     process. */
    bool associated=false;      /**< Indicates if the measurement has been associated with a model. It should only be modified by the Model Initializer.;

};


/** \enum Meas
 * Used to indicate what type of distance between measurements to calculated.
*/
enum DistanceType
{
    kSpatial=0,   
    kTemporal=1,
    kTotal=2
};

/** \enum DerivedMeasurement

 * Used to indicate what type of derived measurement class to create. 
*/
enum DerivedMeasurement
{
    kMeasurement=0
}


/** \class MeasurementBase
 * The class is responsible for performing any calculation that involves only measurements. 
 * All of the member functions in this class should be static. 
*/

class MeasurementBase
{
public:
    /**
     * Transforms a measurement from the previous global frame to the current global frame.
     * @param[in/out] meas The measurement that will be transformed.
     * @param[in] Transformation The transformation object that contains the data necessary to transform the measurement.
     * @param[in] dt The time elapsed since the last transformation was received.
     */
    virtual void TransformMeasurement(Meas& meas,const Transformation& T,const double dt)=0;


    /**
     * Calculates the distance between two measurements depending on the type of distance to calculate. 
     * @param[in] meas1 A measurement.
     * @param[in] meas2 A different measurement.
     * @param[in] type The type of distance to be calculated: spatial, temporal, etc. 
     * @param[in] params Contains all of the user defined parameters. A user can define a weight when calculating the distances. 
     */
    virtual float GetDistance(const Meas& meas1, const Meas& meas2, const DistanceType& type, const Parameters& params)=0;

};



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//          measurement.h           ///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "rransac/common/measurement/measurement_base.h"

/** \class Measurement
 * The default derived measurement class.
*/

class Measurement: public MeasurementBase
{
public:
    /**
     * Transforms a measurement from the previous global frame to the current global frame.
     * @param[in/out] meas The measurement that will be transformed.
     * @param[in] Transformation The transformation object that contains the data necessary to transform the measurement.
     * @param[in] dt The time elapsed since the last transformation was received.
     */
    static void TransformMeasurement(Meas& meas,const Transformation& T,const double dt);


    /**
     * Calculates the distance between two measurements depending on the type of distance to calculate.
     * @param[in] meas1 A measurement.
     * @param[in] meas2 A different measurement.
     * @param[in] type The type of distance to be calculated: spatial, temporal, or total 
     * @param[in] params Contains all of the user defined parameters. A user can define a weight when calculating the distances. 
     */
    static float GetDistance(const Meas& meas1, const Meas& meas2, const DistanceType& type, const Parameters& params);

private:

    /**
     * Calculates the spatial distance between two measurements using the 2-norm.
     * @param[in] meas1 A measurement.
     * @param[in] meas2 A different measurement.
     * @param[in] params Contains all of the user defined parameters. A user can define a weight when calculating the distances.
     */
    static float GetSpatialDistance(Meas& meas1, Meas& meas2, const Parameters& params);

    /**
     * Calculates the temporal distance between two measurements using the 2-norm.
     * @param[in] meas1 A measurement.
     * @param[in] meas2 A different measurement.
     * @param[in] params Contains all of the user defined parameters. A user can define a weight when calculating the distances.
     */
    static float GetTemporalDistance(Meas& meas1, Meas& meas2, const Parameters& params);

    /**
     * Calculates the temporal and spatial distance between two measurements using the 2-norm.
     * @param[in] meas1 A measurement.
     * @param[in] meas2 A different measurement.
     * @param[in] params Contains all of the user defined parameters. A user can define a weight when calculating the distances.
     */
    static float GetTotalDistance(Meas& meas1, Meas& meas2, const Parameters& params);
 

};


//////////////////////////////////////////////////
//          measurement_factory.h
//////////////////////////////////////////////////

#include "rransac/common/measurement/measurement.h

/** \class MeasurementFactory
 * Creates a new instant of a derived measurement class.
*/

class MeasurementFactory 
{
public:

    /**
     * Returns a unique pointer to a derived measurement class.
     * @param[in] type The type of derived measurement class to use. 
     * @see DerivedMeasurement
     */
    std::unique_ptr<MeasurementBase> CreateMeasurementClass(const DerivedMeasurement type);
};


//////////////////////////////////////////////////
//          source.h
//////////////////////////////////////////////////

#include <Eigen/Core>



/** \class Source 
 * A source is an algorithm that takes data from a sensor and extracts measurements from it. Due to the nature of
 * the source, there is measurement noise and false measurements. The measurement noise is assumed to be sampled from
 * a white-noise, zero-mean, Gaussian distribution. If the measurement covariance is fixed, then it should be added to 
 * the source; otherwise, it is supplied with each measurement. We assume that the false measurements are uniformly
 * distributed in the local surveillance region, and that the number
 * of expected measurements per sensor scan from a source can be modeled using a Poisson distribution. 
 * 
 * Since there can be many different sources, we need to be able to distinguish the different sources using a unique ID.
 * When a measurement is received, the measurement indicates which source it came from using the source ID tag @see Meas.
 * We also need to know the dimension of the measurement.
 * 
 */ 

class Source
{

public:

    bool meas_cov_fixed_;       /** < Flag used to indicate if the measurement covariance is the same for every measurement. 
                                      If it is, then the measurement covariance needs to be given to the source. If it isn't, 
                                      then the measurement covariance needs to be given with the measurement. @see Meas */
    Eigen::MatrixXd meas_cov_;  /** < The fixed measurement covariance.*/
    float expected_num_false_meas_;  /** < The expected number of false measurements. We assume that the expected number of false 
                                        measurements from a source per sensor scan can be modeled using a Poisson distribution. */
    unsigned int source_id_;    /** < An identifier that is unique for every source. When a measurement is given, the measurement
                                      indicates which source it came from using the source ID. @see Meas */

    unsigned int meas_dim_; /** < The dimension of the measurement. You can also consider this as the number of generalized coordinates measured.
                                  For example, if the measurement containted the position of a 2D object in the
                                  xy-plane, then the dimension of the measurement would be two. This is used by RANSAC when creating model hypotheses. */

    
};





//////////////////////////////////////////////////
//          consensus_set.h
//////////////////////////////////////////////////

#include "rransac/common/measurement/measurement.h

/**
 * \class ConsensusSet
 * The consensus set of a model is the set of measurements associated with the model.
 * In order to keep the consensus set becoming arbitrarily large. Expired measurements are 
 * pruned from the consensus set. The data member consensus_set is a list of vectors that contain
 * the measurements. Each vector of measurements contains measurements with the same time stamp. 
 */ 

class ConsensusSet
{
public:

/** < 
 * Add a measurement to the consensus set. 
 * @param[in] meas The measurement to be added.
*/
void AddMeasToConsensusSet(Meas meas);

/** <
 * Removes all of the measurements from the consensus_set with a time stamp that occurred before expiration_time. 
 * The value of expiration_time is the current time in seconds minus the time window. I.e. old measurements are removed.
 * @param[in] expiration_time The current time in seconds minus the time window.
 */ 

void PruneConsensusSet(double expiration_time);


std::list<std::vector<Meas>> consensus_set_; /** < Contains the measurements associated with the model that have not expired. Each vector of measurements 
                                                contains measurements with the same time stamp. This allows us to efficiently remove expired measurements
                                                by simply removing an entire vector. */
};




//////////////////////////////////////////////////
//          model_base.h
//////////////////////////////////////////////////

#include <Eigen/Core>

#include "rransac/common/transformation.h"
#include "rransac/data_structures/consensus_set.h"
#include "rransac/common/measurement/measurement_base.h"


/**
 * \struct SimilarSources
 * A phenomenon or target has a number of degrees of freedom (or generalized coordinates). For example, an airplane has 6 generalized coordinates; 3
 * for position and 3 for attitude. The state of the model consists of the generalized coordinates \f$ q_i\f$ and their temporal derivatives. A sensor
 * measures some or all of the generalized coordinates. 
 */ 
struct SimilarSources
{
    std::vector<unsigned int> id;
    unsigned int num_samples;
};



/**
 * \enum DerivedModels
 * Lists the different types of models available.
 */ 
enum DerivedModels
{
    kModelLCV=0, /**< Indicates that the linear constant velocity model should be used.*/
    kModelLCA=1, /**< Indicates that the linear constant acceleration model should be used. */
    kModelSE2=2, /**< Indicates that the SE2 model should be used. Note that this is a constant velocity model.*/
    kModelSE3=3  /**< Indicates that the SE3 model should be used. Note that this is a constant velocity model. */
};


/**
 * \class ModelBase
 * R-RANSAC is designed to be modular and work with a variety of models. However,
 * we assume that the dynamic model is a discrete, time-invariant system of the form 
 * \f[
 * x_k=f\left(x_{k^-},w_{k},\delta_k\right)
 * \f]
 * \f[
 * y_k=h\left(x_k\right)+v_k
 * \f]
 * with \f$x_k\f$ denoting the state of the system at time \f$k\f$, \f$x_{k^-}\f$ the previous state of the system, 
 * \f$y_k\f$ the system output at time \f$k\f$, \f$f\f$ the system model, \f$h\f$ the observation model, \f$w\f$ zero-mean 
 * Gaussian process noise with covariance \f$Q\f$, \f$v\f$ zero-mean Gaussian measurement noise with covariance 
 * \f$R\f$ and \f$\delta_k\f$ denoting the time between the previous state and the current state.
 * 
 * The class is designed to facilitate propagating and updating the model using a method similar to the Kalman filter update.
 * Methods similar to the Kalman filter include the EKF, UKF, IF, EIF, PDA, JPDA, etc. Reglardless of the method used, the
 * system model \f$f\f$ needs to be linearized with respect to the state and the noise. Also, the observation model \f$h\f$
 * needs to be linearized with respect to the state. 
 * 
 * In most methods, the following matrices will be helpful. 
 * Let \f$F\f$ denote the linearized system model with respect to the state calculated as 
 * \f[
 *     F=\frac{\partial f}{\partial x}\bigg |_{x_k,\delta_t},
 * \f]
 * \f$G\f$ denote the linearized system model with respect to the noise calculated as 
 * \f[
 *     G=\frac{\partial f}{\partial w}\bigg |_{x_k,\delta_t},
 * \f] and
 * \f$H\f$ denote the linearized observation model with respect to the noise calculated as 
 * \f[
 *     H=\frac{\partial h}{\partial x}\bigg |_{x_k,\delta_t}.
 * \f]
 * 
 *  
 * 
 * This is a fairly generic model that should work for most cases. It is the responsibility of the user
 * to create a derived model class that works for their specific application. 
 */ 

class ModelBase: ConsensusSet
{

//TODO:: Allow the user to add constraints for the states. 

public:
    Eigen::MatrixXd state_;    /** < The estimated state of the phenomenon or target.*/
    Eigen::MatrixXd err_cov_;  /** < The error covariance. */
    double model_likelihood_;  /** < The likelihood that the model represents an actual phenomenon. This value will be between 0 and 1. */
    std::vector<Meas> new_assoc_meas_; /** < Measurements recently associated with the model. These measurements have not been used to update 
                                                     the model. Once an associated measurement has been used to update the model, they are added to the 
                                                     consensus set and removed from this vector.*/
    bool linear_model_;        /** < If the variable is set to true, it is assumed that the system is linear; otherwise, non-linear.*/
    
    
    

     /**
     * Propagates the state according to the provided time interval. This function can be used to 
     * propagate the model. It might also be used in the RANSAC algorithm.
     * @param[in] dt  The amount of time the state needs to be propagated. The value of dt can be positive or negative. 
     *                a positive value would indicate forward propagation and a negative value would indicate backward propagation.
     * @return Returns the propagated state.
     */ 
    virtual Eigen::MatrixXd PropagateState(const Eigen::MatrixXd& state, const double dt)=0;


    /**
     * A state transition function propagates a state from one point in time to another. In the special 
     * case that the model is linear, the state transition function is a linear function. This means that we can represent
     * it using a matrix such that
     * \f]
     *      x_k = \Phi x_0 
     * \f]
     * Where \f$ x_0\f$ is the given state, \f$ x_k\f$ is the propagated state, and \f$ \Phi \f$ is the linear mapping. 
     * This function generates and returns the linear mapping \f$ \Phi \f$. The mapping is used with R-RANSAC.
     * This function is intended to be used only with linear models. If the user is implementing a non-linear model,
     * this function should not be implemented at the flag ModelBase::linear_model_ should be set to true.
     * @param[in] dt The amount of time the state needs to be propagated. The value of dt can be positive or negative. 
     *                a positive value would indicate forward propagation and a negative value would indicate backward propagation.
     * @return The linear map. 
     */ 
    virtual Eigen::MatrixXd GetLinearStateTransitionMatrix(const double dt)=0;

    /**
     * Propagates the model and error covariance to the current time.
     * @param[in] dt  The amount of time the model needs to be propagated.
     */ 
    virtual void PropagateModel(const double dt)=0; 


    /**
     * Uses the newly associated measurements to update the state of the model 
     * and update the error covariance. 
     * @param[in] param Contains all of the user defined parameters including the sources.
     */ 
    virtual void UpdateModel(const Parameters& param)=0;

    /**
     * Calculates and returns the matrix \f$F\f$.
     * @param[in] dt  The amount of time the model needs to be propagated. The value of dt can be positive or negative. 
     *                a positive value would indicate forward propagation and a negative value would indicate backward propagation. 
     * @return Returns the matrix \f$F\f
     */ 
    virtual Eigen::MatrixXd GetLinSysMatState(const double dt)=0;


    /**
     * Calculates and returns the matrix \f$G\f$.
     * @param[in] dt  The amount of time the model needs to be propagated. The value of dt can be positive or negative. 
     *                a positive value would indicate forward propagation and a negative value would indicate backward propagation.
     * @return Returns the matrix \f$F\f
     */
    virtual Eigen::MatrixXd GetLinSysMatNoise(const double dt)=0;


    /**
     * Calculates the linearized observation matrix \f$H\f$ defined by the source.
     * Since \f$H\f$ is dependent on the model and source, this function might have to 
     * claculate a unique \f$H\f$ for each source.  
     * @param source_ID A unique identifier to identify the source. 
     * @return Returns the matrix \f$H\f$
     */ 
    virtual Eigen::MatrixXd GetLinObsMatState(const unsigned int source_ID)=0;


    /**
     * Calculates an estimated measurement given a state of the model and the source ID.
     * Since there can be multiple sources, the measurement can be different based on the source. 
     * @param source_ID A unique identifier to identify the source. 
     * @return Returns a measurement \f$ y \f$ based on the state and source.
     */ 
    virtual Eigen::MatrixXd GetEstMeas(const Eigen::MatrixXd& state, const unsigned int source_ID)=0

    /**
     * Using the transformation data provided by the user, this function transforms the model from the previous global frame to the current global frame.
     * The previous global frame is the model is expressed in before the transformation is applied.
     * @param[in] T The transformation data provided by the user. The data should contain the information needed to transform the model.
     * @param[in] dt The time interval between the previous global frame and the current global frame. 
     */ 
    virtual void TransformModel(const Transformation& T, const double dt)=0;

};

//////////////////////////////////////////////////
//          model_lcv.h
//////////////////////////////////////////////////

#include "rransac/common/models/model_base.h"

// Of course add comments.

class ModelLCV : ModelBase{

public:

    ModelLCV();
    ~ModelLCV();

    Eigen::MatrixXd PropagateState(const Eigen::MatrixXd& state, const double dt);

    Eigen::MatrixXd GetLinearStateTransitionMatrix(const double dt);

    void PropagateModel(const double dt, const System& sys);

    void UpdateModel(const Parameters& param);

    Eigen::MatrixXd GetLinSysMatState(const double dt, const System& sys);

    Eigen::MatrixXd GetLinSysMatNoise(const double dt, const System& sys);

    void TransformModel(const Transformation& T, const double dt);

    Eigen::MatrixXd GetLinObsMatState(const unsigned int source_ID);

    Eigen::MatrixXd GetEstMeas(const Eigen::MatrixXd& state, const unsigned int source_ID);


};


//////////////////////////////////////////////////
//          model_factory.h
//////////////////////////////////////////////////

#include "rransac/common/models/model_lcv.h"

/** \class ModelFactory
 * Creates a new instant of a derived source class.
*/

class ModelFactory 
{
public:

    /**
     * Returns a unique pointer to a derived model class.
     * @param[in] type The type of derived model class to use. 
     * @see DerivedModel
     */
    std::unique_ptr<SourceBase> CreateMeasurementClass(const DerivedModels type);
};


//////////////////////////////////////////////////
//          model_data_association_base.h
//////////////////////////////////////////////////

enum DerivedModelDataAssociation
{
    kModelDataAssociationPDAF=0 /**< Uses the PDAF to associate measurements to the model and to update the model. */
};

class ModelDataAssociation
{
    /**
     * Associates the new measurements to existing models and updates the models using the newly associated
     * measurements. If a measurement is associated, it must be removed from the std::list of new measurements.
     * We supply this function full access to the system object. Terrible things could happen if this function
     * is not implemented properly.
     * @param sys The system object contains the models and new measurements.
     * @param dt The time elapsed since the last update.
     */ 
    virtual AssociateMeasurements(System& sys, const double dt)=0;
};



//////////////////////////////////////////////////
//          data_manager.h
//////////////////////////////////////////////////

class DataManager
{

public:


    /**
     * Removes expired measurements from the data tree, consensus set and clusters.
     * If removing expired measurements from the cluster causes the cluster to not have enough
     * measurements to be a cluster, the cluster with it measurements are deleted. We delete the measurements 
     * since all of their time stamps are close to being expired as well. 
     * @param sys The system object contains the data tree, models and clusters
     */  
    void RemoveExpiredMeasurements(System& sys);

    /**
     * Transforms all measurements and models from the previous global frame to the current global frame.
     * @param sys The system object contains the models
     */ 
    void TransformMeasAndModels(System& sys, const Transformation T, const double dt);

    /**
     * Associates the new measurements to existing models and updates the models using the newly associated
     * measurements. If a measurement is associated, it must be removed from the std::list of new measurements.
     * In order to associate new measurements, it will use a data association class.
     * @param sys The system object contains the models and new measurements.
     * @param dt The time elapsed since the last update.
     */ 
    virtual AssociateNewMeasurementsToModels(System& sys, const double dt);


    /**
     * After new measurements are associated with existing models, the remaining new measurements are associated
     * with clusters. They are simply added to the cluster.
     * @param sys The system object contains the clusters and remaining new measurements.
     */ 
    void AssociateRemainingNewMeasurementsToClusters(System& sys);

    /**
     * After new measurements are associated with existing models and existing clusters, the remaining new measurements are 
     * used to seed new clusters. 
     * @param sys The system object contains the clusters and remaining new measurements.
     */ 
    void GenereateNewClusters(System& sys);


    /**
     * Creates a data association object based on the mehtod given.
     * @param method The data association method that should be used.  @see DerivedDataAssociation
     */ 
    void CreateDataAssociationObject(const DerivedDataAssociation method);

    std::unique_ptr<DataAssociationBase> data_assoc_; /** < A pointer to the data association class to be used to associate new measurements with a model. */

};

//////////////////////////////////////////////////
//          model_initializer.h
//////////////////////////////////////////////////


/**
 * \class ModelInitializer
 * Uses the RANSAC algorithm to initialize a new model. When the model initializer is ran, it will
 * try to create a new model from every cluster. 
 */ 
class ModelInitializer
{
public:


    /**
     * Depending on the model and whether or not the time steps are fixed, RANSAC can 
     * optimize its performance by initializing variables that will be frequently used. 
     * For example, given a linear system with fixed time steps, we can compute before hand
     * the state transition matrices for the different times steps, the linearized system matrix
     * with respect to noise, and the linearized observation matrix with respect to the system state. 
     * @param System Contains the models and parameters.
     */ 
    void InitializeRANSAC(const System& sys);

    /**
     * For every cluster, RANSAC generates many different model hypotheses, constructs the consensus
     * set for each model hypotheses. The consensus set only contains pointers to the measurements which
     * are in the cluster. The model hypotheses with the largest consensus set (provided that the consensus
     * set has a cardinality greater than some threshold) is kept. The measurements', associated to best model hypotheses, flag Meas::associated
     * is set to true using the consensus set. The measurements are then removed from the consensus set and added to the 
     * model hypotheses. 
     */ 
    void GenerateModelHypotheses(const System& sys);

    /**
     * This function generates new models from the model hypotheses by incorporated all of the measurements
     * in its consensus set into the model using the Kalman Smoothing technique. 
     */ 
    void GenerateModels(System& sys);


    /**
     * When creating a new model from a cluster, associated measurements are removed from the cluster. We must
     * verify that the remaining measurements in the cluster still form a cluster. If it still forms a cluster,
     * then the cluster remains. Otherwise,  if the left over measurements are 
     * older measurements, the cluster and measurements are deleted since they will not be incorporated into any future model, 
     * and if the left over measurements are newer, they are added back onto the tree because they could be incorporated into a new cluster. 
     */ 
    void ManageClusters(System& sys);

private:

    std::vector<std::unique_ptr<ModelBase>> best_model_hypotheses_; /**< During the generation of model hypotheses, the best hypothsis from from each cluster will
                                                                         be added to this data structure provided that the hypothesis is good enough. Once a model 
                                                                        hypothesis is becomes a model, the model hypothesis is removed from this vector.*/
    
};



//////////////////////////////////////////////////
//          merge_models_base.h
//////////////////////////////////////////////////






//////////////////////////////////////////////////
//          model_manager.h
//////////////////////////////////////////////////

#include

class