

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
#include "rransac/common/param.h"


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
//          source_base.h
//////////////////////////////////////////////////

#include <Eigen/Core>


/** \enum DerivedSources
 * Used to indicate what type of derived source class to create. 
*/
enum DerivedSources
{
    kSourceCOF=0 /**< This source detects point measurements using optical flow. */
}


/** \class SourceBase 
 * A source is an algorithm that takes data from a sensor and extracts measurements from it. Due to the nature of
 * the source, there is measurement noise and false measurements. The measurement noise is assumed to be sampled from
 * a white-noise, zero-mean, Gaussian distribution. If the measurement covariance is fixed, then it should be added to 
 * the source; otherwise, it is supplied with each measurement. We assume that the false measurements are uniformly
 * distributed in the local surveillance region, and that the number
 * of expected measurements per sensor scan from a source can be modeled using a Poisson distribution. 
 * 
 * Since there can be many different sources, we need to be able to distinguish the different sources using a unique ID.
 * When a measurement is received, the measurement indicates which source it came from using the source ID tag @see Meas.
 * We also need to know the dimension of the measurement and how to compute the linearized observation matrix that is used
 * in a Kalman update step. This is used to update the model. 
 * 
 */ 

class SourceBase
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

    /** Given the current state of the system, this function returns the linearized observation matrix. This
     * Is the same linearized observation matrix that is used in the Kalman update step. 
     * @param state The current state of the system. 
     * @return The function returns the linearized observation matrix.
     * @see ModelBase
    */
    virtual Eigen::MatrixXd GetLinObsMat(Eigen::MatrixXd state)=0;

    /**
     * When a source is added by the user, we need to verify that the variables were set correctly.
     * This function verifies that they member variables are correct.
     * @return Returns true if the member variables are set correctly; otherwise, returns false.
     */ 
    bool VerifyMemberVariables();
    
};

//////////////////////////////////////////////////
//          source_cof.h
//////////////////////////////////////////////////

#include "rransac/common/sources/source_base.h"

/**
 * \class SourceCOF
 * The name stands for source camera optical flow.
 * This is a derived class of SourceBase. The measurement source
 * uses optical flow to calculate the position and velocity of a moving target
 * in the normalized image plane.
 */ 

class SourceCOF : SourceBase
{

public:

    Eigen::MatrixXd GetLinObsMat(Eigen::MatrixXd state)=0
};

//////////////////////////////////////////////////
//          source_factory.h
//////////////////////////////////////////////////

#include "rransac/common/sources/source_cof.h"

/** \class SourceFactory
 * Creates a new instant of a derived source class.
*/

class SourceFactory 
{
public:

    /**
     * Returns a unique pointer to a derived source class.
     * @param[in] type The type of derived source class to use. 
     * @see DerivedSource
     */
    std::unique_ptr<SourceBase> CreateMeasurementClass(const DerivedSource type);
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

public:
    Eigen::MatrixXd state_;    /** < The estimated state of the phenomenon or target.*/
    Eigen::MatrixXd err_cov_;  /** < The error covariance. */
    double model_likelihood_;  /** < The likelihood that the model represents an actual phenomenon. This value will be between 0 and 1. */
    std::vector<Meas> new_assoc_meas_; /** < Measurements recently associated with the model. These measurements have not been used to update 
                                                     the model. Once an associated measurement has been used to update the model, they are added to the 
                                                     consensus set and removed from this vector.*/


    /**
     * Propagates the model and error covariance to the current time.
     * @param[in] dt  The amount of time the model needs to be propagated.
     * @param[in] sys 
     */ 
    virtual void PropagateModel(const double dt, const System& sys)=0; 


    /**
     * Uses the newly associated measurements to update the state of the model 
     * and update the error covariance. 
     * @param[in] param Contains all of the user defined parameters including the sources.
     */ 
    virtual void UpdateModel(const Parameters& param)=0;

    /**
     * Calculates and returns the matrix \f$F\f$.
     * @param[in] dt  The amount of time the model needs to be propagated.
     * @param[in] sys 
     */ 
    virtual Eigen::MatrixXd GetLinSysMatState(const double dt, const System& sys)=0;


    /**
     * Calculates and returns the matrix \f$G\f$.
     * @param[in] dt  The amount of time the model needs to be propagated.
     * @param[in] sys 
     */
    virtual Eigen::MatrixXd GetLinSysMatNoise(const double dt, const System& sys)=0;


private:

};

