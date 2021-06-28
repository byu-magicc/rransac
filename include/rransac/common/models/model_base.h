#ifndef RRANSAC_COMMON_MODELS_BASE_H_
#define RRANSAC_COMMON_MODELS_BASE_H_
#pragma once



#include <Eigen/Dense>
#include <vector>

#include "lie_groups/state.h"
#include "rransac/common/measurement/measurement_base.h"
#include "rransac/data_containers/consensus_set.h"
#include "rransac/parameters.h"
#include "rransac/common/sources/source_base.h"
#include "rransac/common/sources/source_RN.h"
#include "rransac/common/sources/source_SEN_pos_vel.h"
#include "rransac/common/sources/source_SEN_pose_twist.h"

namespace rransac {


/**
 * \class ModelLikelihoodUpdateInfo
 * Contains the information necessary to update the track's likelihood after data association. * 
 */ 
struct ModelLikelihoodUpdateInfo{

bool in_lsr_and_produced_meas;            /** < Indicates if the track is in the local surveillance region of the source and the source produced measurements. */
int num_assoc_meas;                       /** < Indicates the number of measurements associated to the track from the source */
double delta;                             /** < A value used in the integrated probabilistic data association filter.  */

void Reset() {
    in_lsr_and_produced_meas = false;
    num_assoc_meas = 0;
    delta = 0;
}

};

/**
 * \class ModelBase
 * 
 * R-RANSAC is designed to be modular and work with a variety of models. However, this 
 * version of R-RANSAC is extended to work with any Lie group. We assume that the dynamic model
 * is a discrete, time-invariant system of the form 
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
 * 
 * The class is designed to facilitate propagating and updating the model using a method similar to the Kalman filter update.
 * Methods similar to the Kalman filter include the EKF, UKF, IF, EIF, PDA, JPDA, etc. Reglardless of the method used, the
 * system model \f$f\f$ needs to be linearized with respect to the state and the noise. Also, the observation model \f$h\f$
 * needs to be linearized with respect to the state and measurement noise. 
 * 
 * In most methods, the following matrices will be helpful. 
 * Let \f$F\f$ denote the linearized system model with respect to the state calculated as 
 * \f[
 *     F=\frac{\partial f}{\partial x}\bigg |_{x_k,\delta_t},
 * \f]
 * \f$G\f$ denote the linearized system model with respect to the noise calculated as 
 * \f[
 *     G=\frac{\partial f}{\partial w}\bigg |_{x_k,\delta_t},
 * \f] 
 * \f$H\f$ denote the linearized observation model with respect to the noise calculated as 
 * \f[
 *     H=\frac{\partial h}{\partial x}\bigg |_{x_k,0}.
 * \f] and 
 * \f[
 *     V=\frac{\partial h}{\partial v}\bigg |_{x_k,0}.
 * \f] and 
 * 
 * 
 * This class uses the curiously reccuring template pattern desing methodology and as such specifies the API for all derived classes. 
 * The derived classes must be able to compute the Jacobians \f$F\f$ and \f$G\f$ as well as perform the OPlus operation on a state update, 
 * compute the OMinus operation on two states, and produce a random state. See ModelBase::GetLinTransFuncMatState(), ModelBase::GetLinTransFuncMatNoise(),
 * ModelBase::OPlusEQ(), ModelBase::OMinus(), and ModelBase::GetRandomState().
 * 
 * 
 * See the paper ___ for detailed information regarding the system model. 
 * 
 */ 

template <typename tSourceContainer, int tCovDim,  typename tDerived> 
class ModelBase
{
public:
    typedef typename tSourceContainer::State State;                             /**< The state of the target. @see State. */
    typedef typename State::DataType DataType;                                  /**< The scalar object for the data. Ex. float, double, etc. */
    typedef tSourceContainer SourceContainer;                                   /**< The object type of the source. @see SourceBase. */
    typedef typename tSourceContainer::Transformation Transformation;           /**< The object type of the measurement and track transformation. */
    typedef tDerived Derived;                                                   /**< The object type of the derived class. */
    typedef Eigen::Matrix<DataType, Eigen::Dynamic, Eigen::Dynamic> MatXd;      

    static constexpr unsigned int g_dim_ = State::Group::dim_;                  /**< The dimension of the pose of the state, i.e. the dimension of the group portion of the state. */
    static constexpr unsigned int cov_dim_ = tCovDim;                           /**< The dimension of the error covariance. */
    static constexpr unsigned int num_sources_ = tSourceContainer::num_sources_;/**< The number of sources in the source container. */
    typedef Eigen::Matrix<DataType,cov_dim_,cov_dim_> Mat;                      /**< The object type of the error covariance, Jacobians, and others. */
    typedef Eigen::Matrix<DataType,cov_dim_,1> VecCov;                          /**< The object type of state update. */


    std::vector<ModelLikelihoodUpdateInfo> model_likelihood_update_info_; /**< Contains the information needed to update the model_likelihood */
    std::vector<std::vector<Meas<DataType>>> new_assoc_meas_;             /**< Measurements recently associated with the track. These measurements have not been used to update 
                                                                                the track. Once an associated measurement has been used to update the track, they are added to the 
                                                                                consensus set and removed from this vector. These measurements should have weights assigned to them. 
                                                                                Each vector of measurements corresponds to a unique source ID. */

    State state_;                      /**< The estimated state of the target.*/
    Mat err_cov_;                      /**< The error covariance. */
    ConsensusSet<Meas<DataType>> cs_;  /**< The consensus set. */
    double newest_measurement_time_stamp;  /**< The time stamp associated with the newest measurement. */
    long int label_;                   /**< When the track becomes a good track, it receives a unique label. */   
    double model_likelihood_;          /**< The likelihood that the track represents an actual phenomenon.  */
    Mat F_;                            /**< The Jacobian of the state transition function w.r.t. the states */
    Mat G_;                            /**< The Jacobian of the state transition function w.r.t. the noise  */
    std::vector<MatXd> innovation_covariances_; /**< A vector of innovation covariances. The vector index corresponds to the measurement source index. It is used 
                                                    in order to mitigate how many times the innovation covariance is calculated. */
    std::vector<bool> innov_cov_set_;           /**< Indicates if the member variable innovation_covariances_ have been set and are valid since new measurements were received. This
                                                     variable is reset to false after the update step in UpdateModel and the propagation step in PropagateModel since
                                                     the innovation covariances are no longer valid. The index of the vector corresponds to the source index. */
    Mat Q_;                            /**< Process noise covariance */

#if RRANSAC_VIZ_HOOKS
    std::vector<Eigen::MatrixXd> S_validation_; /**< The innovation covariance used to compute the last validation region. */
#endif

    // Default constructor
    ModelBase()=default;

    // Default destructor
    ~ModelBase()=default;

     /**
     * Initializes the model.
     * @param[in] sources A reference to the sources contained in system
     * @param[in] params  The system parameters specified by the user
     */ 
    void Init(const Parameters& params);
    
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
    static State PropagateState(const State& state, const DataType dt) {
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
    static Mat GetLinTransFuncMatState(const State& state, const DataType dt) {
        return tDerived::DerivedGetLinTransFuncMatState(state, dt);        
    }

    /**
     * Computes the Jacobian of the state transition function with respect to the noise evaluated at the current state estimate.
     * @param[in] state The state to linearize about
     * @param[in] dt  A time interval
     * @return Returns the Jacobian \f$ G_k \f$
     */
    static Mat GetLinTransFuncMatNoise(const State& state, const DataType dt){
        return tDerived::DerivedGetLinTransFuncMatNoise(state, dt);

    }

    /**
     * Propagates the state estimate and error covariance to the current time.
     * @param[in] dt  The amount of time the model needs to be propagated.
     */ 
    void PropagateModel(const DataType dt);

    /**
     * Uses the newly associated measurements to update the state estimate, error covariance, model likelihood and consensus set using a 
     * centralized measurement fusion.
     * @param[in] source_container The container of all of the sources
     * @param[in] param Contains all of the user defined parameters.
     */ 
    void UpdateModel(const SourceContainer& source_container, const Parameters& params);


    /**
     * Performs the OPlus operation on the track's state using the state update provided.
     * @param[in] state_update The state update that will be added to the track's current state.
     */
    void OPlusEQ(const VecCov& state_update){ static_cast<tDerived*>(this)->DerivedOPlusEq(state_update);}


    /**
     * Calculates the Jacobian of the observation matrix with respect to the state estimate
     * evaluated at the current state estimate.
     * The Jacobian  is dependent on the measurement source.
     * @param[in] source_container The container of all of the sources
     * @param[in] source_index The index to the source whose parameters are to be changed.
     * @param[in] state A state of the target.
     * @param[in] transform_state A flag used to indicate if the state needs to be transformed 
     * @param[in] transform_data The data needed to transform the state
     * @return Returns the Jacobian \f$H_k\f$
     */ 
    static MatXd GetLinObsMatState(const SourceContainer& source_container,  const unsigned int source_index, const State& state, const bool transform_state, const MatXd& transform_data){
        return source_container.GetLinObsMatState(source_index,state,transform_state,transform_data);
    }

    /**
     * Calculates the Jacobian of the observation matrix with respect to the measurement noise
     * evaluated at the current state estimate.
     * The Jacobian  is dependent on the measurement source.
     * @param[in] source_container The container of all of the sources
     * @param[in] source_index The index to the source whose parameters are to be changed.
     * @param[in] state A state of the target.
     * @param[in] transform_state A flag used to indicate if the state needs to be transformed 
     * @param[in] transform_data The data needed to transform the state
     * @return Returns the Jacobian \f$V_k\f$
     */ 
    static MatXd GetLinObsMatSensorNoise(const SourceContainer& source_container,  const unsigned int source_index, const State& state, const bool transform_state, const MatXd& transform_data){
        return source_container.GetLinObsMatSensorNoise(source_index,state,transform_state,transform_data);
    }

    // /**
    //  * Calculates an estimated measurement given a state of the model and the source ID.
    //  * Since there can be multiple sources, the measurement can be different based on the source. 
    //  * @param source_ID A unique identifier to identify the source. 
    //  * @return Returns a measurement \f$ y \f$ based on the state and source.
    //  */ 
    // static MatXd GetEstMeas(const std::vector<Source>& sources,const State& state, const unsigned int source_ID){
    //     return sources[source_ID].GetEstMeas(state);
    // }

    /**
     * Returns the innovation covariance associated with a source
     * @param[in] source_container The container of all of the sources
     * @param[in] source_index The index of the source of which we want to compute the innovation covariance.
     * @param[in] transform_state A flag used to indicate if the state needs to be transformed 
     * @param[in] transform_data The data needed to transform the state
     */ 
    MatXd GetInnovationCovariance(const SourceContainer& source_container, const unsigned int source_index, const bool transform_state, const MatXd& transform_data) const ;

    /**
     * Returns the innovation covariance associated with a source. If the innovation covariance has been set in innovation_covariances_
     * then it returns the precomputed value; otherwise, it computes the innovation covariance from the constant method GetInnovationCovariance
     * stores the result in innovation_covariances_ and returns it. 
     * @param[in] source_container The container of all of the sources
     * @param[in] source_index The index of the source of which we want to compute the innovation covariance.
     * @param[in] transform_state A flag used to indicate if the state needs to be transformed 
     * @param[in] transform_data The data needed to transform the state
     */ 
    MatXd GetInnovationCovariance(const SourceContainer& source_container, const unsigned int source_index, const bool transform_state, const MatXd& transform_data);

    /**
     * Using the transformation data provided by the user, this function transforms the state estimate and error covariance
     * from the previous tracking frame to the current tracking frame.
     * @param[in] T The transformation object provided by the user. The object should already have the data it needs to transform the model.
     * @param[in] dt The time interval between the previous tracking frame and the current tracking frame. 
     */ 
    void TransformModel(const Transformation& T){ T.TransformTrack(state_,this->err_cov_); }

    /**
     * Using the transformation data provided by the user, this function transforms the consensus set.
     * @param[in] T The transformation object provided by the user. The object should already have the data it needs to transform the model.
     */  
    void TransformConsensusSet(const Transformation& T) { cs_.TransformConsensusSet(T); }

    /**
     * Computes the OMinus operation for the state. In other words, it computes the geodesic distance between the two track's state estimate.
     * This OMinus operation can be different than the operation defined by the state. 
     * @param[in] track1 A track of this type.
     * @param[in] track2 A track of this type.
     * @return Returns the geodesic distance between the track's state estimate.
     */ 
    static VecCov OMinus(const tDerived& track1, const tDerived& track2) {
        return tDerived::DerivedOMinus(track1, track2);
    }

    /**
     * Returns a Random State
     */ 
    static State GetRandomState(){
        return tDerived::DerivedGetRandomState();
    }

    /**
     * Removes all of the measurements past the expiration time. 
     * @param expiration_time The expiration time of the measurements.
     */
    void PruneConsensusSet(const DataType expiration_time) {
        cs_.PruneConsensusSet(expiration_time);
    } 


    /**
     * Adds a new measurement to the ModelBase::new_assoc_meas_ according to the source index.
     * @param[in] meas The measurement to be added. 
     */ 
    void AddNewMeasurement( const Meas<DataType>& meas);

private:

/**
 * Calculates the state update and covariance update for a single source. The state update and covariance update are given to PerformCentralizedMeasurementFusion in order
 * to perform centralized measurement fusion.
 * @param[in] source_container The container of all of the sources
 * @param[in] meas          All of the new measurements produced by a single source.
 * @param[out] state_update The calculated state update to be applied to the current state estimate.
 * @param[out] cov_update The calculated covariance update to be applied to the current error covariance.
 */ 
void GetStateUpdateAndCovariance(const SourceContainer& source_container, const std::vector<Meas<DataType>>& meas, Eigen::Matrix<DataType,cov_dim_,1>& state_update, Mat& cov_update);

/**
 * Since there can be multiple sources providing measurements, we fuse the measurements together using a centralized measurement fusion method discussed in Tracking and Data Fusion by Bar-Shalom 2011.
 * @param[in] params  The system parameters.
 * @param[in] source_container The container of all of the sources
 * @return The state update. 
 */ 
Eigen::Matrix<DataType,tCovDim,1> PerformCentralizedMeasurementFusion(const SourceContainer& source_container, const Parameters& params);


};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                            Definitions
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename tSourceContainer, int tCovDim,  typename tDerived> 
void ModelBase<tSourceContainer, tCovDim, tDerived>::Init(const Parameters& params) {

    if (params.set_initial_error_covariance_to_id_) {
        err_cov_.setIdentity();
    } else {
        if (err_cov_.rows() != params.initial_error_covariance_.rows()) {
            throw std::runtime_error("ModelBase::Init The initial error covariance is not the right dimension. Its dimension is " + std::to_string(params.initial_error_covariance_.rows()) + ". It should be " + std::to_string(err_cov_.rows()));

        } else {
            err_cov_ = params.initial_error_covariance_;
        }
    }
    
    F_.setIdentity();
    G_.setIdentity();
    SetParameters(params);
    newest_measurement_time_stamp=0;
    model_likelihood_ = 0.5;
    label_ = -1;                    // Indicates that it has not received a proper label.
    innov_cov_set_.resize(num_sources_,false);
    innovation_covariances_.resize(num_sources_);
    new_assoc_meas_.resize(num_sources_);
    model_likelihood_update_info_.resize(num_sources_);
}

//-------------------------------------------------------------------------------------------------------------------

template <typename tSourceContainer, int tCovDim,  typename tDerived> 
void ModelBase<tSourceContainer, tCovDim, tDerived>::PropagateModel(const DataType dt) {

    // Construct matrices to transform covariance.
    F_ = GetLinTransFuncMatState(state_,dt);
    G_ = GetLinTransFuncMatNoise(state_,dt);

    


    // Transform covariance
    err_cov_ = F_*err_cov_*F_.transpose() + G_*Q_*dt*G_.transpose();

    // Propagate state
    state_.g_.OPlusEq(state_.u_.data_*dt);
    std::fill(innov_cov_set_.begin(), innov_cov_set_.end(),false);

}

//-------------------------------------------------------------------------------------------------------------------

template <typename tSourceContainer, int tCovDim,  typename tDerived> 
void ModelBase<tSourceContainer, tCovDim, tDerived>::UpdateModel(const SourceContainer& source_container, const Parameters& params) {

    for(auto& source_meas : new_assoc_meas_) {
        if(source_meas.size() > 0) {
            newest_measurement_time_stamp = source_meas.front().time_stamp;
            break;
        }
    }

    OPlusEQ(PerformCentralizedMeasurementFusion(source_container, params));
    for (auto& new_measurements: new_assoc_meas_) {
        cs_.AddMeasurementsToConsensusSet(new_measurements);
    }
    
    
    
    // UpdateModelLikelihood(sources);
    // Reset member variables
    std::fill(innov_cov_set_.begin(), innov_cov_set_.end(),false);
    for (auto iter = new_assoc_meas_.begin(); iter != new_assoc_meas_.end(); ++iter) {
        iter->clear();
    }

}

//-------------------------------------------------------------------------------------------------------------------

template <typename tSourceContainer, int tCovDim,  typename tDerived> 
Eigen::Matrix<typename tSourceContainer::State::DataType, tCovDim,1> ModelBase<tSourceContainer, tCovDim, tDerived>::PerformCentralizedMeasurementFusion(const tSourceContainer& source_container, const Parameters& params) {

Eigen::Matrix<DataType,cov_dim_,1> state_update_sum;
Eigen::Matrix<DataType,cov_dim_,1> state_update;
Eigen::Matrix<DataType,cov_dim_,1> update;
Mat cov;
Mat cov_sum;
Mat error_cov_inverse = err_cov_.inverse();
state_update_sum.setZero();
cov_sum.setZero();

// loop through the measurements per source
for (auto& meas : new_assoc_meas_) {
    if(meas.size()>0) {
        GetStateUpdateAndCovariance(source_container, meas, state_update, cov);

        state_update_sum+= state_update;
        cov_sum += (cov.inverse() - error_cov_inverse);
    }
}



// Update the error covariance
error_cov_inverse += cov_sum;
err_cov_ = error_cov_inverse.inverse();

#ifdef DEBUG_BUILD
     if (err_cov_.determinant() <=0 )
        throw std::runtime_error("ModelBase::GetStateUpdate The determinant of the error covariance is <=0. It must be positive since it is a positive definite matrix. ") ;
#endif

// Update the state
update = state_update_sum;


return update;

}

//---------------------------------------------------------------------------------------------------------

template <typename tSourceContainer, int tCovDim,  typename tDerived>   
void ModelBase<tSourceContainer, tCovDim, tDerived>::GetStateUpdateAndCovariance(const tSourceContainer& source_container, const std::vector<Meas<DataType>>& meas, Eigen::Matrix<DataType,cov_dim_,1>& state_update, Mat& cov_update) {

state_update.setZero();


Eigen::MatrixXd H = source_container.GetLinObsMatState(meas.front().source_index, state_, meas.front().transform_state, meas.front().transform_data_t_m);      // Jacobian of observation function w.r.t. state
Eigen::MatrixXd V = source_container.GetLinObsMatSensorNoise(meas.front().source_index, state_, meas.front().transform_state, meas.front().transform_data_t_m);                             // Jacobian of observation function w.r.t. noise
Eigen::MatrixXd K;                                                                    // Kalman Gain
Eigen::MatrixXd S_inverse;                                                            // Innovation covariance inverse
Eigen::MatrixXd nu_i;                                                                 
Eigen::MatrixXd nu(V.rows(),1);                                                       // Total innovation term
nu.setZero();
Eigen::MatrixXd covSum(V.rows(),V.rows());
covSum.setZero();

Meas<DataType> estimated_meas = source_container.GetEstMeas(meas.front().source_index, state_, meas.front().transform_state, meas.front().transform_data_t_m); 

S_inverse = GetInnovationCovariance(source_container,meas.front().source_index, meas.front().transform_state, meas.front().transform_data_t_m).inverse();
K = err_cov_*H.transpose()*S_inverse;

DataType B0 = 1;

// Get total weighted innovation and part of the cov_tilde
for (Meas<DataType> m : meas) {
    nu_i = source_container.OMinus(m.source_index,m,estimated_meas);
    nu += m.weight*nu_i;
    covSum+= m.weight*nu_i*nu_i.transpose();
    B0 -= m.weight;
}

// Finish constructing cov_sum
covSum -= nu*nu.transpose();


// construct covariance. It has been verified
cov_update = err_cov_+ K*(covSum*K.transpose() -(1-B0)*H*err_cov_); 

state_update = K*nu;



}

//---------------------------------------------------------------------------------------------------------


template <typename tSourceContainer, int tCovDim,  typename tDerived>  
Eigen::Matrix<typename tSourceContainer::State::DataType,Eigen::Dynamic,Eigen::Dynamic> ModelBase<tSourceContainer, tCovDim, tDerived>::GetInnovationCovariance(const SourceContainer& source_container, const unsigned int source_index, const bool transform_state, const MatXd& transform_data) const {

    // If the innovation covariance has been set this sensor scan return it. 
    if(innov_cov_set_[source_index]) {
       return innovation_covariances_[source_index]; 
    }

    MatXd H = source_container.GetLinObsMatState(source_index, this->state_, transform_state, transform_data);         // Jacobian of observation function w.r.t. state
    MatXd V = source_container.GetLinObsMatSensorNoise(source_index, this->state_, transform_state, transform_data);   // Jacobian of observation function w.r.t. noise

    return H*err_cov_*H.transpose() + V*source_container.GetParams(source_index).meas_cov_*V.transpose();
}

//---------------------------------------------------------------------------------------------------------


template <typename tSourceContainer, int tCovDim,  typename tDerived>  
Eigen::Matrix<typename tSourceContainer::State::DataType,Eigen::Dynamic,Eigen::Dynamic> ModelBase<tSourceContainer, tCovDim, tDerived>::GetInnovationCovariance(const SourceContainer& source_container, const unsigned int source_index, const bool transform_state, const MatXd& transform_data) {


    // If the innovation covariance has not been set this sensor scan, calculate it.
    if(!innov_cov_set_[source_index]) {

        innovation_covariances_[source_index] = const_cast<const ModelBase*>(this)->GetInnovationCovariance(source_container,source_index,transform_state,transform_data);
        innov_cov_set_[source_index] = true;
    }

    return innovation_covariances_[source_index];

}

//---------------------------------------------------------------------------------------------------------

template <typename tSourceContainer, int tCovDim,  typename tDerived> 
void ModelBase<tSourceContainer, tCovDim, tDerived>::AddNewMeasurement( const Meas<DataType>& meas) {

    

#ifdef DEBUG_BUILD
    // The new associated measurements should all have the same time stamp.
    if(new_assoc_meas_.size() != 0)
        if(new_assoc_meas_.front().size() != 0)
            if(meas.time_stamp != new_assoc_meas_.front().front().time_stamp)
                throw std::runtime_error("ModelBase::AddNewMeasurement Not all new measurements have the same time stamp");

#endif

    new_assoc_meas_[meas.source_index].push_back(meas);

}



} // namespace rransac

#endif //RRANSAC_COMMON_MODELS_BASE_H_
