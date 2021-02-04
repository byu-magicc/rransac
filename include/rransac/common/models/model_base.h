#ifndef RRANSAC_COMMON_MODELS_BASE_H_
#define RRANSAC_COMMON_MODELS_BASE_H_
#pragma once



#include <Eigen/Dense>
#include <vector>

#include "state.h"
#include "common/measurement/measurement_base.h"
#include "data_containers/consensus_set.h"
#include "parameters.h"
#include "common/sources/source_base.h"
#include "common/sources/source_RN.h"
#include "common/sources/source_SEN_pos_vel.h"
#include "common/sources/source_SEN_pose_twist.h"

namespace rransac {


struct ModelLikelihoodUpdateInfo{

bool in_local_surveillance_region;        /** < Indicates if the track is in the local surveillance region of the source */
int num_assoc_meas;                       /** < Indicates the number of measurements associated to the track from the source */
unsigned int source_index;                /** < The source index */
double volume;                            /** < The volume of the validation region */



};

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
template <typename tSource, typename tTransformation, int tCovDim,  typename tDerived> 
class ModelBase
{
public:
typedef typename tSource::State State;
typedef typename State::DataType DataType;
typedef tSource Source;
typedef tTransformation Transformation;
typedef tDerived Derived;
typedef Eigen::Matrix<DataType, Eigen::Dynamic, Eigen::Dynamic> MatXd;

static constexpr unsigned int g_dim_ = State::Group::dim_;
// static constexpr unsigned int cov_dim_ = tCovDim;
static constexpr unsigned int cov_dim_ = tCovDim;
typedef Eigen::Matrix<DataType,cov_dim_,cov_dim_> Mat;
typedef Eigen::Matrix<DataType,cov_dim_,1> VecCov;





public:
    State state_;                     /** < The estimated state of the phenomenon or target.*/
    Mat err_cov_;                     /** < The error covariance. */

    std::vector<std::vector<Meas<DataType>>> new_assoc_meas_;   /** < Measurements recently associated with the model. These measurements have not been used to update 
                                                     the model. Once an associated measurement has been used to update the model, they are added to the 
                                                     consensus set and removed from this vector. These measurements should have weights assigned to them. 
                                                     Each vector of measurements corresponds to a unique source ID. */

    ConsensusSet<Meas<DataType>> cs_;           /** < The consensus set. */
    
    double missed_detection_time_;  /** The time elapsed since a measurement was associated with the target. */
    
    long int label_;       /** When the model becomes a good model, it receives a unique label */   

    double model_likelihood_;         /** < The likelihood that the model represents an actual phenomenon. This value will be between 0 and 1. */
    std::vector<ModelLikelihoodUpdateInfo> model_likelihood_update_info_; /**< Contains the information needed to update the model_likelihood */

    // Useful matrices for propagating
    Mat F_;                           /** < The Jacobian of the state transition function w.r.t. the states */
    Mat G_;                           /** < The Jacobian of the state transition function w.r.t. the noise  */
    Mat Q_;                           /** < Process noise covariance */

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
     * @param[in] param Contains all of the user defined parameters.
     */ 
    void UpdateModel(const std::vector<Source>& sources, const Parameters& params) {
        if (new_assoc_meas_.size() > 0) {
            OPlusEQ(GetStateUpdate(sources, params));
            for (auto& new_measurements: new_assoc_meas_) {
                cs_.AddMeasurementsToConsensusSet(new_measurements);
            }
            missed_detection_time_ = 0;      // reset the missed detection time since a measurement was received
        }
        
        new_assoc_meas_.clear();
        UpdateModelLikelihood(sources);

    }

    void OPlusEQ(const VecCov& state_update){
        static_cast<tDerived*>(this)->DerivedOPlusEq(state_update);
    }


    /**
     * Calculates the Jacobian of the observation matrix with respect to the state estimate
     * evaluated at the current state estimate.
     * The Jacobian  is dependent on the measurement source.
     * @param source_ID A unique identifier to identify the source. 
     * @return Returns the Jacobian \f$H_k\f$
     */ 
    static MatXd GetLinObsMatState(const std::vector<Source>& sources, const State& state, const unsigned int source_ID){
        return sources[source_ID].GetLinObsMatState(state);
    }

    /**
     * Calculates the Jacobian of the observation matrix with respect to the measurement noise
     * evaluated at the current state estimate.
     * The Jacobian  is dependent on the measurement source.
     * @param source_ID A unique identifier to identify the source. 
     * @return Returns the Jacobian \f$V_k\f$
     */ 
    static MatXd GetLinObsMatSensorNoise(const std::vector<Source>& sources, const State& state, const unsigned int source_ID){
        return sources[source_ID].GetLinObsMatSensorNoise(state);
    }

    /**
     * Calculates an estimated measurement given a state of the model and the source ID.
     * Since there can be multiple sources, the measurement can be different based on the source. 
     * @param source_ID A unique identifier to identify the source. 
     * @return Returns a measurement \f$ y \f$ based on the state and source.
     */ 
    static MatXd GetEstMeas(const std::vector<Source>& sources,const State& state, const unsigned int source_ID){
        return sources[source_ID].GetEstMeas(state);
    }

    /**
     * Returns the innovation covariance associated with the measurement
     * @param m The measurement
     */ 
    MatXd GetInnovationCovariance(const std::vector<Source>& sources, const unsigned int source_index) const ;

    /**
     * Using the transformation data provided by the user, this function transforms the state estimate and error covariance
     * from the previous global frame to the current global frame.
     * @param[in] T The transformation object provided by the user. The object should already have the data it needs to transform the model.
     * @param[in] dt The time interval between the previous global frame and the current global frame. 
     */ 
    void TransformModel(const Transformation& T){
        if (!T.transform_null_)
            T.TransformTrack(state_,this->err_cov_);
    }

    void TransformConsensusSet(const Transformation& T) {
        if (!T.transform_null_)
            cs_.TransformConsensusSet(T);
    }

    /**
     * Computes the OMinus operation for the state
     */ 
    static VecCov OMinus(const tDerived& model1, const tDerived& model2) {
        return tDerived::DerivedOMinus(model1, model2);
    }

    /**
     * Update the model likelihood using the model_likelihood_update_info_
     */
    void UpdateModelLikelihood(const std::vector<Source>& sources); 

    /**
     * Returns a Random State
     */ 
    static State GetRandomState(){
        return tDerived::DerivedGetRandomState();
    }

    /**
     * Removes all of the measurements past the expiration time. 
     */
    void PruneConsensusSet(const DataType expiration_time) {
        cs_.PruneConsensusSet(expiration_time);
    } 


    /**
     * Adds a new measurement to the new_assoc_meas_ according to the source index.
     */ 
    void AddNewMeasurement( const Meas<DataType>& meas);

private:

void GetInnovationAndCovariance(const std::vector<Source>& sources, const std::vector<Meas<DataType>>& meas, Eigen::Matrix<DataType,cov_dim_,1>& state_update, Mat& cov_update);

Eigen::Matrix<DataType,tCovDim,1> GetStateUpdate(const std::vector<Source>& sources, const Parameters& params);


};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                            Definitions
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename tSource, typename tTransformation, int tCovDim,  typename tDerived>  
void ModelBase<tSource, tTransformation, tCovDim, tDerived>::Init(const Parameters& params) {
    err_cov_.setIdentity();
    F_.setIdentity();
    G_.setIdentity();
    SetParameters(params);
    model_likelihood_ = 0;
    missed_detection_time_ = 0;
    label_ = -1;                    // Indicates that it has not received a proper label.
}

//-------------------------------------------------------------------------------------------------------------------

template <typename tSource, typename tTransformation, int tCovDim,  typename tDerived> 
void ModelBase<tSource, tTransformation, tCovDim, tDerived>::PropagateModel(const DataType dt) {

    // Construct matrices to transform covariance.
    F_ = GetLinTransFuncMatState(state_,dt);
    G_ = GetLinTransFuncMatNoise(state_,dt);

    missed_detection_time_+=dt;


    // Transform covariance
    err_cov_ = F_*err_cov_*F_.transpose() + G_*Q_*G_.transpose();

    // Propagate state
    state_.g_.OPlusEq(state_.u_.data_*dt);

}

//-------------------------------------------------------------------------------------------------------------------

template <typename tSource, typename tTransformation, int tCovDim,  typename tDerived>  
Eigen::Matrix<typename tSource::State::DataType, tCovDim,1> ModelBase<tSource, tTransformation, tCovDim, tDerived>::GetStateUpdate(const std::vector<Source>& sources, const Parameters& params) {

Eigen::Matrix<DataType,cov_dim_,1> state_update_sum;
Eigen::Matrix<DataType,cov_dim_,1> state_update;
Eigen::Matrix<DataType,cov_dim_,1> update;
Mat cov;
Mat cov_sum;
Mat error_cov_inverse = err_cov_.inverse();
state_update_sum.setZero();
cov_sum.setZero();

// std::cout << "state: " << std::endl << state_.g_.data_ << std::endl;

// loop through the measurements per source
for (std::vector<Meas<DataType>> meas : new_assoc_meas_) {

        GetInnovationAndCovariance(sources, meas, state_update, cov);

    // std::cout << "state_update: " << std::endl << state_update << std::endl;
    // std::cout << "cov: " << std::endl << cov << std::endl;


    state_update_sum+= state_update;
    cov_sum += (cov.inverse() - error_cov_inverse);
    // cov_sum += cov;
}



// Update the error covariance
error_cov_inverse += cov_sum;
err_cov_ = error_cov_inverse.inverse();
// err_cov_ = cov;

#ifdef DEBUG_BUILD
     if (err_cov_.determinant() <=0 )
        throw std::runtime_error("ModelBase::GetStateUpdate The determinant of the error covariance is <=0. It must be positive since it is a positive definite matrix. ") ;
#endif

// Update the state
update = err_cov_ * state_update_sum;
// update = cov * state_update_sum;

return update;

}

//---------------------------------------------------------------------------------------------------------

template <typename tSource, typename tTransformation, int tCovDim,  typename tDerived>  
void ModelBase<tSource, tTransformation, tCovDim, tDerived>::GetInnovationAndCovariance(const std::vector<tSource>& sources, const std::vector<Meas<DataType>>& meas, Eigen::Matrix<DataType,cov_dim_,1>& state_update, Mat& cov_update) {

state_update.setZero();
cov_update = err_cov_;


Eigen::MatrixXd H = GetLinObsMatState(sources, state_,meas.front().source_index);              // Jacobian of observation function w.r.t. state
Eigen::MatrixXd V = GetLinObsMatSensorNoise(sources, state_,meas.front().source_index);        // Jacobian of observation function w.r.t. noise
Eigen::MatrixXd K;                                                                    // Kalman Gain
Eigen::MatrixXd S_inverse;                                                            // Innovation covariance inverse
Eigen::MatrixXd nu_i;                                                                 
Eigen::MatrixXd nu(V.rows(),1);                                                       // Total innovation term
nu.setZero();
Eigen::MatrixXd covSum(V.rows(),V.rows());
covSum.setZero();
// std::cerr << "H: " << std::endl << H << std::endl << std::endl;
// std::cerr << "V: " << std::endl << V << std::endl << std::endl;
Meas<DataType> estimated_meas = sources[meas.front().source_index].GetEstMeas(state_);

// std::cerr << "tmp: " << std::endl << tmp << std::endl << std::endl;
// std::cerr << "H*err_cov_*H.transpose() : " << std::endl << H*err_cov_*H.transpose()  << std::endl << std::endl;


S_inverse = (H*err_cov_*H.transpose() + V*sources[meas.front().source_index].params_.meas_cov_*V.transpose()).inverse();
K = err_cov_*H.transpose()*S_inverse;

// std::cerr << "K: " << std::endl << K << std::endl << std::endl;


DataType B0 = 1;

// Get total weighted innovation and part of the cov_tilde
for (Meas<DataType> m : meas) {
    nu_i = sources[m.source_index].OMinus(m, estimated_meas);
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

// cov_update +=  K*H*err_cov_; 
// cov_update = H.transpose()*tmp*H;
state_update = H.transpose()*S_inverse*nu;


}

//---------------------------------------------------------------------------------------------------------

template <typename tSource, typename tTransformation, int tCovDim,  typename tDerived>  
void ModelBase<tSource, tTransformation, tCovDim, tDerived>::UpdateModelLikelihood(const std::vector<tSource>& sources) {
    for (auto& update_info : model_likelihood_update_info_) {
        if ( update_info.in_local_surveillance_region) {
            const tSource& source = sources[update_info.source_index];
            model_likelihood_ += std::log(1 + source.params_.probability_of_detection_*source.params_.gate_probability_*( update_info.num_assoc_meas/(source.params_.spacial_density_of_false_meas_*update_info.volume)-1));
        }
    }
    model_likelihood_update_info_.clear();
} 

//---------------------------------------------------------------------------------------------------------


template <typename tSource, typename tTransformation, int tCovDim,  typename tDerived>  
Eigen::Matrix<typename tSource::State::DataType,Eigen::Dynamic,Eigen::Dynamic> ModelBase<tSource, tTransformation, tCovDim, tDerived>::GetInnovationCovariance(const std::vector<Source>& sources, const unsigned int source_index) const {

    MatXd H = GetLinObsMatState(sources, this->state_,source_index);         // Jacobian of observation function w.r.t. state
    MatXd V = GetLinObsMatSensorNoise(sources, this->state_,source_index);   // Jacobian of observation function w.r.t. noise

    return H*err_cov_*H.transpose() + V*sources[source_index].params_.meas_cov_*V.transpose();
}

//---------------------------------------------------------------------------------------------------------

template <typename tSource, typename tTransformation, int tCovDim,  typename tDerived>
void ModelBase<tSource, tTransformation, tCovDim, tDerived>::AddNewMeasurement( const Meas<DataType>& meas) {

    bool meas_added = false;

#ifdef DEBUG_BUILD
    // The new associated measurements should all have the same time stamp.
    if(new_assoc_meas_.size() != 0)
        if(new_assoc_meas_.front().size() != 0)
            if(meas.time_stamp != new_assoc_meas_.front().front().time_stamp)
                throw std::runtime_error("ModelBase::AddNewMeasurement Not all new measurements have the same time stamp");

#endif

    // See if there is already a measurement with the same source index. If it is, add it to the list
    for(auto iter = this->new_assoc_meas_.begin(); iter != this->new_assoc_meas_.end(); ++iter) {
        if(iter->begin()->source_index == meas.source_index) {
            iter->push_back(meas);
            meas_added = true;
            break;
        }
    }

    // Create a new sub list and add the measurement
    if(!meas_added) {
        this->new_assoc_meas_.emplace_back(std::vector<Meas<DataType>>{meas});
    }

}



} // namespace rransac

#endif //RRANSAC_COMMON_MODELS_BASE_H_
