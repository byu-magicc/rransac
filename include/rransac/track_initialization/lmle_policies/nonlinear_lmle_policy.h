#ifndef RRANSAC_TRACK_INITIALIZATION_LMLE_POLICIES_NONLINEAR_LMLE_POLICY_H_
#define RRANSAC_TRACK_INITIALIZATION_LMLE_POLICIES_NONLINEAR_LMLE_POLICY_H_


#include "data_containers/cluster.h"
#include "system.h"
#include "common/measurement/measurement_base.h"
#include <Eigen/Dense>


namespace rransac
{


/** \class NonLinearLMLEPolicy
 * This policy is used for Nonlinear Time Invariant Models. Using a set of measurements, by which the 
 * model is at least locally observable, it produces a current state estimate. The templates for the class
 * are tModel and tSeed. tModel is the model type and tSeed is a policy to seed
 * the LMLE optimization. 
 */ 

template<typename tModel, template<typename > typename tSeed>    
class NonLinearLMLEPolicy {

public:

typedef typename tModel::State State;
typedef typename tModel::Source Source;
typedef tModel Model;

/**
 * Generates a hypothetical state at the current time step using the provided measurements in meas_subset.
 * @param meas_subset The container of iterators to measurements that will be used to estimate the hypothetical state
 * @param curr_time The current time
 * @param sources The vector of sources used. 
 */ 
static State GenerateStateEstimatePolicy(const std::vector<Cluster::IteratorPair>& meas_subset, const System<tModel>& sys);


template<typename tModel>
struct CostFunctor {

    CostFunctor(Meas m, double dt, const System<tModel>& sys) : m_(m), dt_(dt), sys_(sys) {
        dt_ = m.time_stamp - sys_.current_time_;
        src_index_ = m_.source_index;
    }

    /**
     * The components of x is the Lie algebra of the state
     * 
     */ 
    template <typename T>
    bool operator(const T* const x , const T* const y, T* r) const {

    // Convert array to Eigen vector
    for (int ii = 0; ii < tModel::g_dim_*2; ++ii)
        x_vector_(ii,0) = x[ii];

    // Convert to state, propagate state to the time step of the measurement and get the error
    // between the estimated measurement and the measurement
    state_ = exp(x_vector_);
    state_ = tModel::PropagateState(state_,dt);
    e_ = sys.sources_[src_index_].OMinus(m_, sys.sources_[src_index_].GetEstMeas(state_));

    // Construct innovation covariance
    H_ = tModel::GetLinObsMatState(sys_.sources_,state_,src_index);
    V_ = tModel::GetLinObsMatSensorNoise(sys_.sources_,state_,src_index);
    F_ = tModel::GetLinTransFuncMatState(state_,dt);
    G_ = tModel::GetLinTransFuncMatNoise(state_,dt);
    HF_ = H_*F_;
    HG_ = H_*G_;
    S_inv_sqrt = (V*sys.sources_[src_index].params_.meas_cov_*V.transpose() + HG*sys.params_.process_noise_covariance_ *HG.transpose()).inverse().sqrt();
    
    // Compute Normalized Error
    e_ = S_inv_sqrt*e;

    for (int ii =0; ii < tModel::State::meas_dim_; ++ii) {
        r[ii] = e(ii,0);
    }

    }

    private:
    int src_index_;
    Eigen::Matrix<double,tModel::g_dim_*2,1> x_vector_;  /** < Places x in an Eigen vector */  
    Eigen::MatrixXd e_;         /** < Measurement error */  
    Eigen::MatrixXd S_inv_sqrt; /** < The square root  of the inverse innovation Covariance */
    Meas m_;                    /** < Measurement */
    double dt_;                 /** < Time interval between the current time and the measurement time */
    const System<tModel>& sys_; /** < A reference to the object containing all of R-RANSAC data */
    typename tModel::State state_; /** < The representation of x as a state */

    // Matrices to be used
    Eigen::MatrixXd H_;
    Eigen::MatrixXd V_;
    Eigen::MatrixXd F_;
    Eigen::MatrixXd G_;
    Eigen::MatrixXd HF_;
    Eigen::MatrixXd HG_;

};


};



} // namespace rransac






#endif // RRANSAC_TRACK_INITIALIZATION_LMLE_POLICIES_NONLINEAR_LMLE_POLICY_H_