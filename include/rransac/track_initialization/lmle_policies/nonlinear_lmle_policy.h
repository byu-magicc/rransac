#ifndef RRANSAC_TRACK_INITIALIZATION_LMLE_POLICIES_NONLINEAR_LMLE_POLICY_H_
#define RRANSAC_TRACK_INITIALIZATION_LMLE_POLICIES_NONLINEAR_LMLE_POLICY_H_


#include "data_containers/cluster.h"
#include "system.h"
#include "common/measurement/measurement_base.h"
#include <unsupported/Eigen/MatrixFunctions>
#include <Eigen/Dense>
#include "ceres/ceres.h"


namespace rransac
{


/** \class NonLinearLMLEPolicy
 * This policy is used for Nonlinear Time Invariant Models. Using a set of measurements, by which the 
 * model is at least locally observable, it produces a current state estimate. The templates for the class
 * are tModel and tSeed. tModel is the model type and tSeed is a policy to seed
 * the LMLE optimization. 
 */ 

template<typename tModel, template<typename > typename tSeed>    
class NonLinearLMLEPolicy : public tSeed<tModel> {

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
static State GenerateHypotheticalStateEstimatePolicy(const std::vector<Cluster::IteratorPair>& meas_subset, const System<tModel>& sys);

private:

static void GenerateSeedPolicy(const std::vector<Cluster::IteratorPair>& meas_subset, const System<tModel>& sys, double x[tModel::State::g_type_::dim_*2], const int size) {
    NonLinearLMLEPolicy::GenerateSeed(meas_subset, sys, x, size);
}

// Cost functor used in the nonlinear LMLE


struct CostFunctor {

    CostFunctor(Meas m, const System<tModel>& sys) : m_(m), sys_(sys) {
        dt_ = m.time_stamp - sys_.current_time_;
        src_index_ = m_.source_index;
    }

    /**
     * The components of x is the Lie algebra of the state
     * 
     */ 
    template <typename T>
    bool operator() (const T* const x,  T* r)  const {

    // Convert array to Eigen vector
    const Eigen::Map<const Eigen::Matrix<T,tModel::g_dim_*2,1>> x_vector(x);
    // for (int ii = 0; ii < tModel::g_dim_*2; ++ii)
    //     x_vector_(ii,0) = x[ii];

    // Convert to state, propagate state to the time step of the measurement and get the error
    // between the estimated measurement and the measurement
    typename tModel::State state;
    state.g_.data_ =  state.u_.Exp(x_vector.block(0,0, tModel::g_dim_,1));
    state.u_.data_ = x_vector.block(tModel::g_dim_,0, tModel::State::u_type_::dim_,1);
    state = tModel::PropagateState(state,dt_);
    Eigen::Matrix<T,Eigen::Dynamic,1> e = sys_.sources_[src_index_].OMinus(m_, sys_.sources_[src_index_].GetEstMeas(state));

    // Construct innovation covariance
    Eigen::MatrixXd H = tModel::GetLinObsMatState(sys_.sources_,state,src_index_);
    Eigen::MatrixXd V = tModel::GetLinObsMatSensorNoise(sys_.sources_,state,src_index_);
    Eigen::MatrixXd F = tModel::GetLinTransFuncMatState(state,dt_);
    Eigen::MatrixXd G = tModel::GetLinTransFuncMatNoise(state,dt_);
    Eigen::MatrixXd HF = H*F;
    Eigen::MatrixXd HG = H*G;
    Eigen::MatrixXd S_inv_sqrt = (V*sys_.sources_[src_index_].params_.meas_cov_*V.transpose() + HG*sys_.params_.process_noise_covariance_ *HG.transpose()).inverse().sqrt();
    
    // Compute Normalized Error
    e = S_inv_sqrt*e;
    
    r = e.data();

    return true;

    }

    private:
    int src_index_;
    Meas m_;                    /** < Measurement */
    double dt_;                 /** < Time interval between the current time and the measurement time */
    const System<tModel>& sys_; /** < A reference to the object containing all of R-RANSAC data */



};


};






//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                            Definitions
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename tModel, template<typename > typename tSeed>     
typename tModel::State NonLinearLMLEPolicy<tModel, tSeed>::GenerateHypotheticalStateEstimatePolicy(const std::vector<Cluster::IteratorPair>& meas_subset, const System<tModel>& sys) {

    double x[tModel::State::g_type_::dim_*2];
    GenerateSeedPolicy(meas_subset, sys, x, tModel::State::g_type_::dim_*2);

    ceres::Problem problem;

    for (auto iter = meas_subset.begin(); iter != meas_subset.end(); ++iter) {
        ceres::CostFunction* cost_function = new ceres::AutoDiffCostFunction<CostFunctor,tModel::State::g_type_::dim_*2,tModel::State::g_type_::dim_*2>(new CostFunctor(*iter->inner_it, sys));
        problem.AddResidualBlock(cost_function, nullptr, x);
    }



    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR;
    options.minimizer_progress_to_stdout = true;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);

    std::cout << summary.BriefReport() << "\n";
    std::cout << "x: " << std::endl;

    Eigen::Matrix<double, tModel::g_dim_*2, 1> x_vector;
    typename tModel::State state;

    //     // Convert array to Eigen vector
    // for (int ii = 0; ii < tModel::g_dim_*2; ++ii)
    //     x_vector(ii,0) = x[ii];
    x_vector = Eigen::Map<Eigen::Matrix<double,tModel::g_dim_*2,1>>(x);

    // Convert to state, propagate state to the time step of the measurement and get the error
    // between the estimated measurement and the measurement
    state.g_.data_ =  state.u_.Exp(x_vector.block(0,0, tModel::State::g_type_::dim_,1));
    state.u_.data_ = x_vector.block(tModel::State::g_type_::dim_,0, tModel::State::u_type_::dim_,1);

    std::cout << "g: " << std::endl << state.g_.data_ << std::endl;
    std::cout << "u:" << std::endl << state.u_.data_ << std::endl;

   return state;

    

}


} // namespace rransac
#endif // RRANSAC_TRACK_INITIALIZATION_LMLE_POLICIES_NONLINEAR_LMLE_POLICY_H_