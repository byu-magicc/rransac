#ifndef RRANSAC_TRACK_INITIALIZATION_LMLE_POLICIES_NONLINEAR_LMLE_POLICY_H_
#define RRANSAC_TRACK_INITIALIZATION_LMLE_POLICIES_NONLINEAR_LMLE_POLICY_H_


#include "data_containers/cluster.h"
#include "system.h"
#include "common/measurement/measurement_base.h"
#include <unsupported/Eigen/MatrixFunctions>
#include <Eigen/Dense>
#include "ceres/ceres.h"
#include "state.h"
#include "lie_groups/SE3.h"


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
typedef typename State::DataType DataType;
typedef typename tModel::Source Source;
typedef tModel Model;

/**
 * Generates a hypothetical state at the current time step using the provided measurements in meas_subset.
 * @param meas_subset The container of iterators to measurements that will be used to estimate the hypothetical state
 * @param curr_time The current time
 * @param sources The vector of sources used. 
 */ 
static State GenerateHypotheticalStateEstimatePolicy(const std::vector<typename Cluster<DataType>::IteratorPair>& meas_subset, const System<tModel>& sys, bool& success);

private:

static void GenerateSeedPolicy(const std::vector<typename Cluster<DataType>::IteratorPair>& meas_subset, const System<tModel>& sys, double x[tModel::cov_dim_], const int size) {
    NonLinearLMLEPolicy::GenerateSeed(meas_subset, sys, x, size);
}

// Cost functor used in the nonlinear LMLE


struct CostFunctor {

    CostFunctor(Meas<DataType> m, const System<tModel>& sys) : m_(m), sys_(sys) {
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
    Eigen::Map<const Eigen::Matrix<T,tModel::cov_dim_,1>> x_vector(x);

    // Convert to state, propagate state to the time step of the measurement and get the error
    // between the estimated measurement and the measurement
    // typename tModel::State state;
    typedef typename tModel::template ModelTemplate<T, State::template StateTemplate> ModelT;
    typedef typename ModelT::State StateT;
    typedef typename ModelT::Source SourceT;
    typedef Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> MatXd;

    Meas<T> tmp;
    tmp.type = m_.type;
    tmp.pose = m_.pose.template cast<T>();
    tmp.twist = m_.twist.template cast<T>();
    T dt = static_cast<T>(dt_);


    StateT state;
    state.g_.data_ =  StateT::Algebra::Exp(x_vector.block(0,0, tModel::g_dim_,1));
    if (tModel::cov_dim_ != tModel::g_dim_*2) {
        state.u_.data_.setZero();
        state.u_.data_(0) = x_vector(tModel::g_dim_);
        state.u_.data_.block(StateT::Algebra::dim_t_vel_,0,StateT::Algebra::dim_a_vel_,1) =  x_vector.block(tModel::g_dim_+1,0,StateT::Algebra::dim_a_vel_,1 );

    } else {
        state.u_.data_ = x_vector.block(tModel::g_dim_,0, tModel::State::u_type_::dim_,1);
    }
    

    state = ModelT::PropagateState(state,dt);

    // std::cout << "u data: " << std::endl << state.u_.data_ << std::endl;
    // std::cout << "x vec: " << std::endl << x_vector << std::endl;

    
    Eigen::Matrix<T,Eigen::Dynamic,1> e = SourceT::OMinus(tmp, SourceT::GetEstMeas(state,m_.type));

    if (!sys_.params_.NonLinearInnovCovId_) {
        // Construct innovation covariance
        MatXd meas_cov = sys_.sources_[src_index_].params_.meas_cov_.template cast<T>();
        MatXd process_cov = sys_.params_.process_noise_covariance_.template cast<T>();
        MatXd H = SourceT::GetLinObsMatState(state,m_.type);
        MatXd V = SourceT::GetLinObsMatSensorNoise(state,m_.type);
        MatXd F = ModelT::GetLinTransFuncMatState(state,dt);
        MatXd G = ModelT::GetLinTransFuncMatNoise(state,dt);
        MatXd HF = H*F;
        MatXd HG = H*G;
        MatXd S_inv_sqrt = (V*meas_cov*V.transpose() + HG*process_cov *HG.transpose()).inverse();
        
        // Compute Normalized Error
        e = S_inv_sqrt*e;
    }


    for (unsigned int ii = 0; ii < e.rows(); ++ii) {
        r[ii] = e(ii);

    }

    return true;

    }

    private:
    int src_index_;
    Meas<DataType> m_;                    /** < Measurement */
    double dt_;                 /** < Time interval between the current time and the measurement time */
    const System<tModel>& sys_; /** < A reference to the object containing all of R-RANSAC data */



};


};






//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                            Definitions
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename tModel, template<typename > typename tSeed>     
typename tModel::State NonLinearLMLEPolicy<tModel, tSeed>::GenerateHypotheticalStateEstimatePolicy(const std::vector<typename Cluster<DataType>::IteratorPair>& meas_subset, const System<tModel>& sys, bool& success) {

    success = false;

    double x[tModel::cov_dim_];
    GenerateSeedPolicy(meas_subset, sys, x, tModel::cov_dim_);

    ceres::Problem problem;
    constexpr unsigned int meas_dim = tModel::Source::meas_dim_;

    for (auto iter = meas_subset.begin(); iter != meas_subset.end(); ++iter) {
        if (sys.sources_[iter->inner_it->source_index].params_.has_twist) {
            ceres::CostFunction* cost_function = new ceres::AutoDiffCostFunction<CostFunctor,meas_dim*2,tModel::cov_dim_>(new CostFunctor(*iter->inner_it, sys));
            problem.AddResidualBlock(cost_function, nullptr, x);
        } else {
            ceres::CostFunction* cost_function = new ceres::AutoDiffCostFunction<CostFunctor,meas_dim,tModel::cov_dim_>(new CostFunctor(*iter->inner_it, sys));
            problem.AddResidualBlock(cost_function, nullptr, x);
        }
        
        
    }


    ceres::Solver::Options options;
    // options.linear_solver_type = ceres::DENSE_QR;
    options.minimizer_progress_to_stdout = true;
    options.num_threads =sys.params_.NonLinearLMLECeresThreads_;
    options.max_num_iterations = sys.params_.NonLinearLMLECeresMaxNumIters_;
    options.logging_type = ceres::SILENT;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);

    // std::cout << summary.BriefReport() << "\n";
    // std::cout << "x: " << std::endl;
    // std::cout << "terminatin: " << summary.termination_type << std::endl;

    if (summary.termination_type == ceres::TerminationType::CONVERGENCE || summary.termination_type == ceres::TerminationType::USER_SUCCESS || summary.termination_type == ceres::TerminationType::NO_CONVERGENCE)
        success = true;

    Eigen::Matrix<double, tModel::cov_dim_, 1> x_vector;
    typename tModel::State state;

    // Convert array to Eigen vector
    x_vector = Eigen::Map<Eigen::Matrix<double,tModel::cov_dim_,1>>(x);

    // Convert to state, propagate state to the time step of the measurement and get the error
    // between the estimated measurement and the measurement
    state.g_.data_ =  state.u_.Exp(x_vector.block(0,0, tModel::State::g_type_::dim_,1));
    if (tModel::cov_dim_ != tModel::g_dim_*2) {
        state.u_.data_.setZero();
        state.u_.data_(0) = x_vector(tModel::g_dim_);
        state.u_.data_.block(tModel::State::Algebra::dim_t_vel_,0,tModel::State::Algebra::dim_a_vel_,1) =  x_vector.block(tModel::g_dim_+1,0,tModel::State::Algebra::dim_a_vel_,1 );

    } else {
        state.u_.data_ = x_vector.block(tModel::g_dim_,0, tModel::State::u_type_::dim_,1);
    }



   return state;

    

}


} // namespace rransac
#endif // RRANSAC_TRACK_INITIALIZATION_LMLE_POLICIES_NONLINEAR_LMLE_POLICY_H_