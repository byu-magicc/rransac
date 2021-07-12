#ifndef RRANSAC_TRACK_INITIALIZATION_LMLE_POLICIES_NONLINEAR_LMLE_POLICY_H_
#define RRANSAC_TRACK_INITIALIZATION_LMLE_POLICIES_NONLINEAR_LMLE_POLICY_H_
#pragma once

#include <unsupported/Eigen/MatrixFunctions>
#include <Eigen/Dense>
#include "ceres/ceres.h"
#include <typeinfo>

#include "rransac/data_containers/cluster.h"
#include "rransac/system.h"
#include "rransac/common/measurement/measurement_base.h"
#include "lie_groups/state.h"
#include "lie_groups/lie_groups/SE3.h"


namespace rransac
{


/** \class NonLinearLMLEPolicy
 * This policy uses Ceres to perform nonlinear optimization on a log maximum likelihood estimation problem.
 * The optimization problem seeks to minimize the log maximum likelihood of a subset of measurements conditioned
 * on the current state estimate by optimizing over the current state estimate. Since the optimization problem is nonlinear,
 * only a local solution is guaranteed provided that it exists. The templates for the class
 * are tModel and tSeed. tModel is the model type and tSeed is a policy to seed
 * the LMLE optimization. 
 */ 

template<typename _Model, template<typename > typename _Seed>    
class NonLinearLMLEPolicy : public _Seed<_Model> {

public:

typedef typename _Model::State State;           /**< The state of the target. @see State. */
typedef typename _Model::DataType DataType;     /**< The scalar object for the data. Ex. float, double, etc. */
typedef _Model Model;                           /**< The object type of the model. */
typedef typename Model::Base::Measurement Measurement;
typedef typename Model::Base::TransformDataType TransformDataType;
typedef Cluster<DataType,TransformDataType> ClusterT;
typedef System<Model> Sys;
static constexpr unsigned int cov_dim_ = Model::Base::cov_dim_;
// typedef typename tModel::template ModelTemplate<ceres::Jet<double,tModel::cov_dim_>> ModelT;
// typedef typename ModelT::State StateT;
// typedef typename ModelT::SourceContainer SourceContainerT;





/**
 * Generates a hypothetical state estimate at the current time step using the provided measurements in meas_subset. The nonlinear optimization
 * problem is solved using Ceres. The parameters of the nonlinear optimization problem is the current state estimate in local coordinates. That is, let \f$ x \f$ denote
 * the parameters of the nonlinear optimization problem, then the current state estimate is \f$ \text{Log}\left(x\right)\f$. The optimization problem seeks to minimize the 
 * log maximum likelihood of the measurements in the provided subset conditioned on the current state estimate by optimizing over the current state estimate in local coordinates.
 * @param[in] meas_subset The container of iterators to measurements that will be used to estimate the hypothetical state.
 * @param[in] sys The object that contains the R-RANASAC data. This includes the current time and other information.
 * @param[in,out] success A flag to indicate if the optimization converged to a solution. If and only if the optimization converged will success have a value of true.
 * @return The hypothetical state estimate of the track.
 */ 
static State GenerateHypotheticalStateEstimatePolicy(const std::vector<typename ClusterT::IteratorPair>& meas_subset, const Sys& sys, bool& success);

private:

/**
 * Nonlinear optimization algorithms are only guaranteed to converge to local minimums or a local solution. They can also take a while to converge. The nonlinear optimization problem
 * can be seeded with certain initial conditions in order to speed up convergence and to converge on a desireable local minimum. This function calls the GenerateSeedPolicy member function
 * from the policy class specified by tSeed in order to seed the initial conditions of the optimization problem.
 * @param[in] meas_subset A subset of measurements used to calculate the seed.
 * @param[in] sys The object that contains the R-RANASAC data.
 * @param[in] x The initial conditions of the optimization problem. The initial conditions will change according to the seed policy. 
 * @param[in] size The number parameters the optimization solver is optimizing over which is the size of the input x.
 */ 
static void GenerateSeed(const std::vector<typename ClusterT::IteratorPair>& meas_subset, const Sys& sys, double x[cov_dim_], const int size) {
    NonLinearLMLEPolicy::GenerateSeedPolicy(meas_subset, sys, x, size);
}




/**
 * Ceres builds the optimization problem using residual blocks composed of cost functors. These residual blocks calculate
 * a portion of the error in the total optimization problem. In this specific case. This functor is given one of the measurements
 * from the measurement subset used to generate the estimated state. It takes the current state estimate and propagets it in time 
 * to the same time as the measurement and computes the error between the measurement and state normalized by the innovation covariance. 
 * The error is only normalized by the innovation covariance if the parameter Parameters::nonlinear_innov_cov_id_ is set to false.
 */ 
struct CostFunctor {

    /**
     * Sets the measurement used to compute the error, the time interval between the current time and the
     * time stamp of the measurement, and the source index of the measurement. 
     * @param m The measurement to be used to compute the error.
     * @param sys The object that contains the R-RANSAC data. This includes the current time. 
     */ 
    CostFunctor(Measurement m, const Sys& sys) : m_(m), sys_(sys) {
        dt_ = m.time_stamp - sys_.current_time_;
        src_index_ = m_.source_index;
    }

    /**
     * Computes the error between the measurement and the current state estimate normalized by the innovation covariance. 
     * Since the current state estimate is represented in local coordinates, it maps it back onto the manfold, propagates
     * the state estimate to the same time as the time stamp on the measurement, computes the estimated measurement from 
     * the propagated state estimate, and then uses the estimated measurement and the measurement to compute the error 
     * normalized by the innovation covariance. This error is the residual of the operator. 
     * The error is only normalized by the innovation covariance if the parameter Parameters::nonlinear_innov_cov_id_ is set to false.
     * 
     * This function is templated because Ceres uses a data type called Jet in order to calculate the derivative using
     * automatic differentiation. Because of this, many of the other classes in R-RANSAC have the template parameter DataType
     * so that they are compatible with Ceres. 
     */ 
    template <typename T>
    bool operator() (const T* const x,  T* r)  const {


        // std::cout << "type id: " << typeid(T).name() << std::endl;
        
        // Since Ceres uses the data type Jet, we must create the model, state, and source
        // using the data type Jet in order for the automatic differentiation to work properly.
        typedef typename Model::template ModelTemplate<T> ModelT;
        typedef typename ModelT::State StateT;
        typedef typename ModelT::SourceContainer SourceContainerT;
        static SourceContainerT source_container_t;
        // typedef typename ModelT::Source SourceT;
        typedef Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> MatXd;

        // Convert array to Eigen vector
        Eigen::Map<const Eigen::Matrix<T,cov_dim_,1>> x_vector(x);

        // std::cout << "x vec; " << std::endl << x_vector << std::endl;

        // Convert the measurement and time interval to type Jet.
        Meas<T,typename ModelT::Base::TransformDataType> tmp;
        tmp.type = m_.type;
        tmp.pose = m_.pose.template cast<T>();
        tmp.twist = m_.twist.template cast<T>();
        tmp.transform_data_t_m = m_.transform_data_t_m.template cast<T>();
        T dt = static_cast<T>(dt_);


        // The parameter x is the local coordinate representation of the state. We must map it back to
        // the manifold. 
        StateT state = StateT::Identity(); // State is at identity.    
        ModelT::OPlus(state,x_vector);       

        // std::cout << "state twist: " << std::endl << state.u_.data_ << std::endl;

        // Propagate the state from the current time stamp to the time stamp of the measurement. 
        ModelT::PropagateState(state,dt);

        // Compute the error   
        // std::cout << "src index: " << src_index_ << std::endl; 
        // std::cout << "tmp pose: " << std::endl << tmp.pose << std::endl; 
        // std::cout << "tmp twist: " << std::endl << tmp.twist << std::endl; 
        Eigen::Matrix<T,Eigen::Dynamic,1> e = SourceContainerT::OMinus(src_index_,tmp, SourceContainerT::GetEstMeas(src_index_, state,m_.transform_state,tmp.transform_data_t_m));

        // Compute the inverse innovation covariance and normalize the error. 
        if (!sys_.params_.nonlinear_innov_cov_id_) {
            // Construct innovation covariance
            MatXd meas_cov = sys_.source_container_.GetParams(src_index_).meas_cov_.template cast<T>();
            MatXd H = source_container_t.GetLinObsMatState(src_index_,state,m_.transform_state,tmp.transform_data_t_m);
            MatXd V = source_container_t.GetLinObsMatSensorNoise(src_index_,state,m_.transform_state,tmp.transform_data_t_m);

            MatXd F = ModelT::GetLinTransFuncMatState(state,dt);
            MatXd HF = H*F;
            MatXd S_inv_sqrt = (V*meas_cov*V.transpose()).inverse();

            
            // Compute Normalized Error
            e = S_inv_sqrt*e;
        }

        // std::cout << "error: " << std::endl << e << std::endl;


        // Extract the residual from the error. 
        for (unsigned int ii = 0; ii < e.rows(); ++ii) {
            r[ii] = e(ii);

        }




    return true;

    }





    private:
    int src_index_;             /**< The source index of the measurement. */
    Measurement m_;          /**< Measurement */
    double dt_;                 /**< Time interval between the current time and the measurement time. */
    const Sys& sys_; /**< A reference to the object containing all of R-RANSAC data. */



};


};






//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                            Definitions
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename _Model, template<typename > typename _Seed>     
typename NonLinearLMLEPolicy<_Model, _Seed>::State NonLinearLMLEPolicy<_Model, _Seed>::GenerateHypotheticalStateEstimatePolicy(const std::vector<typename ClusterT::IteratorPair>& meas_subset, const Sys& sys, bool& success) {

    success = false;

    double x[cov_dim_];

    // Construct the seed for the nonlinear optimization problem
    GenerateSeed(meas_subset, sys, x, cov_dim_);

    // Use Ceres to build the optimization problem
    ceres::Problem problem;
    for (auto iter = meas_subset.begin(); iter != meas_subset.end(); ++iter) {
        const int total_meas_dim_ = sys.source_container_.GetParams(iter->inner_it->source_index).total_meas_dim_;
        ceres::CostFunction* cost_function = new ceres::AutoDiffCostFunction<CostFunctor,ceres::DYNAMIC,cov_dim_>(new CostFunctor(*iter->inner_it, sys),total_meas_dim_);
            problem.AddResidualBlock(cost_function, nullptr, x);      
        
    }


    ceres::Solver::Options options;
    options.minimizer_progress_to_stdout = true;
    options.num_threads =sys.params_.nonlinear_LMLE_Ceres_threads_;
    options.max_num_iterations = sys.params_.nonlinear_LMLE_Ceres_max_num_iters_;
    options.logging_type = ceres::SILENT;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);

    // std::cout << summary.BriefReport() << "\n";
    // std::cout << "x: " << std::endl;
    // std::cout << "terminatin: " << summary.termination_type << std::endl;

    if (summary.termination_type == ceres::TerminationType::CONVERGENCE || summary.termination_type == ceres::TerminationType::USER_SUCCESS || summary.termination_type == ceres::TerminationType::NO_CONVERGENCE)
        success = true;

    Eigen::Matrix<double, cov_dim_, 1> x_vector;
    State state;

    // Convert array to Eigen vector
    x_vector = Eigen::Map<Eigen::Matrix<double,cov_dim_,1>>(x);

    Model::OPlus(state,x_vector);



   return state;

    

}


} // namespace rransac
#endif // RRANSAC_TRACK_INITIALIZATION_LMLE_POLICIES_NONLINEAR_LMLE_POLICY_H_