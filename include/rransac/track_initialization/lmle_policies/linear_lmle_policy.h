#ifndef RRANSAC_TRACK_INITIALIZATION_LMLE_POLICIES_LINEAR_LMLE_POLICY_H_
#define RRANSAC_TRACK_INITIALIZATION_LMLE_POLICIES_LINEAR_LMLE_POLICY_H_
#pragma once

#include <Eigen/Dense>

#include "rransac/data_containers/cluster.h"
#include "rransac/system.h"
#include "rransac/common/measurement/measurement_base.h"
#include "rransac/common/utilities.h"

namespace rransac
{

/** \class LinearLMLEPolicy
 * This policy is used for Linear Time Invariant Models. Using a set of measurements, by which the 
 * model is observable, it produces a current state estimate. The templates for the class
 * are Model and tSeed. Model is the model type and tSeed is a policy to seed
 * the LMLE optimization. In the linear case, tSeed is ignored. 
 */ 

template<typename _Model, template<typename > typename _Seed>    
class LinearLMLEPolicy {

public:

typedef typename _Model::State State;                                               /**< The state of the target. @see State. */
typedef typename _Model::DataType DataType;                                         /**< The scalar object for the data. Ex. float, double, etc. */
typedef _Model Model;                                                               /**< The object type of the model. */
typedef utilities::MatXT<typename State::DataType> MatXT;                           /**< Dynamic Matrix. */
typedef typename Model::MatModelCov MatModelCov;                                         /**< The error covariance data type. */
typedef typename Model::Base::VecCov VecCov;                                        /**< The error state data type. */
typedef typename Model::Base::TransformDataType TransformDataType;
typedef typename Model::Base::Measurement Measurement;
typedef Cluster<DataType,TransformDataType> ClusterT;

/**
 * Generates a hypothetical state estimate at the current time step using the provided measurements in meas_subset. The nonlinear optimization
 * problem is solved using Ceres. Ceres uses a object type Jet for the automatic differentiation. It is because of this that the data type (DataType)
 * must be a template parameter.
 * @param[in] meas_subset The container of iterators to measurements that will be used to estimate the hypothetical state.
 * @param[in] sys The object that contains the R-RANASAC data. This includes the current time and other information.
 * @param[in,out] success A flag to indicate if the optimization converged to a solution. If and only if the optimization converged will success have a value of true.
 * @return The hypothetical state estimate of the track.
 */ 
static State GenerateHypotheticalStateEstimatePolicy(const std::vector<typename ClusterT::IteratorPair>& meas_subset, const System<Model>& sys, bool& success);


};



/////////////////////////////////////////////////////////////////////////////////////
//                Definitions
/////////////////////////////////////////////////////////////////////////////////////
template<typename _Model, template<typename > typename _Seed>      
typename _Model::State LinearLMLEPolicy<_Model, _Seed>::GenerateHypotheticalStateEstimatePolicy(const std::vector<typename ClusterT::IteratorPair>& meas_subset, const System<Model>& sys, bool& success){
    
    typename _Model::State x;   // hypothetical state
    success = true;
    MatXT H;
    MatXT V;
    MatXT F;
    MatXT G;
    MatXT HF;
    MatXT HG;
    MatXT S_inv;


    // // Assume to be identity
    // Eigen::MatrixXd V = Model::GetLinObsMatSensorNoise(sources,x,meas_subset.begin()->inner_it->source_index);

    MatModelCov S_sum(MatModelCov::Zero());
    VecCov e_sum(VecCov::Zero());
    

    for (auto iter = meas_subset.begin(); iter != meas_subset.end(); ++iter) {


        int src_index = iter->inner_it->source_index;
        double dt = iter->inner_it->time_stamp - sys.current_time_;
        H = sys.source_container_.GetLinObsMatState(src_index,x,iter->inner_it->transform_state, iter->inner_it->transform_data_t_m);
        V = sys.source_container_.GetLinObsMatSensorNoise(src_index,x,iter->inner_it->transform_state, iter->inner_it->transform_data_t_m);
        F = Model::GetLinTransFuncMatState(x,dt);
        HF = H*F;
        const unsigned int& total_meas_dim = sys.source_container_.GetParams(src_index).total_meas_dim_;
        const unsigned int& meas_pose_rows = sys.source_container_.GetParams(src_index).meas_pose_rows_;
        Eigen::Matrix<typename Model::DataType, Eigen::Dynamic,Eigen::Dynamic> vec_meas(total_meas_dim,1);
        vec_meas.block(0,0,meas_pose_rows,1 ) = iter->inner_it->pose;
        if (sys.source_container_.GetParams(src_index).has_vel_) {
            const unsigned int& meas_twist_rows = sys.source_container_.GetParams(src_index).meas_twist_rows_;
            vec_meas.block(meas_pose_rows,0,meas_twist_rows,1) = iter->inner_it->twist;
        }

        // Builds the inverse innovation covariance 
        S_inv = (V*sys.source_container_.GetParams(src_index).meas_cov_*V.transpose()).inverse();
        e_sum += HF.transpose()*S_inv*vec_meas;     
        S_sum += HF.transpose()*S_inv*HF;

    }

    VecCov tmp = S_sum.inverse()*e_sum;
    x.g_.data_ = tmp.block(0,0,Model::State::Group::dim_,1);
    x.u_.data_ = tmp.block(Model::State::Group::dim_,0,Model::State::Algebra::total_num_dim_,1);

    return x;
    

}


} // namespace rransac

#endif // RRANSAC_TRACK_INITIALIZATION_LINEAR_LMLE_POLICY_H_