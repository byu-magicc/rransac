#ifndef RRANSAC_TRACK_INITIALIZATION_LMLE_POLICIES_LINEAR_LMLE_POLICY_H_
#define RRANSAC_TRACK_INITIALIZATION_LMLE_POLICIES_LINEAR_LMLE_POLICY_H_
#pragma once

#include <Eigen/Dense>

#include "rransac/data_containers/cluster.h"
#include "rransac/system.h"
#include "rransac/common/measurement/measurement_base.h"


namespace rransac
{

/** \class LinearLMLEPolicy
 * This policy is used for Linear Time Invariant Models. Using a set of measurements, by which the 
 * model is observable, it produces a current state estimate. The templates for the class
 * are tModel and tSeed. tModel is the model type and tSeed is a policy to seed
 * the LMLE optimization. In the linear case, tSeed is ignored. 
 */ 

template<typename tModel, template<typename > typename tSeed>    
class LinearLMLEPolicy {

public:

typedef typename tModel::State State;           /**< The state of the target. @see State. */
typedef typename tModel::DataType DataType;     /**< The scalar object for the data. Ex. float, double, etc. */
typedef typename tModel::Source Source;         /**< The object type of the source. @see SourceBase. */
typedef tModel Model;                           /**< The object type of the model. */


/**
 * Generates a hypothetical state estimate at the current time step using the provided measurements in meas_subset. The nonlinear optimization
 * problem is solved using Ceres. Ceres uses a object type Jet for the automatic differentiation. It is because of this that the data type (DataType)
 * must be a template parameter.
 * @param[in] meas_subset The container of iterators to measurements that will be used to estimate the hypothetical state.
 * @param[in] sys The object that contains the R-RANASAC data. This includes the current time and other information.
 * @param[in,out] success A flag to indicate if the optimization converged to a solution. If and only if the optimization converged will success have a value of true.
 * @return The hypothetical state estimate of the track.
 */ 
static State GenerateHypotheticalStateEstimatePolicy(const std::vector<typename Cluster<DataType>::IteratorPair>& meas_subset, const System<tModel>& sys, bool& success);



};



/////////////////////////////////////////////////////////////////////////////////////
//                Definitions
/////////////////////////////////////////////////////////////////////////////////////
template<typename tModel, template<typename > typename tSeed>      
typename tModel::State LinearLMLEPolicy<tModel, tSeed>::GenerateHypotheticalStateEstimatePolicy(const std::vector<typename Cluster<DataType>::IteratorPair>& meas_subset, const System<tModel>& sys, bool& success){
    
    typename tModel::State x;   // hypothetical state
    success = true;

    // Matrices to be used
    Eigen::MatrixXd H;
    Eigen::MatrixXd V;
    Eigen::MatrixXd F;
    Eigen::MatrixXd G;
    Eigen::MatrixXd HF;
    Eigen::MatrixXd HG;
    Eigen::MatrixXd S_inv;
    Eigen::Matrix<double,tModel::State::g_type_::dim_*2,1> y;


    // // Assume to be identity
    // Eigen::MatrixXd V = tModel::GetLinObsMatSensorNoise(sources,x,meas_subset.begin()->inner_it->source_index);

    Eigen::Matrix<double, tModel::State::g_type_::dim_*2,tModel::State::g_type_::dim_*2> S_sum = Eigen::Matrix<double, tModel::State::g_type_::dim_*2,tModel::State::g_type_::dim_*2>::Zero();
    Eigen::Matrix<double, tModel::State::g_type_::dim_*2,1> e_sum = Eigen::Matrix<double, tModel::State::g_type_::dim_*2,1>::Zero();

    for (auto iter = meas_subset.begin(); iter != meas_subset.end(); ++iter) {


        int src_index = iter->inner_it->source_index;
        double dt = iter->inner_it->time_stamp - sys.current_time_;
        H = tModel::GetLinObsMatState(sys.sources_,x,src_index);
        V = tModel::GetLinObsMatSensorNoise(sys.sources_,x,src_index);
        F = tModel::GetLinTransFuncMatState(x,dt);
        G = tModel::GetLinTransFuncMatNoise(x,dt);
        HF = H*F;
        HG = H*G;
        

        // Builds the inverse innovation covariance 
        S_inv = (V*sys.sources_[src_index].params_.meas_cov_*V.transpose() + HG*sys.params_.process_noise_covariance_ *HG.transpose()).inverse();
        
        
        if (iter->inner_it->type == MeasurementTypes::RN_POS)
            e_sum += HF.transpose()*S_inv*iter->inner_it->pose;
        else {
            y.block(0,0,tModel::State::g_type_::dim_,1) = iter->inner_it->pose; 
            y.block(tModel::State::g_type_::dim_,0,tModel::State::g_type_::dim_,1) = iter->inner_it->twist; 
            e_sum += HF.transpose()*S_inv*y;
        }

        S_sum += HF.transpose()*S_inv*HF;

    }

    Eigen::MatrixXd tmp = S_sum.inverse()*e_sum;
    x.g_.data_ = tmp.block(0,0,tModel::State::g_type_::dim_,1);
    x.u_.data_ = tmp.block(tModel::State::g_type_::dim_,0,tModel::State::g_type_::dim_,1);

    return x;
    

}


} // namespace rransac

#endif // RRANSAC_TRACK_INITIALIZATION_LINEAR_LMLE_POLICY_H_