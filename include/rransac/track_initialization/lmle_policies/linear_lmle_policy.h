#ifndef RRANSAC_TRACK_INITIALIZATION_LMLE_POLICIES_LINEAR_LMLE_POLICY_H_
#define RRANSAC_TRACK_INITIALIZATION_LMLE_POLICIES_LINEAR_LMLE_POLICY_H_

#include "data_containers/cluster.h"
#include "system.h"
#include "common/measurement/measurement_base.h"
#include <Eigen/Dense>

namespace rransac
{

template<typename tModel>    
class LinearLMLEPolicy {

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
static State GenerateStateEstimatePolicy(const std::vector<Cluster::ConstIteratorPair>& meas_subset, const System<tModel>& sys);



};



/////////////////////////////////////////////////////////////////////////////////////
//                Definitions
/////////////////////////////////////////////////////////////////////////////////////
template<typename tModel>    
typename tModel::State LinearLMLEPolicy<tModel>::GenerateStateEstimatePolicy(const std::vector<Cluster::ConstIteratorPair>& meas_subset, const System<tModel>& sys){
    
    typename tModel::State x;   // hypothetical state
    

    // Matrices to be used
    Eigen::MatrixXd H;
    Eigen::MatrixXd V;
    Eigen::MatrixXd F;
    Eigen::MatrixXd G;
    Eigen::MatrixXd HF;
    Eigen::MatrixXd HG;
    Eigen::MatrixXd S_inv;
    Eigen::Matrix<double,tModel::State::g_type_::dim_*2,1> y


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
        

        // Build innovation covariance depending on if the measurement covariance is fixed or not. 
        if (sys.sources_[src_index].params_.meas_cov_fixed_) {
            S_inv = (V*sys.sources_[src_index].params_.meas_cov_*V.transpose() + HG*sys.params_.process_noise_covariance_ *HG.transpose()).inverse();
        } else {
            S_inv = (V*iter->inner_it->meas_cov*V.transpose()+ HG*sys.params_.process_noise_covariance_ *HG.transpose()).inverse();
        }
        
        
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