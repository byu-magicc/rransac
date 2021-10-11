#ifndef RRANSAC_COMMON_MODELS_CENTRALIZED_MEASUREMENT_FUSION_H_
#define RRANSAC_COMMON_MODELS_CENTRALIZED_MEASUREMENT_FUSION_H_

#include <vector>
#include "rransac/common/sources/source_container.h"
#include "rransac/common/measurement/measurement_base.h"


namespace rransac
{
    

template< typename _ModelBase>
class CentralizedMeasurementFusion {
public:

typedef _ModelBase ModelBase;
typedef typename _ModelBase::DerivedModel Model;                              /**< The model type **/
typedef typename ModelBase::State State;                                      /**< The state of the target. @see State. */
typedef typename ModelBase::DataType DataType;                                /**< The scalar object for the data. Ex. float, double, etc. */
typedef typename ModelBase::SourceContainer SourceContainer;                  /**< The object type of the source. @see SourceBase. */
typedef typename ModelBase::MatModelCov MatModelCov;                          /**< The object type of the error covariance, Jacobians, and others. */
typedef typename ModelBase::VecCov VecCov;                                    /**< The object type of state update. */
typedef typename ModelBase::MatS MatS;                                        /**< The object type of the innovation covariance. */
typedef typename SourceContainer::Source0::MatH MatH;                         /**< The object type of the Jacobian of the observation function w.r.t. to the state. */
typedef typename SourceContainer::Source0::MatV MatV;                         /**< The object type of the Jacobian of the observation function w.r.t. the measurement noise. */
typedef typename SourceContainer::Measurement Measurement;


static constexpr unsigned int cov_dim_ = Model::cov_dim_;                            /**< The dimension of the error covariance. */
static constexpr unsigned int num_sources_ = SourceContainer::num_sources_;  /**< The number of sources in the source container. */


void PerformCentralizedMeasurementFusion(const SourceContainer& source_container,  const bool  sequential_else_parallel_fusion,  Model& model) {
    if (sequential_else_parallel_fusion)
    {
        PerformSequentialCentralizedMeasurementFusion(source_container,model);
    } else {
        PerformParallelCentralizedMeasurementFusion(source_container,model);
    }
    
}


const std::vector<unsigned int>& GetSourcesProducedMeasurements() const {return sources_produced_measurements_;};

unsigned int GetNumSourcesProducedMeasurements() const {return num_sources_produced_measurements_;}

const std::vector<unsigned int>& GetNumSingleAssociationEvents() const {return num_single_association_events_;}

unsigned int GetNumJointAssociationEvents() const {return num_joint_association_events_;}

const std::vector<std::vector<double>>& GetSingleAssociationEventsWeights() const {return single_association_events_weights_;}

const std::vector<MatModelCov>& GetJointAssociationCov()const {return joint_association_cov_;} 


private:

    std::vector<unsigned int> sources_produced_measurements_;
    unsigned int num_sources_produced_measurements_;
    std::vector<unsigned int> num_single_association_events_;
    unsigned int num_joint_association_events_;
    std::vector<std::vector<double>> single_association_events_weights_;


    std::vector<unsigned int> joint_association_index_splice_;

    std::vector<MatModelCov> joint_association_cov_;
    std::vector<std::vector<VecCov>> HT_VRVTinv_nu_;



    void PerformParallelCentralizedMeasurementFusion(const SourceContainer& source_container,  Model& model);


    //
    void PerformSequentialCentralizedMeasurementFusion(const SourceContainer& source_container,  Model& model) {
        for(auto& meas: model.new_assoc_meas_){
            if(meas.size() > 0)
            {
                PerformUpdateFromSingleSensor(source_container,meas,model);
            }
        }
    }



    void ConstructCovAndHT_VRVTinv_nu(const SourceContainer& source_container, const Model& model);


    std::vector<unsigned int> JointToSingleAssociationIndex(const unsigned int joint_association_index);

    std::vector<std::vector<double>> GetSingleAssociationEventWeights(const Model& model);

    unsigned int SingleAssociationIndexToAssociationCovIndex(const std::vector<unsigned int>& single_association_index);

    std::vector<bool> SourcesAssociatedWithAssociationCov(const unsigned int cov_index);


    /**
     * Updates the state estimate and error covariance using the associated measurements from a single sensor.
     * @param[in] source_container The container of all of the sources
     * @param[in] meas          All of the new measurements produced by a single source.
     */ 
    void PerformUpdateFromSingleSensor(const SourceContainer& source_container, const std::vector<Measurement>& meas, Model& model);

};

//------------------------------------------------------------------------------------------------------------------------------
//                                  Definitions
//------------------------------------------------------------------------------------------------------------------------------
template<typename _Model>
void CentralizedMeasurementFusion<_Model>::PerformParallelCentralizedMeasurementFusion(const SourceContainer& source_container,  Model& model) {

        sources_produced_measurements_.clear();
        num_single_association_events_.clear();
        single_association_events_weights_ = GetSingleAssociationEventWeights(model);
        joint_association_index_splice_.clear();
        joint_association_cov_.clear();
        num_joint_association_events_ = 1;
        HT_VRVTinv_nu_.clear();



        // Get the sources that produced measurements, the number of measurements and association events for each source that produced measurements
        for(auto& meas: model.new_assoc_meas_){
            if (meas.size() > 0) {
                sources_produced_measurements_.push_back(meas.front().source_index);
                num_single_association_events_.push_back(meas.size()+1);
                num_joint_association_events_ *= (meas.size()+1);
            } else {
                num_single_association_events_.push_back(0);
            }
        }
        num_sources_produced_measurements_ = sources_produced_measurements_.size();

        for (unsigned int ii = 0; ii < num_sources_produced_measurements_; ++ii) {
            if (ii ==0) {
                joint_association_index_splice_.push_back(1);
            }else {
                joint_association_index_splice_.push_back(joint_association_index_splice_[ii-1]*num_single_association_events_[sources_produced_measurements_[ii-1]]);
            }
        }

        ConstructCovAndHT_VRVTinv_nu(source_container, model);


        VecCov mu_sum, mu_a, sum_err;
        mu_sum.setZero();
        double weight = 0;
        std::vector<unsigned int> single_association_indicies;
        unsigned int association_cov_index;
        size_t source_index = 0;
        size_t meas_index = 0;
        // MatModelCov& P;
        MatModelCov P_sum, mu_cov_sum;
        P_sum.setZero();
        mu_cov_sum.setZero();
        // std::cout << "sources produced measurements: " << std::endl;
        // for(int ii = 0; ii < sources_produced_measurements_.size(); ++ii) {
        //     std::cout << sources_produced_measurements_[ii] << " ";
        // }
        // std::cout << std::endl;
        for(size_t joint_association_index = 0; joint_association_index< num_joint_association_events_; ++ joint_association_index) {

            // std::cerr << "joint_association_index: " << joint_association_index << std::endl;

            weight = 1;
            sum_err.setZero();
            // std::cout << "here1" << std::endl;
            single_association_indicies = JointToSingleAssociationIndex(joint_association_index);
            // std::cout << "single association indicies " << std::endl;
            // for (int ii = 0; ii < single_association_indicies.size(); ++ii) {
            //     std::cout << single_association_indicies[ii] << ", ";
            // }
            // std::cout << std::endl;
            // std::cout << "here2" << std::endl;

            association_cov_index = SingleAssociationIndexToAssociationCovIndex(single_association_indicies);
            // std::cout << "cov index" << association_cov_index << std::endl;
            // std::cout << "here3" << std::endl;

            MatModelCov&P = joint_association_cov_[association_cov_index];

            for(size_t sai = 0; sai< num_sources_produced_measurements_; ++sai) {
                source_index = sources_produced_measurements_[sai];
                weight *= single_association_events_weights_[source_index][single_association_indicies[sai]];
                if(single_association_indicies[sai] > 0) {
                    meas_index = single_association_indicies[sai] - 1;
                    sum_err += HT_VRVTinv_nu_[source_index][meas_index];
                }
            }

            // std::cout << "here4" << std::endl;


            mu_a = P*sum_err;
            mu_sum += mu_a*weight;
            P_sum += P*weight;
            mu_cov_sum += weight*mu_a*mu_a.transpose();



        }

            // std::cout << "here5" << std::endl;


        model.err_cov_ = P_sum + mu_cov_sum - mu_sum*mu_sum.transpose();     

        // Reset mean to zero   
        model.OPlusEQ(mu_sum);

        MatModelCov Jr = Model::Jr(mu_sum);
        model.err_cov_ = Jr * model.err_cov_ *Jr.transpose();
        

}

//------------------------------------------------------------------------------------------------------------------------------
template<typename _Model>
void CentralizedMeasurementFusion<_Model>::ConstructCovAndHT_VRVTinv_nu(const SourceContainer& source_container, const Model& model) {

std::vector<MatModelCov> HT_VRVTinv_H(num_sources_);
Eigen::Matrix<DataType,Eigen::Dynamic,Eigen::Dynamic> H, V, VRVTinv, HT_VRVTinv;
Measurement estimated_meas;
MatModelCov err_cov_inv = model.err_cov_.inverse();

HT_VRVTinv_nu_.clear();
joint_association_cov_.clear();
HT_VRVTinv_nu_.resize(num_sources_);

// Construct HTRinvH and the measurement weights HT_VRVTinv_nu where nu is the innovation term
for(unsigned int source_index: sources_produced_measurements_) {
    H = source_container.GetLinObsMatState(source_index, model.state_, model.new_assoc_meas_[source_index].front().transform_state, model.new_assoc_meas_[source_index].front().transform_data_t_m);
    V = source_container.GetLinObsMatSensorNoise(source_index, model.state_, model.new_assoc_meas_[source_index].front().transform_state, model.new_assoc_meas_[source_index].front().transform_data_t_m);
    VRVTinv = (V*source_container.GetParams(source_index).meas_cov_*V.transpose()).inverse();
    HT_VRVTinv = H.transpose()*VRVTinv;
    HT_VRVTinv_H[source_index] = HT_VRVTinv*H;

    estimated_meas = source_container.GetEstMeas(source_index, model.state_, model.new_assoc_meas_[source_index].front().transform_state, model.new_assoc_meas_[source_index].front().transform_data_t_m); 
    for (auto& meas : model.new_assoc_meas_[source_index]) {

        HT_VRVTinv_nu_[source_index].push_back( HT_VRVTinv*source_container.OMinus(source_index,meas,estimated_meas) );
    }

}

std::vector<bool> source_associated_with_association_cov;
MatModelCov association_cov;
double num_cov = powf(2,num_sources_produced_measurements_);
for (unsigned int ii = 0; ii < num_cov; ++ii) {
    association_cov = err_cov_inv;
    source_associated_with_association_cov = SourcesAssociatedWithAssociationCov(ii);

    for (unsigned int jj = 0; jj < num_sources_produced_measurements_; ++jj) {
        if ( source_associated_with_association_cov[jj]) {
            association_cov += HT_VRVTinv_H[sources_produced_measurements_[jj]];
        }
    }

    joint_association_cov_.push_back( association_cov.inverse());

}

}

//------------------------------------------------------------------------------------------------------------------------------
template<typename _Model>
std::vector<unsigned int> CentralizedMeasurementFusion<_Model>::JointToSingleAssociationIndex(const unsigned int joint_association_index) {

    std::vector<unsigned int> single_association_index(num_sources_produced_measurements_);

    for (unsigned int ii = 0; ii < num_sources_produced_measurements_; ++ii) {
        // std::cout << "ii: " << ii << std::endl;
        // std::cout << " joint_association_index_splice_[ii]  " <<  joint_association_index_splice_[ii] << std::endl;
        single_association_index[ii] = ( joint_association_index / joint_association_index_splice_[ii] ) % num_single_association_events_[sources_produced_measurements_[ii]];
    }

    return single_association_index;
}

//------------------------------------------------------------------------------------------------------------------------------
template<typename _Model>
std::vector<std::vector<double>> CentralizedMeasurementFusion<_Model>::GetSingleAssociationEventWeights(const Model& model) {

    std::vector<std::vector<double>> single_association_events_weights;
    std::vector<double> weights;

    double B0 = 1;

    for(auto& measurements : model.new_assoc_meas_) {
        
        weights.clear();
        if(measurements.size() > 0) {
            B0 = 1;
            weights.push_back(B0);

            for (auto& meas : measurements) {
                weights.push_back(meas.weight);
                B0 -= meas.weight;
            }
            weights[0] = B0;
        }

        single_association_events_weights.push_back(weights);
    }

    return single_association_events_weights;

}

//------------------------------------------------------------------------------------------------------------------------------

template<typename _Model>
unsigned int CentralizedMeasurementFusion<_Model>::SingleAssociationIndexToAssociationCovIndex(const std::vector<unsigned int>& single_association_index) {
    unsigned int index = 0;

    for (int ii = 0; ii < single_association_index.size(); ++ ii) {

        if (single_association_index[ii] > 0) {
            index += powf(2,ii);
        }

    }

    return index;
}

//------------------------------------------------------------------------------------------------------------------------------

template<typename _Model>
std::vector<bool> CentralizedMeasurementFusion<_Model>::SourcesAssociatedWithAssociationCov(const unsigned int cov_index) {
    
    std::vector<bool> produced_measurement_source_indicies( num_sources_produced_measurements_,false);

    unsigned int tmp = 0x1;

    for (int ii = 0; ii < num_sources_produced_measurements_; ++ii ) {

        if( cov_index & (tmp << ii) ) {
            produced_measurement_source_indicies[ii] = true;
        }
    }

    return produced_measurement_source_indicies;
}


//---------------------------------------------------------------------------------------------------------

template<typename _Model>  
void CentralizedMeasurementFusion<_Model>::PerformUpdateFromSingleSensor(const SourceContainer& source_container, const std::vector<Measurement>& meas, Model& model){



MatH H = source_container.GetLinObsMatState(meas.front().source_index, model.state_, meas.front().transform_state, meas.front().transform_data_t_m);      // Jacobian of observation function w.r.t. state
MatV V = source_container.GetLinObsMatSensorNoise(meas.front().source_index, model.state_, meas.front().transform_state, meas.front().transform_data_t_m);                             // Jacobian of observation function w.r.t. noise
Eigen::MatrixXd K;                                                                    // Kalman Gain
MatS S_inverse;                                                            // Innovation covariance inverse
Eigen::MatrixXd nu_i;                                                                 
Eigen::MatrixXd nu(V.rows(),1);                                                       // Total innovation term
nu.setZero();
Eigen::MatrixXd covSum(V.rows(),V.rows());
covSum.setZero();
VecCov mu_minus;

unsigned int source_index = meas.front().source_index;

Measurement estimated_meas = source_container.GetEstMeas(meas.front().source_index, model.state_, meas.front().transform_state, meas.front().transform_data_t_m); 

S_inverse = (H*model.err_cov_*H.transpose() + V*source_container.GetParams(source_index).meas_cov_ *V.transpose()).inverse();
K = model.err_cov_*H.transpose()*S_inverse;

// std::cout << "K " << std::endl << K << std::endl;

DataType B0 = 1;

// Get total weighted innovation and part of the cov_tilde
for (Measurement m : meas) {
    nu_i = source_container.OMinus(m.source_index,m,estimated_meas);
    nu += m.weight*nu_i;
    covSum+= m.weight*nu_i*nu_i.transpose();
    B0 -= m.weight;
    // std::cout << "weight: " << m.weight << std::endl;
    // std::cout << "nu_i: " << std::endl << nu_i << std::endl;
}


// std::cout << "nu" << std::endl << nu << std::endl;

// Finish constructing cov_sum
covSum -= nu*nu.transpose();


// construct covariance. It has been verified
model.err_cov_ = model.err_cov_+ K*(covSum*K.transpose() -(1-B0)*H*model.err_cov_); 
mu_minus = K*nu;


// Reset mean to zero   
model.OPlusEQ(mu_minus);

MatModelCov Jr = Model::Jr(mu_minus);
// std::cout << "Jr: " << std::endl << Jr << std::endl;
model.err_cov_ = Jr * model.err_cov_ *Jr.transpose();


// std::cout << "H: " << std::endl << H << std::endl;
// std::cout << "S: " << std::endl << S_inverse << std::endl;
// std::cout << "K: " << std::endl << K << std::endl;


}





} // namespace rransac

#endif //RRANSAC_COMMON_MODELS_CENTRALIZED_MEASUREMENT_FUSION_H_