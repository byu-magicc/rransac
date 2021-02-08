#ifndef RRANSAC_COMMON_DATA_ASSOCIATION_DATA_ASSOCIATION_HOST_H_
#define RRANSAC_COMMON_DATA_ASSOCIATION_DATA_ASSOCIATION_HOST_H_
#pragma once



#include "rransac/system.h"

namespace rransac
{

/**
 * \class DataAssociationHost
 * 
 * Associates the new measurements System::new_measurements_ to the tracks, clusters, and data tree. 
 * This class is host to two policies: the model data association policy and the cluster and data tree data association policy. 
 *  
 * The model data association is a policy that must associate the new measurements to tracks and assign model associated measurement a weight. The sum of
 * the weights of associated measurement from every unique measurement source must be 1. For example: if a model is associated with five measurements, three 
 * from one source and two from another, then the wights of the three measurements from the first source must sum to one and the weights of the 
 * two measurements from the other source must sum to one. If a measurement is associated to a model, it is removed from System::new_meas_ and added to ModelBase::new_assoc_meas_
 * using the method ModelBase::AddNewMeasurement. It must also update the model member variable Model::model_likelihood_update_info_ with the proper information.
 * The policy must expose the function static void PolicyDataAssociationModel(System& sys) which is called by the host class.
 * 
 * The model data association policy must also expose the policy CalculateMeasurmentAndLikelihoodDataPolicy which calculates the weights for the measurements and the model likelihood.
 * This policy is not used by this class. but it is used by RANSAC. 
 * 
 * 
 * The cluster and data tree policy associates new measurements to clusters and the data tree. When measurements are associated, they are removed from 
 * System::new_measurements_ and added to the corresponding cluster and data tree. The policy class must expose a function 
 * static void PolicyDataAssociationClusterDataTree(System<tModel>& sys).
 */ 
template<typename tModel, template<class> typename tModelDataAssociationPolicyClass, template<class> typename tClusterDataTreeAssociationPolicyClass>
class DataAssociationHost : public tModelDataAssociationPolicyClass<tModel>, tClusterDataTreeAssociationPolicyClass<tModel>{

public:

    /**
     * Associates new measurements to the tracks, clusters and the data tree. 
     * @param[in,out] sys The object that contains all of the data of RRANSAC. This includes the new measurements, tracks, clusters and data tree. 
     */ 
    static void AssociateNewMeasurements(System<tModel>& sys ) {
        DataAssociationModel(sys);
        DataAssociationClusterDataTree(sys);

        if(sys.new_meas_.size() != 0)   
            throw std::runtime_error("DataAssociationHost::AssociateNewMeasurements: All new measurements should have been copied to a model, cluster, or data tree and removed from System<Model>::new_meas_");
    }

private:

    /**
     * Calls the policy class member function void PolicyDataAssociationModel to
     * associate new measurements to the model.
     * @param[in,out] sys The object that contains all of the data of RRANSAC. This includes the new measurements, tracks, clusters and data tree. 
     */ 
    static void DataAssociationModel(System<tModel>& sys){
        DataAssociationHost::PolicyDataAssociationModel(sys);
    }

    /**
     * Calls the policy class memver function void PolicyDataAssociationClusterDataTree
     * to associate new measurements to clusters and the data tree
     * @param[in,out] sys The object that contains all of the data of RRANSAC. This includes the new measurements, tracks, clusters and data tree.
     */
    static void  DataAssociationClusterDataTree(System<tModel>& sys) {
        DataAssociationHost::PolicyDataAssociationClusterDataTree(sys);
    }

};



} // namespace rransac



#endif // RRANSAC_COMMON_DATA_ASSOCIATION_DATA_ASSOCIATION_HOST_H_