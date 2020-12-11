#ifndef RRANSAC_COMMON_DATA_ASSOCIATION_DATA_ASSOCIATION_HOST_H_
#define RRANSAC_COMMON_DATA_ASSOCIATION_DATA_ASSOCIATION_HOST_H_

#include "system.h"

namespace rransac
{

/**
 * \class DataAssociationHost
 * 
 * Associates the new measurements System<Model>::new_measurements_ to the models, clusters, and data tree. 
 * This class is a host to two policies: model data association policy and cluster and data tree data association policy. 
 *  
 * The model data association is a policy that must associate the new measurements to models and assign model associated measurement a weight even if the
 * weight is 1. If a measurement is associated to a model, it is removed from System<Model>::new_meas_ and added to Model<State>::new_assoc_meas_. 
 * It must also update the model member variable Model<State>::model_likelihood_update_info_ with the proper information.
 * The policy must expose the function static void PolicyDataAssociationModel(System<Model>& sys) which is called by the host class.
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
     * Associates new measurements to the models, clusters, and data tree
     */ 
    static void AssociateNewMeasurements(System<tModel>& sys ) {
        DataAssociationModel(sys);
        DataAssociationClusterDataTree(sys);
    }

private:

    /**
     * Calls the policy class member function void PolicyDataAssociationModel(System<Model>& sys) to
     * associate new measurements to the model.
     */ 
    static void DataAssociationModel(System<tModel>& sys){
        DataAssociationHost::PolicyDataAssociationModel(sys);
    }

    /**
     * Calls the policy class memver function void PolicyDataAssociationClusterDataTree(System<tModel>& sys)
     * to associate new measurements to clusters and the data tree
     */
    static void  DataAssociationClusterDataTree(System<tModel>& sys) {
        DataAssociationHost::PolicyDataAssociationClusterDataTree(sys);
    }

};



} // namespace rransac



#endif // RRANSAC_COMMON_DATA_ASSOCIATION_DATA_ASSOCIATION_HOST_H_