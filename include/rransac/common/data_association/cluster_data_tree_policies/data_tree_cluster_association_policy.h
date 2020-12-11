#ifndef RRANSAC__COMMON__DATA_ASSOCIATION__DATA_CLUSTER_DATA_TREE_POLICIES__DATA_TREE_CLUSTER_ASSOCIATION_POLICY_H_
#define RRANSAC__COMMON__DATA_ASSOCIATION__DATA_CLUSTER_DATA_TREE_POLICIES__DATA_TREE_CLUSTER_ASSOCIATION_POLICY_H_

#include "system.h"

namespace rransac
{

/**
 * \class
 * This policy assumes that the data tree is of type DataTreeClusters. The remaining new measurements are added
 * the the data tree and then all of the new measurements are removed from System<tModel>::new_meas_
 * 
 */ 

template<typename tModel>
class DataTreeClusterAssociationPolicy {

public:

    static void  PolicyDataAssociationClusterDataTree(System<tModel>& sys) {
        for(auto outer_iter = sys.new_meas_.begin(); outer_iter != sys.new_meas_.end(); ++outer_iter) {
            for(auto inner_iter = outer_iter->begin(); inner_iter != outer_iter->end(); ++inner_iter) {
                sys.data_tree_.AddMeasurement(sys, *inner_iter);
            }
        }
        sys.new_meas_.clear();
    }



};

} // namespace rransac



#endif // RRANSAC__COMMON__DATA_ASSOCIATION__DATA_CLUSTER_DATA_TREE_POLICIES__DATA_TREE_CLUSTER_ASSOCIATION_POLICY_H_