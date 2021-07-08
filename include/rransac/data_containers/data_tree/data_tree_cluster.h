#ifndef RRANSAC_DATA_CONTAINERS_DATA_TREE_CLUSTER_H_
#define RRANSAC_DATA_CONTAINERS_DATA_TREE_CLUSTER_H_
#pragma once

#include "rransac/data_containers/data_tree/data_tree_base.h"
#include "rransac/data_containers/cluster.h"

namespace rransac
{

/**
 * \class DataTreeClusters
 * This data tree organizes the measurements in clusters. A cluster is a group of neighboring measurements. For more information on a cluster, see the class Cluster. 
 * Since track initialization works on clusters to improve speed performance, one technique to organize the measurements is in clusters. Since the measurements are already 
 * in clusters, we do not have to construct clusters which saves times. 
 */ 
template<typename tDataType, typename _TransformDataType>
class DataTreeClusters : public DataTreeBase<std::list<Cluster<tDataType,_TransformDataType>>, tDataType, _TransformDataType, DataTreeClusters> {

public:

typedef DataTreeBase<std::list<Cluster<tDataType,_TransformDataType>>, tDataType, _TransformDataType, DataTreeClusters> Base;
typedef tDataType DataType;  /**< The scalar object for the data. Ex. float, double, etc. */
typedef typename Base::TransformDataType TransformDataType;
typedef typename Base::Measurement Measurement;
typedef typename Base::TreeDataType TreeDataType;

/**
 * \struct MeasurementLocationInfo
 * Contains the location information of a measurement. This includes an iterator to the cluster it pertains to, and the Cluster::IteratorPair to locate
 * the measurement in the cluster.
 */ 
struct MeasurementLocationInfo {

    typename TreeDataType::iterator cluster_iter;  /**< An iterator to the cluster a measurement is in. */
    typename Cluster<DataType,TransformDataType>::IteratorPair iter_pair;            /**< The information needed to locate the measurement in the cluster specified by MeasurementLocationInfo::cluster_iter. */

};

/**
 * Adds a measurement to the data tree. It is assumed that the measurement has a valid time stamp, type and data. A measurement is added by first
 * checking to see if it is a neighbor to an existing cluster, if it is not, a new cluster will be formed with the measurement. 
 * @param[in] sys The object that contains all of the data of RRANSAC. Thus it contains all of the measurements and data tree. 
 * @param[in] meas The measurement to be added.
 */
template <typename tSystem>
void DerivedAddMeasurement(const tSystem& sys, const Measurement& meas);

/**
 * Removes a measurement from the data tree using meas_info. If the measurement to be removed is the only measurement in the cluster, the cluster will
 * also be removed. 
 * @param[in] meas_info Contains the information necessary to remove a measurement. 
 */
void DerivedRemoveMeasurement(const MeasurementLocationInfo& meas_info); 

/**
 * The data tree has a list of clusters. This function finds every cluster that meets the minimum subset requirements as defined by
 * Parameters::RANSAC_minimum_subset_ and gives a reference to System::clusters_ to later be used by the track initializer.
 * @param[in,out] sys The complete system information.
 * @return returns true if a cluster is found.
 */
template <typename tSystem>
void DerivedConstructClusters(tSystem& sys);

/**
 * Removes all of the measurements with a time stamp before or equal to
 * the expiration time. After removing the measurements, if a cluster is empty, the cluster is also removed. 
 * @param[in] sys The complete system information.
 * @param[in] expiration_time measurements before or equal to this time will be removed from the data set.
 */ 
template <typename tSystem>
void DerivedPruneDataTree(const tSystem sys, const double expiration_time);

/**
 * Transforms the measurements using the transform provided. 
 * @param[in] transform The transformation object provided by the user. The object should already have the data it needs to transform the model.
 * @see TransformBase
 */ 
template< typename tTransform>
void DerivedTransformMeasurements(const tTransform& transform) {
    for(auto iter = this->data_.begin(); iter!= this->data_.end(); ++iter)
        iter->TransformMeasurements(transform);
}

private:

/**
 * Adds all of the measurements from iter2 to iter1 and removes iter2
 * @param[in] iter1 An iterator to a cluster;
 * @param[in] iter2 An iterator to a cluster
 */ 
void MergeClusters(const typename TreeDataType::iterator iter1, const typename TreeDataType::iterator iter2);



};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                            Definitions
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename tDataType, typename _TransformDataType>
template<typename tSystem>
void DataTreeClusters<tDataType,_TransformDataType>::DerivedAddMeasurement(const tSystem& sys, const Measurement& meas) {

    std::vector<typename TreeDataType::iterator> neighbor_clusters;

    // Search through all of the clusters and see which ones the measurement is a neighbor to
    for(auto iter = this->data_.begin(); iter != this->data_.end(); ++iter) {

        // see if the measurement is a neighbor
        if (iter->IsNeighboringMeasurement(sys.source_container_,0, sys.params_, meas)) {
            neighbor_clusters.push_back(iter);
        }
    }
    
    if (neighbor_clusters.size() == 0) {           // The measurement is not a neighbor to any cluster, so make a new one
        this->data_.emplace_back(meas);
    } else if (neighbor_clusters.size() == 1) {    // The measurement is only a neighbor to one cluster, so add it

        neighbor_clusters.front()->AddMeasurement(meas);
    } else {                                       // The measurement is a neighbor to multiple clusters, so add it to one and merge the others

        auto iter_first = neighbor_clusters.begin();
        (*iter_first)->AddMeasurement(meas);
        for(auto iter_other = std::next(iter_first); iter_other != neighbor_clusters.end(); ++ iter_other) {
            MergeClusters(*iter_first, *iter_other);
        }
    }

    ++this->size_;
}

//-----------------------------------------------------------------------------------------------------------------------
template<typename tDataType, typename _TransformDataType>
void DataTreeClusters<tDataType,_TransformDataType>::MergeClusters(const typename TreeDataType::iterator iter1, const typename TreeDataType::iterator iter2) {

    for(auto outer_iter = iter2->data_.begin(); outer_iter != iter2->data_.end(); ++outer_iter ) {
        for (auto inner_iter = outer_iter->begin(); inner_iter != outer_iter->end(); ++inner_iter) {
            iter1->AddMeasurement(*inner_iter);
        }
    }

    this->data_.erase(iter2);

}

//-----------------------------------------------------------------------------------------------------------------------
template<typename tDataType, typename _TransformDataType>
void DataTreeClusters<tDataType,_TransformDataType>::DerivedRemoveMeasurement(const MeasurementLocationInfo& meas_info) {

    meas_info.cluster_iter->RemoveMeasurement(meas_info.iter_pair);
    if (meas_info.cluster_iter->Size() == 0)
        this->data_.erase(meas_info.cluster_iter);

    --this->size_;
}

//-----------------------------------------------------------------------------------------------------------------------
template<typename tDataType, typename _TransformDataType>
template<typename tSystem>
void DataTreeClusters<tDataType,_TransformDataType>::DerivedConstructClusters(tSystem& sys) {

    sys.clusters_.clear();
    for(auto cluster_iter = this->data_.begin(); cluster_iter != this->data_.end(); ++cluster_iter) {
        if(cluster_iter->data_.size() >= sys.params_.cluster_min_size_requirement_) {
            if (cluster_iter->cluster_label_ < 0) {
                cluster_iter->cluster_label_ = sys.cluster_label_;
                sys.cluster_label_++;
            }
            sys.clusters_.push_back(cluster_iter); }

    }

}

//-----------------------------------------------------------------------------------------------------------------------

template<typename tDataType, typename _TransformDataType>
template <typename tSystem>
void DataTreeClusters<tDataType,_TransformDataType>::DerivedPruneDataTree(const tSystem sys, const double expiration_time) {
    auto iter = this->data_.begin();
    while(iter != this->data_.end()) {
        double size_diff = iter->Size();
        iter->PruneCluster(expiration_time);
        size_diff -= iter->Size();
        this->size_ -= size_diff;

        // After pruning the cluster, if it doesn't have any elements, remove it
        if(iter->Size() == 0) {
            iter = this->data_.erase(iter); // removes the element at iter and returns the next iterator
        } 
        // If the latest measurement in the cluster has a time stamp before the current time minus the time threshold,
        // It won't receive any other measurements and it hasn't been made into a model. Thus it won't ever become a model
        // so it should be removed. 
        else if (iter->data_.back().front().time_stamp <= (sys.current_time_ - sys.params_.cluster_time_threshold_)) {
            this->size_ -= iter->Size();
            iter = this->data_.erase(iter); // removes the element at iter and returns the next iterator
        } else {
            ++iter; 
        }
 
    }
}

} // namespace rransac

#endif //RRANSAC_DATA_CONTAINERS_DATA_TREE_CLUSTER_H_