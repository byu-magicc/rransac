#ifndef RRANSAC_DATA_CONTAINERS_DATA_TREE_BASE_H_
#define RRANSAC_DATA_CONTAINERS_DATA_TREE_BASE_H_

#include "common/measurement/measurement_base.h"
#include <list>

namespace rransac
{
    
/**
 * \class DataTreeBase
 * This container holds all of the measurements that are not associated with a model or cluster. The data object is
 * std::list<tContainer> where tContainer is another container with an iterator object defined. 
 * 
 */

template <typename tContainer, typename tTransform, typename tDerived>
class DataTreeBase {

typedef tContainer Container;
typedef tTransform transform;
typedef tDerived Derived;
   
public:

/**
 * \struct IteratorPair
 * Contains a match of iterators. The outer iterator pertains to the outer container for the data object, and
 * the inner container pertains to the inner container for the data object. Note that the data object is type
 * std::list<tContainer>. So the outer iterator points to an element in the list and the inner iterator points to and 
 * element of tContainer.
 */ 
struct IteratorPair {
    typename std::list<tContainer>::iterator& outer_it;
    typename tContainer::iterator& inner_it;
};

/**
 * Adds a measurement to the data tree. It is assumed that the measurement has a valid time stamp and valid data.
 * @param[in] meas The measurement to be added.
 */
void AddMeas(const Meas& meas) {
    ++size_;
    static_cast<tDerived*>(this)->DerivedAddMeas(meas);
} 

/**
 * Adds measurements from a container that has an iterator object. The measurements can have different time stamps.
 * It is assumed that the measurements have a valid time stamp and valid data. This function will call AddMeas for every
 * measurement.
 */ 
template <typename tContainerMeas>
void AddMeasurements(const tContainerMeas& measurements ){
    for (auto iter = measurements.begin(); iter != measurements.end(); ++iter) {
        AddMeas(*iter);
    }
}

/**
 * Removes a measurement from the data tree using iterators defined in iter_pair
 * @param[in] iter_pair Contains the matched iterator pair necessary to remove an element. * 
 */
void RemoveMeas(const IteratorPair& iter_pair){
    static_cast<tDerived*>(this)->DerivedRemoveMeas(iter_pair);
} 

/**
 * Removes all of the measurements indicated by iter_pairs by calling RemoveMeas
 * @param[in] iter_pair Contains all of the matched iterator pairs that indicate which measurements to remove.
 */ 
template <typename tContainerIteratorPair>
void RemoveMeasurements(const tContainerIteratorPair& iter_pairs){
    for(auto iter = iter_pairs.begin(); iter != iter_pairs.end(); ++iter){
        RemoveMeas(*iter);
    }
}

/**
 * Attempts to find a cluster from the measurements in the data tree. If a cluster was not found, the size of the vector container of iterator pairs will be 
 * zero.
 * @param[in] params The system parameters. It contains parameters that define the cluster.
 * @param[out] iter_pairs Contains iterators to the elements that make up the cluster.
 * @return returns true if a cluster is found.
 */
bool FindCluster(const Parameters& params, std::vector<IteratorPair>& iter_pairs) const {
    return static_cast<tDerived*>(this)->DerivedFindCluster(params, iter_pairs);
}

/**
 * Removes all of the measurements with a time stamp before or equal to
 * the expiration time. 
 */ 
void PruneDataTree(const double expiration_time);

/**
 * Transforms the measurements using the transform provided. 
 */ 
void TransformMeasurements(const tTransform& transform);

unsigned int Size() {return size_;};

std::list<tContainer> data_;

private:
DataTreeBase() : size_{0} {}
~DataTreeBase()=default;
friend tDerived;



unsigned long int size_;


};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                            Definitions
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename tContainer, typename tTransform, typename tDerived>
void DataTreeBase<tContainer,tTransform, tDerived>::PruneDataTree(const double expiration_time) {

    auto iter = data_.begin();
    while(iter != data_.end() && iter->begin()->time_stamp <= expiration_time) {
        size_ -= iter->size();
        iter = data_.erase(iter); // Erases the element and returns an iterator to the next element which is now the first element of the list
    }
   
}

// -----------------------------------------------------------------------------------------------------------------

template <typename tContainer, typename tTransform, typename tDerived>
void DataTreeBase<tContainer,tTransform, tDerived>::TransformMeasurements(const tTransform& transform) {

    if (!transform.transform_null_) {
        for (auto outer_iter = data_.begin(); outer_iter != data_.end(); ++outer_iter) {
            for(auto inner_iter = outer_iter->begin(); inner_iter != outer_iter->end(); ++ inner_iter) {
                transform.TransformMeasurement(*inner_iter);
            }
        }
    }

}

} // namespace rransac



#endif // RRANSAC_DATA_CONTAINERS_DATA_TREE_BASE_H_