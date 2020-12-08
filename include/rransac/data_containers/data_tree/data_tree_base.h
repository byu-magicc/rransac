#ifndef RRANSAC_DATA_CONTAINERS_DATA_TREE_BASE_H_
#define RRANSAC_DATA_CONTAINERS_DATA_TREE_BASE_H_

#include "common/measurement/measurement_base.h"
#include "system.h"
#include <list>

namespace rransac
{
    
/**
 * \class DataTreeBase
 * This container holds all of the measurements that are not associated with a model or cluster. The data object is
 * std::list<tContainer> where tContainer is another container with an iterator object defined. 
 * 
 */

template <typename tData, typename tDerived>
class DataTreeBase {

typedef tData Data;
typedef tDerived Derived;
   
public:

/**
 * Adds a measurement to the data tree. It is assumed that the measurement has a valid time stamp, type and valid data.
 * @param[in] meas The measurement to be added.
 */
template <typename tSystem>
void AddMeasurement(const tSystem& sys, const Meas& meas) {
    static_cast<tDerived*>(this)->DerivedAddMeasurement(sys, meas);
} 

/**
 * Adds measurements from a container that has an iterator object. The measurements can have different time stamps.
 * It is assumed that the measurements have a valid time stamp and valid data. This function will call AddMeasurement for every
 * measurement.
 */ 
template <typename tSystem, typename tContainerMeas>
void AddMeasurements(const tSystem& sys, const tContainerMeas& measurements ){
    for (auto iter = measurements.begin(); iter != measurements.end(); ++iter) {
        AddMeasurement(sys, *iter);
    }
}

/**
 * Removes a measurement from the data tree using meas_info
 * @param[in] meas_info Contains the information necessary to remove a measurement 
 */
template <typename tMeasurementLocationInfo>
void RemoveMeasurement(const tMeasurementLocationInfo& meas_info){
    static_cast<tDerived*>(this)->DerivedRemoveMeasurement(meas_info);
} 

/**
 * Removes all of the measurements indicated by iter_pairs by calling RemoveMeasurement
 * @param[in] container_meas_info A container that contains information to remove multiple measurements
 */ 
template <typename tContainerMeasurementLocationInfo>
void RemoveMeasurements(const tContainerMeasurementLocationInfo& container_meas_info){
    for(auto iter = container_meas_info.begin(); iter != container_meas_info.end(); ++iter){
        RemoveMeasurement(*iter);
    }
}

/**
 * Attempts to find a cluster from the measurements in the data tree. If a measurement is found, it will 
 * construct a cluster and add it to the System. A cluster must have enough measurements to form a minimum 
 * subset as defined by Parameters::RANSAC_minimum_subset_
 * @param[in,out] sys The complete system information
 * @return returns true if a cluster is found.
 */
template <typename tSystem>
void ConstructClusters(tSystem& sys) {
    return static_cast<tDerived*>(this)->DerivedConstruct(sys);
}

/**
 * Removes all of the measurements with a time stamp before or equal to
 * the expiration time. 
 * @param[in] sys The complete system information
 * @param[in] expiration_time measurements before or equal to this time will be removed from the data set
 */ 
template <typename tSystem>
void PruneDataTree(const tSystem& sys, const double expiration_time) {
    static_cast<tDerived*>(this)->DerivedPruneDataTree(sys, expiration_time);
}

/**
 * Transforms the measurements using the transform provided. 
 */ 
template< typename tTransform>
void TransformMeasurements(const tTransform& transform) {
    static_cast<tDerived*>(this)->DerivedTransformMeasurements(transform);
}

unsigned int Size() {return size_;};

tData data_;

private:
DataTreeBase()=default;
~DataTreeBase()=default;
friend tDerived;



unsigned long int size_=0;


};






} // namespace rransac



#endif // RRANSAC_DATA_CONTAINERS_DATA_TREE_BASE_H_