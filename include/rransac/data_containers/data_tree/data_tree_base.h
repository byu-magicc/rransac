#ifndef RRANSAC_DATA_CONTAINERS_DATA_TREE_BASE_H_
#define RRANSAC_DATA_CONTAINERS_DATA_TREE_BASE_H_
#pragma once

#include <list>

#include "rransac/common/measurement/measurement_base.h"
#include "rransac/system.h"


namespace rransac
{
    
/**
 * \class DataTreeBase
 * This class is the base class for the containers that hold all of the measurements that are not associated with a track or cluster. 
 * This base class uses the curiously recurring template parameter design methodoology to specify the API for all derived classes. 
 * The object type of the data is std::list<tContainer> where tContainer is another container such ast std::list or std::vector with 
 * an iterator object defined. 
 * 
 */

template <typename _TreeDataType, typename tDataType, typename _TransformDataType, template<typename, typename> typename tDerived>
class DataTreeBase {

typedef _TreeDataType TreeDataType;    /**< The object type of the data. */
typedef tDerived<tDataType,_TransformDataType> Derived;   /**< The derived class. */
typedef tDataType DataType;            /**< The scalar object for the data. Ex. float, double, etc. */
typedef _TransformDataType TransformDataType;
typedef Meas<DataType,TransformDataType> Measurement;

public:

/**
 * Adds a measurement to the data tree. It is assumed that the measurement has a valid time stamp, type and data.
 * @param[in] sys The object that contains all of the data of RRANSAC. Thus it contains all of the measurements and data tree. 
 * @param[in] meas The measurement to be added.
 */
template <typename tSystem>
void AddMeasurement(const tSystem& sys, const Measurement& meas) {
    static_cast<Derived*>(this)->DerivedAddMeasurement(sys, meas);
} 

/**
 * Adds measurements from a container that has an iterator object such as std::list or std::vector. The measurements can have different time stamps.
 * It is assumed that the measurements have a valid time stamp and valid data. This function will call AddMeasurement for every
 * measurement.
 * @param[in] sys The object that contains all of the data of RRANSAC. Thus it contains all of the measurements and data tree. 
 * @param[in] measurements The measurements to be added.
 */ 
template <typename tSystem, typename tContainerMeas>
void AddMeasurements(const tSystem& sys, const tContainerMeas& measurements ){
    for (auto iter = measurements.begin(); iter != measurements.end(); ++iter) {
        AddMeasurement(sys, *iter);
    }
}

/**
 * Removes a measurement from the data tree using meas_info. The object type of 
 * the measurement location info is specified by the derived class. 
 * @param[in] meas_info Contains the information necessary to remove a measurement. The necessary information is determined by the derived class. 
 */
template <typename tMeasurementLocationInfo>
void RemoveMeasurement(const tMeasurementLocationInfo& meas_info){
    static_cast<Derived*>(this)->DerivedRemoveMeasurement(meas_info);
} 

/**
 * Removes all of the measurements indicated by container_meas_info by calling RemoveMeasurement for each measurement. 
 * The object type of the measurement location info is specified by the derived class. 
 * @param[in] container_meas_info A container that contains information to remove multiple measurements. The necessary information is determined by the derived class.
 */ 
template <typename tContainerMeasurementLocationInfo>
void RemoveMeasurements(const tContainerMeasurementLocationInfo& container_meas_info){
    for(auto iter = container_meas_info.begin(); iter != container_meas_info.end(); ++iter){
        RemoveMeasurement(*iter);
    }
}

/**
 * Attempts to find a cluster from the measurements in the data tree. If a cluster is found, it will 
 * construct a cluster and add it to System::cluster_. A cluster must have enough measurements to form a minimum 
 * subset as defined by Parameters::RANSAC_minimum_subset_
 * @param[in,out] sys The complete system information
 * @return returns true if a cluster is found.
 */
template <typename tSystem>
void ConstructClusters(tSystem& sys) {
    static_cast<Derived*>(this)->DerivedConstructClusters(sys);
}

/**
 * Removes all of the measurements with a time stamp before or equal to
 * the expiration time. 
 * @param[in] sys The complete system information. Thus it contains all of the measurements.
 * @param[in] expiration_time Measurements before or equal to this time will be removed from the data set.
 */ 
template <typename tSystem>
void PruneDataTree(const tSystem& sys, const double expiration_time) {
    static_cast<Derived*>(this)->DerivedPruneDataTree(sys, expiration_time);
}

/**
 * Transforms the measurements using the transform provided. 
 * @param[in] T The transformation object provided by the user. The object should already have the data it needs to transform the model.
 * @see TransformBase
 */ 
template< typename tTransform>
void TransformMeasurements(const tTransform& transform) {
    static_cast<Derived*>(this)->DerivedTransformMeasurements(transform);
}

unsigned int Size() const {return size_;}; /**< The number of measurements in the data tree. */

TreeDataType data_; /**< The data that contains all of the measurements on the data tree. */

private:
DataTreeBase()=default;
~DataTreeBase()=default;
friend Derived;            /**< The derived class. */



unsigned long int size_=0;   /**< The number of measurements in the data tree. */


};






} // namespace rransac



#endif // RRANSAC_DATA_CONTAINERS_DATA_TREE_BASE_H_