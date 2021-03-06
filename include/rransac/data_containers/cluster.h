#ifndef RRANSAC_DATA_CONTAINERS_CLUSTER_BASE_H_
#define RRANSAC_DATA_CONTAINERS_CLUSTER_BASE_H_
#pragma once

#include <list>

#include <Eigen/Core>

#include "rransac/common/measurement/measurement_base.h"
#include "rransac/common/sources/source_container.h"

namespace rransac
{


/**
 * \class Cluster
 * This data structure is used to contain neighboring measurements. Two measurements of the same time stamp are considered neighbors 
 * if the distance between them is withing the threshold Parameters::cluster_position_threshold_, or if the distance in pose 
 * normalized by the time difference between the two measurements are within the threshold Parameters::cluster_velocity_threshold_ and 
 * distance of their time stamps is within Parameters::cluster_time_threshold_.
 * 
 * The data is organized is a list of lists, i.e. std::list<std::list<Measurement>. The outer list organized the measurements in chronological 
 * order, and the inner list contains all of the measurements given to the cluster with the same time stamp.  
 */ 

template<typename _DataType, typename _TransformDataType>    
class Cluster {

public:

typedef _TransformDataType TransformDataType;          /**< The transformation data type. */
typedef _DataType DataType;                            /**< The scalar object for the data. Ex. float, double, etc. */
typedef Meas<DataType,TransformDataType>  Measurement; /**< The measurement type. */

/**
 * \struct ConstIteratorPair
 * Contains a match of iterators that cannot change the data they point to. 
 * The outer iterator pertains to the outer container for the data object, and
 * the inner container pertains to the inner container for the data object. So the outer iterator points to a
 * list, and the inner iterator points to a measurement.
 */ 
struct ConstIteratorPair {
    typename std::list<std::list<Measurement>>::const_iterator outer_it;
    typename std::list<Measurement>::const_iterator inner_it;
};

/**
 * \struct IteratorPair
 * Contains a match of iterators. The outer iterator pertains to the outer container for the data object, and
 * the inner container pertains to the inner container for the data object. So the outer iterator points to a
 * list, and the inner iterator points to a measurement.
 */ 
struct IteratorPair {
    typename std::list<std::list<Measurement>>::iterator outer_it;
    typename std::list<Measurement>::iterator inner_it;
};

/**
 * Default constructor.
 */ 
Cluster() : cluster_label_(-1), size_(0) {}

/**
 * A cluster can be constructed with a measurement.
 * @param[in] meas The first measurement to be added to the cluster.
 */ 
explicit Cluster(const Measurement& meas) : Cluster() { AddMeasurement(meas); }

/**
 * The copy constructor
 */ 
Cluster(const Cluster& other) : cluster_label_(other.cluster_label_), size_(other.size_), data_(other.data_) {}

/**
 * The assignment constructor
 */ 
Cluster& operator=(const Cluster& other) {
    cluster_label_ = other.cluster_label_;
    size_ = other.size_;
    data_ = other.data_;
}

/**
 * Adds a measurement to the cluster. It is assumed that the measurement has a valid time stamp and valid data.
 * @param[in] meas The measurement to be added.
 */
void AddMeasurement(const Measurement& meas);

/**
 * Adds measurements from a container that has an iterator object. For example, a list or vector of of measurements.
 * It is assumed that the measurements have a valid time stamp and valid data. 
 * The measurements can have different time stamps.This function will call AddMeasurement for every
 * measurement.
 * @param[in] measurements The measurements to be added. 
 */ 
template<typename _ContainerMeas>
void AddMeasurements(const _ContainerMeas& measurements){    
    for (auto iter = measurements.begin(); iter != measurements.end(); ++iter) {
        AddMeasurement(*iter);
    }
}


/**
 * Removes a measurement from the data tree using an object of type IteratorPair.
 * @param[in] iter_pair Contains the matched iterator pair necessary to remove an element. 
 */
void RemoveMeasurement(const IteratorPair& iter_pair);


/**
 * Removes multiple measurements by calling RemoveMeasurement for each measurement. 
 * @param[in] iter_pair Contains all of the matched iterator pairs that indicate which measurements to remove.
 */ 
template<typename _ContainerIteratorPair>
void RemoveMeasurements(const _ContainerIteratorPair& iter_pairs){
    for(auto iter = iter_pairs.begin(); iter != iter_pairs.end(); ++iter) {
        RemoveMeasurement(*iter);
    }
}

/**
 * Removes all of the measurements with a time stamp before or equal to
 * the expiration time.
 * @param[in] expiration_time The expiration time stamp.
 */ 
void PruneCluster(const double expiration_time);

/**
 * Transforms the measurements using the transform provided. 
 * The transform provided must have a member function called
 * TransformMeasurement.
 * @param[in] transform An object of the transformation class
 * @see TransformBase
 */ 
template< typename _Transform>
void TransformMeasurements(const _Transform& transform);

/**
 * Returns the number of measurements in the cluster
 */
unsigned int Size() const {return size_;};

/**
 * Returns true if the measurement is a neighboring measurement to a measurement in the cluster.
 * Two measurements of the same time stamp are considered neighbors 
 * if the distance between them is withing the threshold Parameters::cluster_position_threshold_, or if the distance in pose 
 * normalized by the time difference between the two measurements are within the threshold Parameters::cluster_velocity_threshold_ and 
 * distance of their time stamps is within Parameters::cluster_time_threshold_.
 * @param[in] source The measurement source that produced the measurement.
 * @param[in] param The system parameters.
 * @param[in] meas The measurement that is tested to see if it is a neighboring measurement.
 * @see SourceBase
 */ 
template<typename _SourceContainer>
bool IsNeighboringMeasurement(const _SourceContainer& source_container, const unsigned int source_index, const Parameters& params, const Measurement& meas) const;


std::list<std::list<Measurement>> data_; /**< Contains all of the measurements. The outer container organizes the measurements in chronological order. */

int64_t cluster_label_;  /**< When a cluster is elevated to a good cluster, it will receive a unique label whose numerical value is non negative. */ 

private:

unsigned int size_; /**< The total number of measurements in the cluster. */

};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                            Definitions
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename _DataType,typename _TransformDataType> 
void Cluster<_DataType,_TransformDataType>::AddMeasurement(const Measurement& meas) {

++size_;

// There are no measurements. Just add it.
if (this->data_.begin() == this->data_.end()) {
    this->data_.emplace_back(std::list<Measurement>{meas});
} 
// All measurements occurred before the new one. So add it to the back.
else if (this->data_.back().begin()->time_stamp < meas.time_stamp) {
    this->data_.emplace_back(std::list<Measurement>{meas});
} 
// The new measurement occurred before all the other measurements
else if (this->data_.front().begin()->time_stamp > meas.time_stamp) {
    this->data_.emplace_front(std::list<Measurement>{meas});
}
else {
    // Search from the end to the beginning to find where to place it
    for (auto iter = this->data_.rbegin(); iter != this->data_.rend(); ++iter) {
        
        

        if ((*iter).begin()->time_stamp == meas.time_stamp) {
            (*iter).push_back(meas);                                  
            break;
        } else if ((*iter).begin()->time_stamp < meas.time_stamp) {
            this->data_.insert(iter.base(),std::list<Measurement>{meas});
            break;
        } 
    }
}

}

//--------------------------------------------------------------------------------------------
template<typename _DataType,typename _TransformDataType> 
void Cluster<_DataType,_TransformDataType>::RemoveMeasurement(const IteratorPair& iter_pair) {

    --size_;

    iter_pair.outer_it->erase(iter_pair.inner_it);
    if (iter_pair.outer_it->size() == 0) {
        this->data_.erase(iter_pair.outer_it);
    }
}

//---------------------------------------------------------------------------------------------
template<typename _DataType,typename _TransformDataType> 
void Cluster<_DataType,_TransformDataType>::PruneCluster(const double expiration_time) {

    auto iter = data_.begin();
    while(iter != data_.end() && iter->begin()->time_stamp <= expiration_time) {
        size_ -= iter->size();
        iter = data_.erase(iter); // Erases the element and returns an iterator to the next element which is now the first element of the list
    }

}

//---------------------------------------------------------------------------------------------
template<typename _DataType,typename _TransformDataType> 
template< typename _Transform>
void Cluster<_DataType,_TransformDataType>::TransformMeasurements(const _Transform& transform) {
    for (auto outer_iter = data_.begin(); outer_iter != data_.end(); ++outer_iter) {
        for(auto inner_iter = outer_iter->begin(); inner_iter != outer_iter->end(); ++ inner_iter) {
            transform.TransformMeasurement(*inner_iter);
        }
    }
}

//---------------------------------------------------------------------------------------------
template<typename _DataType,typename _TransformDataType> 
template<typename _SourceContainer>
bool Cluster<_DataType,_TransformDataType>::IsNeighboringMeasurement(const _SourceContainer& source_container, const unsigned int source_index, const Parameters& params, const Measurement& meas) const {

    for (auto outer_iter = std::prev(data_.end()); outer_iter != data_.end(); --outer_iter) {

        // New measurement is too far away from any recent measurement
        if( source_container.GetTemporalDistance(meas, outer_iter->front(), params) > params.cluster_time_threshold_) {
            return false;
        }

        for(auto inner_iter = outer_iter->begin(); inner_iter != outer_iter->end(); ++ inner_iter){

            // If same time stamp, use the position distance
            if(inner_iter->time_stamp == meas.time_stamp) {
                if(source_container.GetSpatialDistance(*inner_iter, meas, params) <= params.cluster_position_threshold_) {
                    return true;
                }
            } else { // Else use the velocity distance
                if(source_container.GetVelocityDistance(*inner_iter,meas,params) <= params.cluster_velocity_threshold_) {
                    return true;
                }
            }

        }

    }

    return false;
}

} // namespace rransac
#endif // RRANSAC_DATA_CONTAINERS_CLUSTER_BASE_H_