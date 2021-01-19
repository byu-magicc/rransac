#ifndef RRANSAC_DATA_CONTAINERS_CLUSTER_BASE_H_
#define RRANSAC_DATA_CONTAINERS_CLUSTER_BASE_H_

#include <list>

#include "common/measurement/measurement_base.h"

namespace rransac
{

template<typename tDataType = double>    
class Cluster {

public:

typedef tDataType DataType;

/**
 * \struct ConstIteratorPair
 * Contains a match of iterators that cannot change the data they point to. 
 * The outer iterator pertains to the outer container for the data object, and
 * the inner container pertains to the inner container for the data object. So the outer iterator points to an 
 * element in the list and the inner iterator points to and element of tContainer.
 */ 
struct ConstIteratorPair {
    typename std::list<std::list<Meas<DataType>>>::const_iterator outer_it;
    typename std::list<Meas<DataType>>::const_iterator inner_it;
};

/**
 * \struct IteratorPair
 * Contains a match of iterators. The outer iterator pertains to the outer container for the data object, and
 * the inner container pertains to the inner container for the data object. So the outer iterator points to an 
 * element in the list and the inner iterator points to and element of tContainer.
 */ 
struct IteratorPair {
    typename std::list<std::list<Meas<DataType>>>::iterator outer_it;
    typename std::list<Meas<DataType>>::iterator inner_it;
};

Cluster(){};
Cluster(const Meas<DataType>& meas) { AddMeasurement(meas); }

/**
 * Adds a measurement to the data tree. It is assumed that the measurement has a valid time stamp and valid data.
 * @param[in] meas The measurement to be added.
 */
void AddMeasurement(const Meas<DataType>& meas);

/**
 * Adds measurements from a container that has an iterator object. The measurements can have different time stamps.
 * It is assumed that the measurements have a valid time stamp and valid data. This function will call AddMeas for every
 * measurement.
 */ 
template<typename tContainerMeas>
void AddMeasurements(const tContainerMeas& measurements){    
    for (auto iter = measurements.begin(); iter != measurements.end(); ++iter) 
        AddMeasurement(*iter);
}


/**
 * Removes a measurement from the data tree using iterators defined in iter_pair
 * @param[in] iter_pair Contains the matched iterator pair necessary to remove an element. * 
 */
void RemoveMeasurement(const IteratorPair& iter_pair);


/**
 * Removes all of the measurements indicated by iter_pairs by calling RemoveMeas
 * @param[in] iter_pair Contains all of the matched iterator pairs that indicate which measurements to remove.
 */ 
template<typename tContainerIteratorPair>
void RemoveMeasurements(const tContainerIteratorPair& iter_pairs){
    for(auto iter = iter_pairs.begin(); iter != iter_pairs.end(); ++iter)
        RemoveMeasurement(*iter);
}

/**
 * Removes all of the measurements with a time stamp before or equal to
 * the expiration time. 
 */ 
void PruneCluster(const double expiration_time);

/**
 * Transforms the measurements using the transform provided. 
 */ 
template< typename tTransform>
void TransformMeasurements(const tTransform& transform);

/**
 * Returns the number of measurements in the cluster
 */
unsigned int Size() const {return size_;};

/**
 * Returns true if the measurement is a neighboring measurement to one of the more recent measurements.
 * A recent measurement is a measurements whose time stamp is within Parameters::cluster_time_threshold_ of the latest measurement
 */ 
template<typename tSource>
bool IsNeighboringMeasurement(const tSource& source, const Parameters& param, const Meas<DataType>& meas) const;


std::list<std::list<Meas<DataType>>> data_; /** Contains all of the measurements. The outer container separates the measurements according to time */


private:

unsigned int size_=0; /** The total number of measurements in the cluster */

};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                            Definitions
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename tDataType> 
void Cluster<tDataType>::AddMeasurement(const Meas<DataType>& meas) {

++size_;

// There are no measurements. Just add it.
if (this->data_.begin() == this->data_.end()) {
    this->data_.emplace_back(std::list<Meas<DataType>>{meas});
} 
// All measurements occurred before the new one. So add it to the back.
else if (this->data_.back().begin()->time_stamp < meas.time_stamp) {
    this->data_.emplace_back(std::list<Meas<DataType>>{meas});
} 
// The new measurement occurred before all the other measurements
else if (this->data_.front().begin()->time_stamp > meas.time_stamp) {
    this->data_.emplace_front(std::list<Meas<DataType>>{meas});
}
else {
    // Search from the end to the beginning to find where to place it
    for (auto iter = this->data_.rbegin(); iter != this->data_.rend(); ++iter) {
        
        

        if ((*iter).begin()->time_stamp == meas.time_stamp) {
            (*iter).push_back(meas);                                  // This will need to be different for the R*tree
            break;
        } else if ((*iter).begin()->time_stamp < meas.time_stamp) {
            // std::cout << "add to middle" << std::endl;
            // std::vector<M> tmp{meas};
            this->data_.insert(iter.base(),std::list<Meas<DataType>>{meas});
            break;
        } 
    }
}

}

//--------------------------------------------------------------------------------------------
template<typename tDataType> 
void Cluster<tDataType>::RemoveMeasurement(const IteratorPair& iter_pair) {

    --size_;

    iter_pair.outer_it->erase(iter_pair.inner_it);
    if (iter_pair.outer_it->size() == 0) {
        this->data_.erase(iter_pair.outer_it);
    }
}

//---------------------------------------------------------------------------------------------
template<typename tDataType> 
void Cluster<tDataType>::PruneCluster(const double expiration_time) {

    auto iter = data_.begin();
    while(iter != data_.end() && iter->begin()->time_stamp <= expiration_time) {
        size_ -= iter->size();
        iter = data_.erase(iter); // Erases the element and returns an iterator to the next element which is now the first element of the list
    }

}

//---------------------------------------------------------------------------------------------
template<typename tDataType> 
template< typename tTransform>
void Cluster<tDataType>::TransformMeasurements(const tTransform& transform) {
    if (!transform.transform_null_) {
        for (auto outer_iter = data_.begin(); outer_iter != data_.end(); ++outer_iter) {
            for(auto inner_iter = outer_iter->begin(); inner_iter != outer_iter->end(); ++ inner_iter) {
                transform.TransformMeasurement(*inner_iter);
            }
        }
    }
}

//---------------------------------------------------------------------------------------------
template<typename tDataType> 
template<typename tSource>
bool Cluster<tDataType>::IsNeighboringMeasurement(const tSource& source, const Parameters& params, const Meas<tDataType>& meas) const {

    for (auto outer_iter = std::prev(data_.end()); outer_iter != data_.end(); --outer_iter) {

        // New measurement is too far away from any recent measurement
        if( source.GetTemporalDistance(meas, outer_iter->front(), params) > params.cluster_time_threshold_)
            return false;

        for(auto inner_iter = outer_iter->begin(); inner_iter != outer_iter->end(); ++ inner_iter){

            // If same time stamp, use the position distance
            if(inner_iter->time_stamp == meas.time_stamp) {
                if(source.GetSpatialDistance(*inner_iter, meas, params) <= params.cluster_position_threshold_)
                    return true;
            } else { // Else use the velocity distance
                if(source.GetVelocityDistance(*inner_iter,meas,params) <= params.cluster_velocity_threshold_)
                    return true;
            }

        }

    }

    return false;
}

} // namespace rransac
#endif // RRANSAC_DATA_CONTAINERS_CLUSTER_BASE_H_