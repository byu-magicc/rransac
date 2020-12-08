#ifndef RRANSAC_DATA_CONTAINERS_DATA_TREE_LIST_H_
#define RRANSAC_DATA_CONTAINERS_DATA_TREE_LIST_H_

#include "data_containers/data_tree/data_tree_base.h"
#include<list>

namespace rransac
{

/**
 * \class DataTreeList
 * Data tree that uses the STL std::list as the inner container.
 */ 

 
class DataTreeList : public DataTreeBase<std::list<Meas>, DataTreeList> {

public:

typedef DataTreeBase<std::list<Meas>, DataTreeList> ParentClass;

using ParentClass::IteratorPair;

/**
 * Default constructor
 */ 
DataTreeList() : ParentClass() {}


/**
 * Adds a measurement to the data tree. It is assumed that the measurement has a valid time stamp and valid data.
 * @param[in] meas The measurement to be added.
 */
void DerivedAddMeas(const Meas& meas);

/**
 * Removes a measurement from the data tree using iterators defined in iter_pair
 * @param[in] iter_pair Contains the matched iterator pair necessary to remove an element. * 
 */
void DerivedRemoveMeas(const typename ParentClass::IteratorPair& iter_pair) {
    iter_pair.outer_it->erase(iter_pair.inner_it);
    if (iter_pair.outer_it->size() == 0) {
        this->data_.erase(iter_pair.outer_it);
    }
}

/**
 * Attempts to find a cluster from the measurements in the data tree. If a cluster was not found, the size of the vector container of iterator pairs will be 
 * zero.
 * @param[in] source Any source that is defined in system. It is used to calculate different distances
 * @param[in] params The system parameters. It contains parameters that define the cluster.
 * @param[out] iter_pairs Contains iterators to the elements that make up the cluster. The outer container distinguishes the iterators by time
 * @return returns true if a cluster was found.
 */
template<typename tSource>
bool DerivedFindCluster(const tSource& source, const Parameters& params, std::list<std::list<typename ParentClass::IteratorPair>>& iter_pairs) const;

private:

template<typename tSource>
static bool IsNeighboringMeasurement(const tSource source, const Parameters& params, const Meas& meas, const std::list<std::list<typename ParentClass::IteratorPair>>& iter_pairs);


};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                            Definitions
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void DataTreeList::DerivedAddMeas(const Meas& meas) {


// There are no measurements. Just add it.
if (this->data_.begin() == this->data_.end()) {
    this->data_.emplace_back(std::list<Meas>{meas});
} 
// All measurements occurred before the new one. So add it to the back.
else if (this->data_.back().begin()->time_stamp < meas.time_stamp) {
    this->data_.emplace_back(std::list<Meas>{meas});
} 
// The new measurement occurred before all the other measurements
else if (this->data_.front().begin()->time_stamp > meas.time_stamp) {
    this->data_.emplace_front(std::list<Meas>{meas});
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
            this->data_.insert(iter.base(),std::list<Meas>{meas});
            break;
        } 
    }
}

}

//-------------------------------------------------------------------------------------------------
template<typename tSource>
bool DataTreeList::DerivedFindCluster(const tSource& source, const Parameters& params, std::list<std::list<typename ParentClass::IteratorPair>>& iter_pairs) const {

    bool cluster_found = false;
    typename ParentClass::IteratorPair iter_pair;

    // Create a list of iterators that point to the all of the measurements of the latest time step. They will be used to seed the find cluster method
    std::list<std::list::iterator> iter_lasts;
    for (auto iter = this->data_.back().begin(); iter!= this->data_.back().end(), ++iter)
        iter_lasts.push_back(iter);

    
    while (iter_lasts.size() > 0 ) {

        auto iter_last = iter_lasts.begin();
        iter_pair.outer_it = this->data_.back().begin();
        iter_pair.inner_it = *iter_last;
        iter_pairs.clear();
        iter_pairs.emplace_front(std::list<typename ParentClass::IteratorPair>{iter_pair});
        iter_lasts.erase(iter_last);

        for (auto p_iter = iter_pairs->begin()->begin(); p_iter != iter_paris->begin()->end(); ++p_iter) {
            for (auto l_iter = iter_lasts.begin(); l_iter != iter_lasts.end(); ++l_iter) {
                if( source.GetSpatialDistance(*p_iter, *l_iter, params) <= params.cluster_position_threshold_ ) {
                    iter_pair.inner_it = l_iter;
                    iter_pairs.begin()->push_back(*l_iter);
                    iter_lasts.erase(l_iter);
                }
            }
        }

    }


    
    auto data_outer_iter_previous = this->data_.end();
    auto data_outer_iter_current = std::prev(this->data_.end(),1);

    

    // neighbors same time

    // neighbors previous time

    // seed the cluster using a measurement from the latest time step
    for (auto data_inner_iter_last = this->data_.back().begin(); iter_last != this->data_.back().end(); ++iter_last) {
        iter_pairs.clear();
        iter_pair.outer_it = data_outer_iter_current;
        iter_pair.inner_it = data_inner_iter_last;
        iter_pairs.emplace_front(std::list<typename ParentClass::IteratorPair>{iter_pair});
        auto pair_outer_iter = iter_pairs.begin();

        // Add neighborhing measurements from the current time step
        for (auto data_inner_iter_last_other = std::next(data_inner_iter_last,1); data_inner_iter_last_other!= his->data_.back().end(); ++data_inner_iter_last_other) {

            // See if it is a neghboring measurements
            if (source.GetSpatialDistance(*data_inner_iter_last, *data_inner_iter_last_other, params) <= params.cluster_position_threshold_ ) {
                iter_pair.inner = data_inner_iter_last_other;
                iter_pairs.begin()->push_back(iter_pair);
            }
        }


        data_outer_iter_current = std::prev(data_outer_iter_current,1); // Go to the list element of the previous time stamp

        // Add neighborhing measurements to the measurements from different time steps
        while (data_outer_iter_current != this->data_.end()){

            // See if the measurement time stamp is within the specified distance. If it is, continue as normal; otherwise end
            if(source.GetTemporalDistance(*(data_outer_iter_current->begin()), *(pair_outer_iter->begin())) >= params.cluster_time_threshold_) {
                break; // The rest of the measurements are too far away temporally so end search
            } else {

                for (auto data_inner_iter = data_outer_iter_current.begin(); data_inner_iter != data_outer_iter_current; ++ data_inner_iter) {

                }

            }

        }

    }

}

//-------------------------------------------------------------------------------------------------------------------------------------------------

template<typename tSource>
bool DataTreeList::IsNeighboringMeasurement(const tSource source, const Parameters& params, const Meas& meas, const std::list<std::list<typename ParentClass::IteratorPair>>& iter_pairs) {
    

    for (auto outer_iter = iter_pairs.begin(); outer_iter!= iter_paris.end(); ++outer_iter) {

        // If the measurement is too far away temporally, return false
        if(source.GetTemporalDistance(meas, *(outer_iter->begin()), params) > params.cluster_time_threshold_)
            return false;

        // Search the measurements of the current outer iterator to see if it is a neighbor
        for (auto inner_iter = outer_iter->begin(); inner_iter != outer_iter->end(); ++inner_iter) {

            if(source.GetVelocityDistance(meas, *inner_iter,params) <= params.cluster_velocity_threshold_)
                return true;

        }
    }

    return false;

}



} // namespace rransac


#endif // RRANSAC_DATA_CONTAINERS_DATA_TREE_LIST_H_