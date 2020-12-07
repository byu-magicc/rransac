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

template<typename tTransformation> 
class DataTreeList : public DataTreeBase<std::list<Meas>, tTransformation, DataTreeList<tTransformation>> {

public:

typedef DataTreeBase<std::list<Meas>, tTransformation, DataTreeList<tTransformation>> ParentClass;

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
 * @param[in] params The system parameters. It contains parameters that define the cluster.
 * @param[out] iter_pairs Contains iterators to the elements that make up the cluster.
 * @return returns true if a cluster was found.
 */
template<typename tContainerIteratorPairs>
bool DerivedFindCluster(const Parameters& params, tContainerIteratorPairs& iter_pairs) const;

private:

};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                            Definitions
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename tTransformation>
void DataTreeList<tTransformation>::DerivedAddMeas(const Meas& meas) {


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

// template<typename tTransformation>
// void DataTreeList<tTransformation>::DerivedRemoveMeas(const IteratorPair& iter_pair) {

    

// }


} // namespace rransac


#endif // RRANSAC_DATA_CONTAINERS_DATA_TREE_LIST_H_