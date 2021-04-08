#ifndef RRANSAC_SYSTEM_H_
#define RRANSAC_SYSTEM_H_
#pragma once

#include <vector>
#include <list>

#include "rransac/parameters.h"
#include "rransac/common/measurement/measurement_base.h"
#include "rransac/data_containers/cluster.h"
#include "rransac/data_containers/data_tree/data_tree_cluster.h"

// // Sources
// #include "common/sources/source_base.h"
// #include "common/sources/source_RN.h"
// #include "common/sources/source_SEN_pos_vel.h"
// #include "common/sources/source_SEN_pose_twist.h"

// // Models
// #include "common/models/model_base.h"
// #include "common/models/model_RN.h"
// #include "common/models/model_SEN_pos_vel.h"
// #include "common/models/model_SEN_pose_twist.h"

// // Transformations
// #include "common/transformations/transformation_base.h"
// #include "common/transformations/transformation_null.h"
// #include "common/transformations/trans_homography.h"

namespace rransac {


/**
 * \class System
 * This class holds all of the system information including parameters, sources, and models. 
 * 
 */ 
template <typename tModel> 
class System {

typedef tModel Model;
typedef typename tModel::State State; 
typedef typename State::DataType DataType;
typedef typename tModel::Source Source; 

public:

    Parameters params_;                                        /** < System parameters */
    std::vector<typename tModel::Source> sources_;             /** < Contains all of the instantiated sources. */
    std::list<Meas<DataType>> new_meas_;                       /** < Contains all of the new measurements. */
    typename tModel::Transformation transformaion_;            /** < The transformation for the measurements and tracks */
    std::list<tModel> models_;                                 /** < The models created by rransac */
    std::vector<tModel*> good_models_;                         /** < A list of pointers to the good models */
    long int model_label_ =0;                                  /** < The label incrementer for good models */
    double current_time_ =0;                                   /** < The current system time */
    bool time_set_ = false;                                    /** < Indicates that the current time has not be set. The current time is set to the time stamp of the latest measurement. */
    double dt_ = 0;                                            /** < The time elapsed since the previous sensor scan. This is set in RRANSAC::AddMeasurements. */
    DataTreeClusters<DataType> data_tree_;                               /** < Contains measurements that are not in a consensus set */
    std::vector<typename std::list<Cluster<DataType>>::iterator> clusters_;       /** < Iterators to clusters. RANSAC tries to form measurements from each clusters */
    long int cluster_label_ =0; /**< When a cluster is elevated to a good cluster, it will receive a unique label. A good cluster is a cluster that 
                                     meets the minimum subset requirement as specified by Parameters::RANSAC_minimum_subset_. This label is for visualization purposes only */
    long int source_index_counter_ =0;                         /**< Used to ensure that sources are added with sequential source indices. */
    long int accumulative_number_of_tracks_ =0;                /** < The total number of tracks created by R-RANSAC. Each time a track is initialized, this value is increased by 1. */
    
    System() {
        transformaion_.Init();
    }

private:

};



} // rransac

#endif // RRANSAC_SYSTEM_H_
