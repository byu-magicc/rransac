#ifndef RRANSAC_TRACK_INITIALIZATION_SEED_POLICIES__RADAR_R2_R3_SEED_POLICY_H_
#define RRANSAC_TRACK_INITIALIZATION_SEED_POLICIES__RADAR_R2_R3_SEED_POLICY_H_
#pragma once

#include "rransac/system.h"
#include "rransac/data_containers/cluster.h"

namespace rransac
{

/**
 * \class RadarR2R3SeedPolicy
 * The seed policies are used to seed the LMLE algorithm that finds a current hypothetical state estimate using a minimum subset.
 */ 

template<typename _Model>
class RadarR2R3SeedPolicy {

public:

typedef _Model Model;
typedef typename Model::Measurement Measurement;
typedef typename _Model::State State;
typedef typename State::DataType DataType;
typedef typename _Model::Base::TransformDataType TransformDataType;
typedef Cluster<DataType,TransformDataType> ClusterT;
static constexpr unsigned int meas_pose_dim_ = Model::SourceContainer::Source0::meas_pose_dim_;
typedef Eigen::Matrix<DataType,meas_pose_dim_,1> VecPoseMeas;
typedef typename Model::Transformation Transformation;


/**
 * This function sets the initial guess of the state estimate.
 * @param meas_subset The subset used to generate the seed for the LMLE algorithm
 * @param sys The object that contains the R-RANASAC data
 * @param x The initial state for the LMLE algorithm
 * @param size The size of the pointer x
 * 
 */ 
    static void GenerateSeedPolicy(const std::vector<typename ClusterT::IteratorPair>& meas_subset, const System<_Model>& sys, double x[_Model::cov_dim_], const int size) {
        
        if (meas_subset.size() < 3)
            throw std::runtime_error("SE2PosSeedPolicy: Minimum subset must be at least 3");

        std::vector<typename ClusterT::IteratorPair> meas_subset_ordered = meas_subset;
        
        VecPoseMeas pose1, pose2, pose3;
        VecPoseMeas pos1, pos2, pos3, vel, vel1, vel2;

        DataType t1,t2,t3;
        DataType dt1;
        DataType dt2;
     

        

        // Sort the measurements in chronological order
        std::sort(meas_subset_ordered.begin(), meas_subset_ordered.end(), SortChronological);
        unsigned int oldest_index = 0;
        unsigned int newest_index = meas_subset_ordered.size()-1;
        unsigned int middle_index = meas_subset_ordered.size()/2;
        if (middle_index <= oldest_index)
            middle_index = oldest_index +1;
        else if (middle_index >= newest_index)
            middle_index = newest_index -1;

        pose1 = GetPose(*meas_subset_ordered[newest_index].inner_it);
        pose2 = GetPose(*meas_subset_ordered[middle_index].inner_it);
        pose3 = GetPose(*meas_subset_ordered[oldest_index].inner_it);      

        t1 = meas_subset_ordered[newest_index].inner_it->time_stamp;
        t2 = meas_subset_ordered[middle_index].inner_it->time_stamp;
        t3 = meas_subset_ordered[oldest_index].inner_it->time_stamp;

        dt1 = t1 - t2;
        dt2 = t2 - t3;

        pos1 = ToCartesian(pose1);
        pos2 = ToCartesian(pose2);
        pos3 = ToCartesian(pose3);

        if(t1 != t2 && t2 != t3) {
            vel1 = (pos1 - pos2)/dt1;
            vel2 = (pos2 - pos3)/dt2;
            vel = 0.5*vel1 + 0.5*vel2;
        } else if (t1 != t2) {
            vel = (pos1 - pos2)/dt1;
        } else if (t2 != t3) {
            vel = vel2 = (pos2 - pos3)/dt2;
        } else {
            vel.setZero();
        }

       typename Model::VecCov state;
       state.setZero();
       state.block(0,0,meas_pose_dim_,1) = pos1;
       state.block(meas_pose_dim_,0,meas_pose_dim_,1) = vel;

       for (int ii = 0; ii < size; ++ii) {
           x[ii] = state(ii);
       }


    }


    private:

    /** 
     * Used to sort the measurements in the provided measurement subset in chronological order by comparing the time stamps of
     * two measurements.
     * @param[in] meas_iter1 An iterator to the first measurement.
     * @param[in] meas_iter2 An iterator to the seconde measurement.
     */ 
    static bool SortChronological(typename ClusterT::IteratorPair& meas_iter1, typename ClusterT::IteratorPair& meas_iter2) {
        return meas_iter1.inner_it->time_stamp < meas_iter2.inner_it->time_stamp;
    }


    /**
     * Gets the pose in the tracking frame.
     */ 
    static VecPoseMeas GetPose(const Measurement& meas) {

        VecPoseMeas pose;
        if(meas.transform_state) {
            pose = Transformation::TransformMeasurement(meas,meas.transform_data_t_m.inverse()).pose;
        } else {
            pose = meas.pose;
        }
        return pose;
    }

    /**
     * Converts the measurement pose from polar or spherical coordinates to cartesian.
     * 
     */ 
    static VecPoseMeas ToCartesian(const VecPoseMeas& pose) {

        VecPoseMeas position;
        DataType r = pose(0);
        DataType azimuth = pose(1);

        if(meas_pose_dim_ == 2) {
            position << r*cos(azimuth), r*sin(azimuth);

        } else if(meas_pose_dim_ == 3) {
            DataType zenith = pose(2);
            position << r*cos(azimuth)*sin(zenith),r*sin(azimuth)*sin(zenith),r*cos(zenith);

        }

        return position;
    }


};



} // namespace rransac


#endif //RRANSAC_TRACK_INITIALIZATION_SEED_POLOCIES_NULL_POLICY_H_