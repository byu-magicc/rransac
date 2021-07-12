#ifndef RRANSAC_TRACK_INITIALIZATION_SEED_POLOCIES_SE2_POS_POLICY_H_
#define RRANSAC_TRACK_INITIALIZATION_SEED_POLOCIES_SE2_POS_POLICY_H_
#pragma once

#include <algorithm>
#include <math.h>
#include <functional>

#include "rransac/system.h"
#include "rransac/data_containers/cluster.h"
#include "rransac/common/models/model_SEN_pos_vel.h"
#include "rransac/common/sources/source_SE3_cam_depth.h"

#include "lie_groups/lie_groups/SO3.h"
#include "lie_groups/lie_algebras/se3.h"
#include "lie_groups/state.h"


namespace rransac
{

/**
 * \class SE3CamDepthSeedPolicy
 * This seeding policy is used with the model ModelSENPosVel, the source SourceSE3CamDepth, and the state SE3. It
 * assumes that the target has only translational velocity in the forward direction of the body frame, that the heading
 * is aligned with the velocity vector in the inertial frame, that the pitch of the target is between -pi/2 and pi/2, and
 * that the roll can be approximated as zero. 
 * 
 */ 

template<typename _Model>
class SE3CamDepthSeedPolicy {

public:

typedef _Model Model;
typedef typename Model::Base::SourceContainer SourceContainer;
typedef typename Model::State State;
typedef typename State::DataType DataType;
typedef typename Model::Base::TransformDataType TransformDataType;
typedef Cluster<DataType,TransformDataType> ClusterT;
typedef System<Model> Sys;
typedef Eigen::Matrix<DataType,3,1> VecPos; 
typedef Eigen::Matrix<DataType,3,3> MatRot;
typedef Eigen::Matrix<DataType,4,1> VecMeasPose;
typedef typename Model::Base::Measurement Measurement;
typedef typename Model::Transformation Transformation;

// Perform compatibility checks
static_assert(std::is_same< ModelSENPosVel<SourceContainer> , Model>::value, "SE3CamDepthSeedPolicy: The source is not compatible with the seed policy");
static_assert(std::is_same<State,lie_groups::SE3_se3::template StateTemplate<typename State::DataType>>::value, "SE3CamDepthSeedPolicy: The state is not compatible with the seed policy");


/**
 * This function sets the initial guess of the state estimate. This seed
 * @param[in] meas_subset A subset of measurements used to calculate the seed.
 * @param[in] sys The object that contains the R-RANASAC data.
 * @param[in] x The initial conditions of the optimization problem. The initial conditions will change according to the seed policy. 
 * @param[in] size The number parameters the optimization solver is optimizing over which is the size of the input x.
 */ 
    static void GenerateSeedPolicy(const std::vector<typename ClusterT::IteratorPair>& meas_subset, const Sys& sys, DataType x[Model::cov_dim_], const int size) {
        
        if (meas_subset.size() < 3)
            throw std::runtime_error("SE3CamDepthSeedPolicy::GenerateSeedPolicy: Minimum subset must be at least 3");


        // Sort the measurements in chronological order
        std::vector<typename ClusterT::IteratorPair> meas_subset_ordered = meas_subset;
        std::sort(meas_subset_ordered.begin(), meas_subset_ordered.end(), SortChronological);



        // Smooth the measurements with the same time stamp
        std::vector<VecPos> positions;
        std::vector<DataType> times;
        VecPos pos;
        pos.setZero();
        int num_added = 0;
        DataType current_time = meas_subset_ordered[0].inner_it->time_stamp;
        for (int ii = 0; ii < meas_subset_ordered.size(); ++ii){           


            if(meas_subset_ordered[ii].inner_it->time_stamp >  current_time) {
                positions.push_back(pos/num_added);
                times.push_back(current_time);
                // std::cout << "pos: " << std::endl << pos/num_added << std::endl;
                pos.setZero();
                num_added = 0;
                current_time = meas_subset_ordered[ii].inner_it->time_stamp;
            } 

            pos += ToCartesian( GetPoseTrackingFrame(*meas_subset_ordered[ii].inner_it));
            num_added++;            

        }
        positions.push_back(pos/num_added);
        times.push_back(current_time);
        // std::cout << "pos: " << std::endl << pos/num_added << std::endl;
        pos.setZero();
        num_added = 0;


        if (positions.size() < 3) {
            throw std::runtime_error("SE3CamDepthSeedPolicy::GenerateSeedPolicy: There must be measurements from at least three different time steps.");

        }

        // Get velocities and rotation matrices, and body frame velocities;
        std::vector<MatRot> rotations;
        std::vector<DataType> rho_xs;
        std::vector<DataType> time_intervals;

        for (int ii = 0; ii < positions.size()-1; ++ii) {
            DataType dt = times[ii+1] - times[ii];

            #ifdef DEBUG_BUILD
            if (dt< 1e-10){
                throw std::runtime_error("SE3CamDepthSeedPolicy::GenerateSeedPolicy: The time step between measurements is too small");
            }
            #endif

            time_intervals.push_back(dt);
            VecPos vel = (positions[ii+1] - positions[ii])/dt;
            rho_xs.push_back(vel.norm());
            MatRot R = ConstructRotation(vel);

            #ifdef DEBUG_BUILD
            // Make sure that the rotation matrix is unitary
            if ( (R*R.transpose() - MatRot::Identity()).norm() > 1e-8 ){
                throw std::runtime_error("SE3CamDepthSeedPolicy::GenerateSeedPolicy: Rotation Matrix poorly formed. ");
            }
            #endif

            rotations.push_back(R);

        }

        // Get angular velocities
        std::vector<VecPos> angular_velocities;
        for(int ii=0; ii < rotations.size()-1; ++ii) {
            VecPos omega = lie_groups::SO3<DataType>::OMinus(rotations[ii+1],rotations[ii]);
            angular_velocities.push_back(omega);
        }

        // get averages
        DataType rho_x = std::accumulate(rho_xs.begin(), rho_xs.end(), 0.0)/rho_xs.size();
        VecPos omega(VecPos::Zero());
        for (auto& w : angular_velocities){
            omega +=w;
        } 
        omega /= angular_velocities.size();

        // Construct the estimated position
        Eigen::Matrix<DataType,4,4> state;
        state.setIdentity();
        state.block(0,0,3,3) = rotations.back();
        state.block(0,3,3,1) = positions.back();
        Eigen::Matrix<DataType,6,1> vec_state = lie_groups::se3<DataType>::Log(state);
        // std::cout << "seed state: " << std::endl << state << std::endl;
        // std::cout << "seed vel: " << std::endl << rho_x << std::endl << omega << std::endl;

        x[0] = vec_state(0);
        x[1] = vec_state(1);
        x[2] = vec_state(2);
        x[3] = vec_state(3);
        x[4] = vec_state(4);
        x[5] = vec_state(5);
        x[6] = rho_x;
        x[7] = omega(0);
        x[8] = omega(1);
        x[9] = omega(2);

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
     * Converts the pose of the measurement to Cartesian coordinates. 
     */ 
    static VecPos ToCartesian(const VecMeasPose& meas) {
        // std::cout << "meas: " << std::endl << meas << std::endl;
        // std::cout << "pos: " << std::endl << meas(0,0)*meas.block(1,0,3,1) << std::endl;
        return meas(0,0)*meas.block(1,0,3,1);
        
    }

    /**
     * Extracts the pose of the measurement in the tracking frame.
     */
    static VecMeasPose GetPoseTrackingFrame(const Measurement& meas) {
        VecMeasPose pose;
        // std::cout << "meas: " << std::endl << meas.pose << std::endl;

        if(meas.transform_meas) {
            pose = Transformation::TransformMeasurement(meas,meas.transform_data_m_t).pose;
        } else {
            pose = meas.pose;
        }
        // std::cout << "pose: " << std::endl << pose << std::endl;

        return pose;
    } 

    /**
     * Uses the derivative of the position, i.e. velocity in the inertial frame,
     * to construct the rotation matrix.
     */ 
    static MatRot ConstructRotation(const VecPos& td) {
        MatRot R;
        R.block(0,0,3,1) = td.normalized();
        DataType C_th = R.block(0,0,2,1).norm(); // assume that -pi/2 < th < pi/2
        if (C_th > 1e-4) {
            R(0,1) = -R(1,0)/C_th;
            R(1,1) = R(0,0)/C_th;
        } else {
            R(0,1) = -R(1,0);
            R(1,1) = R(0,0);
        }
        R(2,1) = 0;
        VecPos v1 = R.block(0,0,3,1);
        VecPos v2 = R.block(0,1,3,1);
        R.block(0,2,3,1) = v1.cross(v2);
        return R;
    }



};



} // namespace rransac


#endif //RRANSAC_TRACK_INITIALIZATION_SEED_POLOCIES_SE2_POS_POLICY_H_