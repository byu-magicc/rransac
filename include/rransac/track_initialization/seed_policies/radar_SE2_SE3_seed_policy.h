#ifndef RRANSAC_TRACK_INITIALIZATION_SEED_POLICIES__RADAR_SE2_SE3_SEED_POLICY_H_
#define RRANSAC_TRACK_INITIALIZATION_SEED_POLICIES__RADAR_SE2_SE3_SEED_POLICY_H_
#pragma once

#include "rransac/system.h"
#include "rransac/data_containers/cluster.h"

namespace rransac
{
    
/**
 * \class RadarSE2SE3SeedPolicy
 * The seed policies are used to seed the LMLE algorithm that finds a current hypothetical state estimate using a minimum subset.
 * If the target evolves on SE2 the measurements are assumed to be in polar coordinates. If the target evolves on SE3, the measurements
 * are assumed to be in spherical coordinates. @see SourceRadarSE2SE3.
 */ 

template<typename _Model>
class RadarSE2SE3SeedPolicy {

public:

typedef _Model Model;
typedef typename Model::Measurement Measurement;
typedef typename _Model::State State;
typedef typename State::DataType DataType;
typedef typename _Model::Base::TransformDataType TransformDataType;
typedef Cluster<DataType,TransformDataType> ClusterT;
static constexpr unsigned int meas_pose_dim_ = Model::SourceContainer::Source0::meas_pose_dim_;
typedef Eigen::Matrix<DataType,meas_pose_dim_,1> VecPoseMeas;
typedef Eigen::Matrix<DataType,State::Group::dim_pos_,1> VecPos; 
typedef Eigen::Matrix<DataType,State::Group::dim_rot_,1> VecAngularVel;

typedef typename Model::Transformation Transformation;
typedef Eigen::Matrix<DataType,State::Group::dim_pos_, State::Group::dim_pos_> MatRot;

static constexpr bool polar_coordinates_ = Model::SourceContainer::Source0::polar_coordinates_;

static_assert(lie_groups::utilities::StateIsSEN_seN<State>::value, "RadarSE2SE3SeedPolicy: The state is not compatible with the model");
static_assert(State::Group::dim_ != 3 || State::Group::dim_ != 6, "RadarSE2SE3SeedPolicy: The dimensions of the group element of the state must be either 2 or 3.");
static_assert(Model::SourceContainer::Source0::measurement_type_ == MeasurementTypes::SE2_SE3_RADAR, "RadarSE2SE3SeedPolicy: The measurement type is not compatible with the source."    );



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
        throw std::runtime_error("RadarSE2SE3SeedPolicy: Minimum subset must be at least 3");

    std::vector<typename ClusterT::IteratorPair> meas_subset_ordered = meas_subset;
    

    // Sort the measurements in chronological order
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

        pos += ToCartesian( GetPose(*meas_subset_ordered[ii].inner_it));
        num_added++;            

    }
    positions.push_back(pos/num_added);
    times.push_back(current_time);

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
    std::vector<VecAngularVel> angular_velocities;
    for(int ii=0; ii < rotations.size()-1; ++ii) {
        VecAngularVel omega = State::Group::RotGroup::OMinus(rotations[ii+1],rotations[ii])/time_intervals[ii];
        angular_velocities.push_back(omega);
    }

    // get averages
    DataType rho_x = std::accumulate(rho_xs.begin(), rho_xs.end(), 0.0)/rho_xs.size();
    VecAngularVel omega(VecAngularVel::Zero());
    for (auto& w : angular_velocities){
        omega +=w;
    } 
    omega /= angular_velocities.size();

    // Construct the estimated position
    State state;
    state.g_.data_.setIdentity();
    state.g_.data_.block(0,0,State::Group::dim_pos_,State::Group::dim_pos_) = rotations.back()*State::Group::RotAlgebra::Exp(omega*time_intervals.back());
    // state.g_.data_.block(0,0,State::Group::dim_pos_,State::Group::dim_pos_) = rotations.back();
    state.g_.data_.block(0,State::Group::dim_pos_,State::Group::dim_pos_,1) = positions.back();
    state.u_.data_.setZero();
    state.u_.p_(0) = rho_x;
    state.u_.th_ = omega;
    typename State::Vec_SC vec_state = State::Log(state);
    // std::cout << "seed state: " << std::endl << state << std::endl;
    // std::cout << "seed vel: " << std::endl << rho_x << std::endl << omega << std::endl;

    size_t state_index = 0;

    for(size_t ii = 0; ii < State::Group::dim_; ++ii) {
        x[state_index] = vec_state(ii);
        state_index++;
    }

    x[state_index] = rho_x;
    state_index++;

    for(size_t jj=State::Group::dim_ + State::Group::dim_pos_; jj < vec_state.size(); ++jj) {
        x[state_index] = vec_state(jj);
        state_index++;
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


    /**
     * Uses the derivative of the position, i.e. velocity in the inertial frame,
     * to construct the rotation matrix.
     */ 
    static MatRot ConstructRotation(const VecPos& td) {

        MatRot R;
        R.block(0,0,State::Group::dim_pos_,1) = td.normalized();

        if(polar_coordinates_) {

            R.block(0,1,State::Group::dim_pos_,1);
            R(0,1) =  -R(1,0);
            R(1,1) = R(0,0);

        } else {
            DataType C_th = R.block(0,0,2,1).norm(); // assume that -pi/2 < th < pi/2
            if (C_th > 1e-4) {
                R(0,1) = -R(1,0)/C_th;
                R(1,1) = R(0,0)/C_th;
            } else {
                R(0,1) = -R(1,0);
                R(1,1) = R(0,0);
            }
            R(2,1) = 0;
            Eigen::Matrix<DataType,3,1> v1 = R.block(0,0,3,1);
            Eigen::Matrix<DataType,3,1> v2 = R.block(0,1,3,1);
            R.block(0,2,State::Group::dim_pos_,1) = v1.cross(v2);
        }

        return R;
    }

};


} // namespace rransac



#endif // RRANSAC_TRACK_INITIALIZATION_SEED_POLICIES__RADAR_SE2_SE3_SEED_POLICY_H_