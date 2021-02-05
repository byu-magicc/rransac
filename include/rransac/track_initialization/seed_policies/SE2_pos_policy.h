#ifndef RRANSAC_TRACK_INITIALIZATION_SEED_POLOCIES_SE2_POS_POLICY_H_
#define RRANSAC_TRACK_INITIALIZATION_SEED_POLOCIES_SE2_POS_POLICY_H_
#pragma once

#include "system.h"
#include "data_containers/cluster.h"
#include <algorithm>
#include <math.h>
#include "lie_algebras/so2.h"

namespace rransac
{

/**
 * \class SE2PosSeedPolicy
 * The seed policies are used to seed the LMLE algorithm that finds a current hypothetical state estimate using a minimum subset.
 * This policy is used when the state is SE2 and only the position and velocity of the body w.r.t. the tracking frame expressed in the tracking frame
 * are observed.
 */ 

template<typename tModel>
class SE2PosSeedPolicy {

public:

typedef typename tModel::State State;
typedef typename State::DataType DataType;

/**
 * This function sets the initial guess of the state estimate to zero.
 * @param meas_subset The subset used to generate the seed for the LMLE algorithm. There must be at least three measurements from different time steps.
 * @param sys The object that contains the R-RANASAC data
 * @param x The initial state for the LMLE algorithm
 * @param size The size of the pointer x
 * 
 */ 
    static void GenerateSeedPolicy(const std::vector<typename Cluster<DataType>::IteratorPair>& meas_subset, const System<tModel>& sys, DataType x[tModel::cov_dim_], const int size) {
        
        if (meas_subset.size() < 3)
            throw std::runtime_error("SE2PosSeedPolicy: Minimum subset must be at least 3");

        std::vector<typename Cluster<DataType>::IteratorPair> meas_subset_ordered = meas_subset;
        Eigen::Matrix<DataType,2,1> td1, td2, rho, z; // Position derivatives
        Eigen::Matrix<DataType,3,1> se2;
        Eigen::Matrix<DataType,2,2> R1, R2;
        Eigen::Matrix<DataType,3,3> G = Eigen::Matrix<DataType,3,3>::Identity();
        bool velocity_set = false;

        // Sort the measurements in chronological order
        std::sort(meas_subset_ordered.begin(), meas_subset_ordered.end(), SortChronological);
        unsigned int oldest_index = 0;
        unsigned int newest_index = meas_subset_ordered.size()-1;
        unsigned int middle_index = meas_subset_ordered.size()/2;
        if (middle_index <= oldest_index)
            middle_index = oldest_index +1;
        else if (middle_index >= newest_index)
            middle_index = newest_index -1;
        


        if(meas_subset_ordered[newest_index].inner_it->type == MeasurementTypes::SEN_POS_VEL) {
            td1 = meas_subset_ordered[newest_index].inner_it->twist;
        } else {
            td1 = (meas_subset_ordered[newest_index].inner_it->pose - meas_subset_ordered[middle_index].inner_it->pose)/(meas_subset_ordered[newest_index].inner_it->time_stamp - meas_subset_ordered[middle_index].inner_it->time_stamp);
        }

        if(meas_subset_ordered[oldest_index].inner_it->type == MeasurementTypes::SEN_POS_VEL) {
            td2 = meas_subset_ordered[oldest_index].inner_it->twist;
        } else {
            td2 = (meas_subset_ordered[middle_index].inner_it->pose - meas_subset_ordered[oldest_index].inner_it->pose)/(meas_subset_ordered[middle_index].inner_it->time_stamp - meas_subset_ordered[oldest_index].inner_it->time_stamp);
        }



        DataType mag1, mag2;
        mag1 = td1.norm();
        mag2 = td2.norm();
        rho << mag1, 0;
        R1 << td1(0)/mag1, -td1(1)/mag1, td1(1)/mag1, td1(0)/mag1;
        R2 << td2(0)/mag2, -td2(1)/mag2, td2(1)/mag2, td2(0)/mag2;
        G.block(0,0,2,2) = R1;
        G.block(0,2,2,1) = meas_subset_ordered[newest_index].inner_it->pose;
        se2 = State::Algebra::Log(G);
        
        DataType dt = (meas_subset_ordered[oldest_index].inner_it->time_stamp - meas_subset_ordered[newest_index].inner_it->time_stamp );

        Eigen::Matrix<double,1,1> th_tmp = lie_groups::so2<DataType>::Log(R1.transpose()*R2);
        DataType thd = th_tmp(0,0)/dt;

        if (fabs(thd) < 1e-9)
            thd = 0;
        

        x[0] = se2(0);
        x[1] = se2(1);
        x[2] = se2(2);
        x[3] = (td1.norm()+td2.norm())/2.0;
        x[4] = thd;


    }

private:

    static bool SortChronological(typename Cluster<DataType>::IteratorPair& iter1, typename Cluster<DataType>::IteratorPair& iter2) {
        return iter1.inner_it->time_stamp < iter2.inner_it->time_stamp;
    }

    static double sgn(double val) {
    return (double(0) < val) - (val < double(0));
    }

};



} // namespace rransac


#endif //RRANSAC_TRACK_INITIALIZATION_SEED_POLOCIES_SE2_POS_POLICY_H_