#ifndef RRANSAC_TRACK_INITIALIZATION_SEED_POLOCIES_SE2_POS_POLICY_H_
#define RRANSAC_TRACK_INITIALIZATION_SEED_POLOCIES_SE2_POS_POLICY_H_
#pragma once

#include <algorithm>
#include <math.h>

#include "rransac/system.h"
#include "rransac/data_containers/cluster.h"
#include "lie_groups/lie_algebras/so2.h"

namespace rransac
{

/**
 * \class SE2PosSeedPolicy
 * The seed policy is used to seed the nonlinear log maximum likelihood estimation problem defined in NonLinearLMLEPolicy when the 
 * state of the target is \f$ SE2\times \mathfrak{se}2 \f$ and the measurement space is \f$ \mathbb{R}^2 \f$. It calculates the
 * state estimate under the assumption that the heading of the target is aligned with the velocity vector. This assumption means
 * that the first component of the body velocity is positive and the second is zero. 
 */ 

template<typename tModel>
class SE2PosSeedPolicy {

public:

typedef typename tModel::State State;
typedef typename State::DataType DataType;

/**
 * This function sets the initial guess of the state estimate. TODO:: This algorithm is detailed in the paper _____________
 * @param[in] meas_subset A subset of measurements used to calculate the seed.
 * @param[in] sys The object that contains the R-RANASAC data.
 * @param[in] x The initial conditions of the optimization problem. The initial conditions will change according to the seed policy. 
 * @param[in] size The number parameters the optimization solver is optimizing over which is the size of the input x.
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

    /** 
     * Used to sort the measurements in the provided measurement subset in chronological order by comparing the time stamps of
     * two measurements.
     * @param[in] meas_iter1 An iterator to the first measurement.
     * @param[in] meas_iter2 An iterator to the seconde measurement.
     */ 
    static bool SortChronological(typename Cluster<DataType>::IteratorPair& meas_iter1, typename Cluster<DataType>::IteratorPair& meas_iter2) {
        return meas_iter1.inner_it->time_stamp < meas_iter2.inner_it->time_stamp;
    }

};



} // namespace rransac


#endif //RRANSAC_TRACK_INITIALIZATION_SEED_POLOCIES_SE2_POS_POLICY_H_