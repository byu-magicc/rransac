#include <gtest/gtest.h>
#include <math.h>
#include <Eigen/Core>
#include <iostream>

#include "lie_groups/state.h"
#include "rransac/data_containers/cluster.h"
#include "rransac/common/transformations/trans_homography.h"
#include "rransac/parameters.h"
#include "rransac/common/sources/source_RN.h"

using namespace rransac;
using namespace lie_groups;

class ClusterTestObject : public ::testing::Test {
public:


// DataTreeListTest(unsigned int num_meas) : num_meas_{num_meas_} {}

void SetUp() override {

Eigen::Matrix3d homography;
homography.setZero();
homography.block(2,2,1,1)<< 1;

trans_.Init();
trans_.SetData(homography);
Meas<double> m;

Eigen::Matrix<double,1,1> time_;

for (unsigned int ii = 0; ii < num_meas_; ++ii) {

// Give it a random time stamp [-10,10] and pose [-20,20]
time_.setRandom();
m.time_stamp = std::round(time_(0,0)*10);
m.pose = Eigen::Matrix<double,2,1>::Random()*20;
m.source_index = source_id;
++source_id;
cluster_.AddMeasurement(m);

}



}

typedef TransformHomography<R2_r2> Transform;

Transform trans_;
Cluster<double> cluster_;
Parameters params_;

unsigned int num_meas_ = 1000;
unsigned int source_id = 0;        // used to uniquely identify a measurement

};


// Test the add remove prune and transform measurements function

TEST_F(ClusterTestObject, AddRemovePruneTransformMeasurementsTest) {

unsigned int size=0;
unsigned int num_additional_meas = 100;


/////////////////////////////////
//    Transform Measurement Test
/////////////////////////////////

// After the transformation all of the measurements should be zero
cluster_.TransformMeasurements(trans_);
for(auto outer_iter = cluster_.data_.begin(); outer_iter != cluster_.data_.end(); ++outer_iter) {
    for(auto inner_iter = outer_iter->begin(); inner_iter != outer_iter->end(); ++inner_iter) {
        ASSERT_DOUBLE_EQ(inner_iter->pose.norm(), 0); // Should be zero after the transformation
    }
}


////////////////////////////////////
//      Add measurement test
////////////////////////////////////


// Construct a vector of measurements to add.
Eigen::Matrix<double,1,1> time;
std::vector<Meas<double>> vec_meas(num_additional_meas);
for (auto& m : vec_meas) {
    time.setRandom();
    m.time_stamp = std::round(time(0,0)*10);
    m.source_index = source_id;
    source_id++;
    m.pose = Eigen::Matrix<double,2,1>::Random();
}

cluster_.AddMeasurements(vec_meas);


// Construct a list of measurements to add.
std::list<Meas<double>> list_meas(num_additional_meas);
Meas<double> list_m;
for (auto iter = list_meas.begin(); iter != list_meas.end(); ++iter) {
    time.setRandom();
    list_m.time_stamp = std::round(time(0,0)*10);
    list_m.pose = Eigen::Matrix<double,2,1>::Random();
    list_m.source_index = source_id;
    source_id++;
    (*iter) = list_m;
}

cluster_.AddMeasurements(list_meas);


for (auto iter = cluster_.data_.begin(); iter != cluster_.data_.end(); ++ iter) {
    auto iter2 = std::next(iter,1);
    size += iter->size();
    double time_stamp = iter->front().time_stamp;


    // Make sure the previous element's time stamps come before the other ones
    if (iter2 != cluster_.data_.end()) {
        ASSERT_LT( iter->front().time_stamp , iter2->front().time_stamp) << "Elements not sorted properly according to time";
    }

    for (auto iter3 = iter->begin(); iter3 != iter->end(); ++ iter3) {
        ASSERT_DOUBLE_EQ(iter3->time_stamp, time_stamp) << "Elements of the same list element don't have the same time stamp";
    }

}

ASSERT_EQ(cluster_.Size(), size);

/////////////////////////////////
//    Remove Measurement Test
/////////////////////////////////

// Remove random measurements
typename Cluster<double>::IteratorPair iter_pair;
std::vector<Cluster<double>::IteratorPair> iter_pairs;
std::vector<Meas<double>> meas_to_remove;
for(auto outer_iter = cluster_.data_.begin(); outer_iter != cluster_.data_.end(); ++outer_iter) {

    time.setRandom();
    auto inner_iter = outer_iter->begin();
    inner_iter = std::next(inner_iter, std::round(fabs(time(0,0))*10));

    // We dont want the last iterator since it contains nothing
    if (inner_iter == outer_iter->end())
        inner_iter++;

    iter_pair.outer_it = outer_iter;
    iter_pair.inner_it = inner_iter;
    iter_pairs.push_back(iter_pair);
    meas_to_remove.push_back(*inner_iter);

}

// Remove measurements
cluster_.RemoveMeasurements(iter_pairs);
size -= iter_pairs.size();

// Make sure measurements are removed
for (auto& m: meas_to_remove) {
    for(auto outer_iter = cluster_.data_.begin(); outer_iter != cluster_.data_.end(); ++outer_iter) {
        for(auto inner_iter = outer_iter->begin(); inner_iter != outer_iter->end(); ++inner_iter) {
            ASSERT_NE(inner_iter->source_index, m.source_index); // The source index uniquely identifies the measurement in this test. So we make sure we cannot
                                                                 // find it after removing it.
        }
    }
}

ASSERT_EQ(cluster_.Size(), size);


// Remove all of the elements from a single time step;
iter_pairs.clear();
meas_to_remove.clear();

if (cluster_.data_.begin() != cluster_.data_.end()) {

    // grab a random element
    auto outer_iter_sts = cluster_.data_.begin();
    time.setRandom();
    outer_iter_sts = std::next(outer_iter_sts,std::round(time(0,0)*20));

    // We dont want the last one so we just iterate again
    if(outer_iter_sts == cluster_.data_.end())
        outer_iter_sts++;

    for (auto inner_iter = outer_iter_sts->begin(); inner_iter != outer_iter_sts->end(); ++ inner_iter ) {
        iter_pair.outer_it = outer_iter_sts;
        iter_pair.inner_it = inner_iter;
        iter_pairs.push_back(iter_pair);
        meas_to_remove.push_back(*inner_iter);
    }

    cluster_.RemoveMeasurements(iter_pairs);

    // Make sure measurements are removed
    for (auto& m: meas_to_remove) {
        for(auto outer_iter = cluster_.data_.begin(); outer_iter != cluster_.data_.end(); ++outer_iter) {
            for(auto inner_iter = outer_iter->begin(); inner_iter != outer_iter->end(); ++inner_iter) {
                ASSERT_NE(inner_iter->source_index, m.source_index); // The source index uniquely identifies the measurement in this test. So we make sure we cannot
                                                                    // find it after removing it.
            }
        }
    }

}



/////////////////////////////////
//    Prune Consensus Test
/////////////////////////////////

// Test prune consensus test
if (cluster_.data_.begin() != cluster_.data_.end()) {
double time_stamp_first = cluster_.data_.begin()->begin()->time_stamp;
cluster_.PruneCluster(time_stamp_first);
ASSERT_NE(cluster_.data_.begin()->begin()->time_stamp, time_stamp_first);

}

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
TEST(ClusterTest, IsNeighborTest) {

Parameters params;
params.cluster_velocity_threshold_ = 1;
params.cluster_position_threshold_ = 1;
params.cluster_time_threshold_ = 2;

Cluster<double> cluster;

SourceR2 source; // We need the source for calculating distances
SourceParameters source_params;
source_params.spacial_density_of_false_meas_ = 0.1;
source_params.gate_probability_ = 0.8;
source_params.probability_of_detection_ = 0.8;
source_params.meas_cov_ = Eigen::Matrix2d::Identity();
source_params.type_ = MeasurementTypes::RN_POS;
source.Init(source_params);

Meas<double> m, new_meas;
unsigned int num_meas = 10;
m.type = MeasurementTypes::RN_POS;
new_meas.type = MeasurementTypes::RN_POS;

////////////////////////////////////////////////////////
// The cluster has no measurements so their are no neighbor measurements
///////////////////////////////////////////////////////

new_meas.time_stamp = 0;
ASSERT_FALSE(cluster.IsNeighboringMeasurement(source,params,new_meas));


////////////////////////////////////////////////////////
// Measurements outside of time threshold wont be added
///////////////////////////////////////////////////////

for (int ii = 0; ii < num_meas; ++ii) {
    m.time_stamp = 0;
    m.pose = Eigen::Matrix<double,2,1>::Random()*10;
    cluster.AddMeasurement(m);
}

// Set the new measurement time stamp to be outside the time threshold
new_meas.time_stamp = m.time_stamp + params.cluster_time_threshold_ +0.1;
new_meas.pose = m.pose;

ASSERT_FALSE(cluster.IsNeighboringMeasurement(source,params,new_meas));

new_meas.time_stamp = m.time_stamp - params.cluster_time_threshold_ - 0.1;

ASSERT_FALSE(cluster.IsNeighboringMeasurement(source,params,new_meas));


////////////////////////////////////////////////////////
// measurement within time threshold
///////////////////////////////////////////////////////

new_meas.time_stamp =  m.time_stamp - params.cluster_time_threshold_;
ASSERT_TRUE(cluster.IsNeighboringMeasurement(source,params,new_meas));

// add more measurements
for (int ii = 0; ii < num_meas; ++ii) {
    m.time_stamp = 1;
    m.pose = Eigen::Matrix<double,2,1>::Random()*10;
    cluster.AddMeasurement(m);
}
new_meas.pose = m.pose;

// add more measurements
for (int ii = 0; ii < num_meas; ++ii) {
    m.time_stamp = 2;
    m.pose = Eigen::Matrix<double,2,1>::Random()*10;
    cluster.AddMeasurement(m);
}

// It should no longer be a neighbor due to recent time measurement
ASSERT_FALSE(cluster.IsNeighboringMeasurement(source,params,new_meas));

// Update time stamp so it is a neighbor
new_meas.time_stamp = 1;
ASSERT_TRUE(cluster.IsNeighboringMeasurement(source,params,new_meas));

new_meas.pose = m.pose;
ASSERT_TRUE(cluster.IsNeighboringMeasurement(source,params,new_meas));

// Add more measurements so that the new measurement is no longer a neighbor

// add more measurements
for (int ii = 0; ii < num_meas; ++ii) {
    m.time_stamp = 3.1;
    m.pose = Eigen::Matrix<double,2,1>::Random()*10;
    cluster.AddMeasurement(m);
}

ASSERT_FALSE(cluster.IsNeighboringMeasurement(source,params,new_meas));


}

