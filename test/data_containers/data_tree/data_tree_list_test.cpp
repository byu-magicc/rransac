#include <gtest/gtest.h>

#include "data_containers/data_tree/data_tree_list.h"
#include "common/transformations/trans_homography.h"
#include "state.h"
#include <math.h>
#include "parameters.h"
#include <Eigen/Core>
#include <iostream>

using namespace rransac;
using namespace lie_groups;

class DataTreeListTest : public ::testing::Test {
public:


// DataTreeListTest(unsigned int num_meas) : num_meas_{num_meas_} {}

void SetUp() override {

Eigen::Matrix3d homography;
homography.setZero();
homography.block(2,2,1,1)<< 1;

trans_.Init();
trans_.SetData(homography);
Meas m;

Eigen::Matrix<double,1,1> time_;

for (unsigned int ii = 0; ii < num_meas_; ++ii) {

// Give it a random time stamp [-10,10] and pose [-20,20]
time_.setRandom();
m.time_stamp = std::round(time_(0,0)*10);
m.pose = Eigen::Matrix<double,2,1>::Random()*20;
m.source_index = source_id;
++source_id;
data_tree_.AddMeas(m);

}



}

typedef TransformHomography<R2_r2, Eigen::Matrix<double,4,4>> Transform;

Transform trans_;
DataTreeList<Transform> data_tree_;
Parameters params_;

unsigned int num_meas_ = 1000;
unsigned int source_id = 0;        // used to uniquely identify a measurement

};


// Test the add remove prune and transform measurements function

TEST_F(DataTreeListTest, AddRemovePruneTransformMeasurementsTest) {

unsigned int size=0;
unsigned int num_additional_meas = 100;


/////////////////////////////////
//    Transform Measurement Test
/////////////////////////////////

// After the transformation all of the measurements should be zero
data_tree_.TransformMeasurements(trans_);
for(auto outer_iter = data_tree_.data_.begin(); outer_iter != data_tree_.data_.end(); ++outer_iter) {
    for(auto inner_iter = outer_iter->begin(); inner_iter != outer_iter->end(); ++inner_iter) {
        ASSERT_DOUBLE_EQ(inner_iter->pose.norm(), 0); // Should be zero after the transformation
    }
}


////////////////////////////////////
//      Add measurement test
////////////////////////////////////


// Construct a vector of measurements to add.
Eigen::Matrix<double,1,1> time;
std::vector<Meas> vec_meas(num_additional_meas);
for (auto& m : vec_meas) {
    time.setRandom();
    m.time_stamp = std::round(time(0,0)*10);
    m.source_index = source_id;
    source_id++;
    m.pose = Eigen::Matrix<double,2,1>::Random();
}

data_tree_.AddMeasurements(vec_meas);


// Construct a list of measurements to add.
std::list<Meas> list_meas(num_additional_meas);
Meas list_m;
for (auto iter = list_meas.begin(); iter != list_meas.end(); ++iter) {
    time.setRandom();
    list_m.time_stamp = std::round(time(0,0)*10);
    list_m.pose = Eigen::Matrix<double,2,1>::Random();
    list_m.source_index = source_id;
    source_id++;
    (*iter) = list_m;
}

data_tree_.AddMeasurements(list_meas);


for (auto iter = data_tree_.data_.begin(); iter != data_tree_.data_.end(); ++ iter) {
    auto iter2 = std::next(iter,1);
    size += iter->size();
    double time_stamp = iter->front().time_stamp;


    // Make sure the previous element's time stamps come before the other ones
    if (iter2 != data_tree_.data_.end()) {
        ASSERT_LT( iter->front().time_stamp , iter2->front().time_stamp) << "Elements not sorted properly according to time";
    }

    for (auto iter3 = iter->begin(); iter3 != iter->end(); ++ iter3) {
        ASSERT_DOUBLE_EQ(iter3->time_stamp, time_stamp) << "Elements of the same list element don't have the same time stamp";
    }

}

ASSERT_EQ(data_tree_.Size(), size);

/////////////////////////////////
//    Remove Measurement Test
/////////////////////////////////

// Remove random measurements
typename DataTreeList<Transform>::IteratorPair iter_pair;
std::vector<DataTreeList<Transform>::IteratorPair> iter_pairs;
std::vector<Meas> meas_to_remove;
for(auto outer_iter = data_tree_.data_.begin(); outer_iter != data_tree_.data_.end(); ++outer_iter) {

    time.setRandom();
    auto inner_iter = outer_iter->begin();
    inner_iter = std::next(inner_iter, std::round(time(0,0)*10));

    // We dont want the last iterator since it contains nothing
    if (inner_iter == outer_iter->end())
        inner_iter++;

    iter_pair.outer_it = outer_iter;
    iter_pair.inner_it = inner_iter;
    iter_pairs.push_back(iter_pair);
    meas_to_remove.push_back(*inner_iter);

}

// Remove measurements
data_tree_.RemoveMeasurements(iter_pairs);

// Make sure measurements are removed
for (auto& m: meas_to_remove) {
    for(auto outer_iter = data_tree_.data_.begin(); outer_iter != data_tree_.data_.end(); ++outer_iter) {
        for(auto inner_iter = outer_iter->begin(); inner_iter != outer_iter->end(); ++inner_iter) {
            ASSERT_NE(inner_iter->source_index, m.source_index); // The source index uniquely identifies the measurement in this test. So we make sure we cannot
                                                                 // find it after removing it.
        }
    }
}


// Remove all of the elements from a single time step;
iter_pairs.clear();
meas_to_remove.clear();

if (data_tree_.data_.begin() != data_tree_.data_.end()) {

    // grab a random element
    auto outer_iter_sts = data_tree_.data_.begin();
    time.setRandom();
    outer_iter_sts = std::next(outer_iter_sts,std::round(time(0,0)*20));

    // We dont want the last one so we just iterate again
    if(outer_iter_sts == data_tree_.data_.end())
        outer_iter_sts++;

    for (auto inner_iter = outer_iter_sts->begin(); inner_iter != outer_iter_sts->end(); ++ inner_iter ) {
        iter_pair.outer_it = outer_iter_sts;
        iter_pair.inner_it = inner_iter;
        iter_pairs.push_back(iter_pair);
        meas_to_remove.push_back(*inner_iter);
    }

    data_tree_.RemoveMeasurements(iter_pairs);

    // Make sure measurements are removed
    for (auto& m: meas_to_remove) {
        for(auto outer_iter = data_tree_.data_.begin(); outer_iter != data_tree_.data_.end(); ++outer_iter) {
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
if (data_tree_.data_.begin() != data_tree_.data_.end()) {
double time_stamp_first = data_tree_.data_.begin()->begin()->time_stamp;
data_tree_.PruneDataTree(time_stamp_first);
ASSERT_NE(data_tree_.data_.begin()->begin()->time_stamp, time_stamp_first);

}

}

