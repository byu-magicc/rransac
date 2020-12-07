#include <gtest/gtest.h>

#include "data_containers/data_tree/data_tree_list.h"
#include "common/transformations/trans_homography.h"
#include "state.h"
#include <math.h>
#include "parameters.h"
#include <Eigen/Core>

using namespace rransac;
using namespace lie_groups;

class DataTreeListTest : public ::testing::Test {
public:


// DataTreeListTest(unsigned int num_meas) : num_meas_{num_meas_} {}

void SetUp() override {

Eigen::Matrix3d homography;
homography.setRandom();

trans_.Init();
trans_.SetData(homography);
Meas m;

Eigen::Matrix<double,1,1> time_;

for (unsigned int ii = 0; ii < num_meas_; ++ii) {

// Give it a random time stamp [-10,10] and pose [-20,20]
time_.setRandom();
m.time_stamp = std::round(time_(0,0)*10);
m.pose = Eigen::Matrix<double,2,1>::Random()*20;
data_tree_.AddMeas(m);

}



}

typedef TransformHomography<R2_r2, Eigen::Matrix<double,4,4>> Transform;

Transform trans_;
DataTreeList<Transform> data_tree_;
Parameters params_;

unsigned int num_meas_ = 1000;

};


// Test the add measurements function

TEST_F(DataTreeListTest, AddMeasurementsTest) {

unsigned int size=0;

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

}