#include <gtest/gtest.h>
#include "data_containers/data_tree/data_tree_cluster.h"
#include "system.h"
#include "common/sources/source_RN.h"
#include "common/models/model_RN.h"
#include "common/transformations/trans_homography.h"

using namespace lie_groups;
using namespace rransac;

class DataTreeClustersTestObject : public ::testing::Test {

public:

typedef R2_r2 State;
typedef SourceR2 Source;
typedef TransformHomography<State, Eigen::Matrix4d> Transform;
typedef ModelRN<State, Transform> Model;
typedef Eigen::Matrix<double,2,1> MatData;

void SetUp() override {

// Setup system with source and parameters
Source source;
SourceParameters source_params;
source_params.expected_num_false_meas_ = 0.1;
source_params.gate_probability_ = 0.8;
source_params.probability_of_detection_ = 0.8;
source_params.meas_cov_fixed_ = false;
source_params.type_ = MeasurementTypes::RN_POS;
source.Init(source_params);

Parameters params;
params.cluster_time_threshold_ = 2;
params.cluster_position_threshold_ = 2.2;
params.cluster_velocity_threshold_ = 1;

Transform trans;
Eigen::Matrix3d homography;
homography.setZero();
homography.block(2,2,1,1)<< 1;
trans.Init();
trans.SetData(homography);

sys_.sources_.push_back(source);
sys_.params_ = params;
sys_.transformaion_ = trans;

// setup measurement
m_.type = MeasurementTypes::RN_POS;
m_.pose = MatData::Zero();

}

System<Model> sys_;
DataTreeClusters data_tree_;
Meas m_;

};


// This tests the AddMeasurement, AddMeasurrement, and MergeClusters functions
TEST_F(DataTreeClustersTestObject, AddMeasurementsTest) {

double size = 0;
ASSERT_EQ(data_tree_.Size(), size);  // Make sure size is initialized to zero

// Add one measurement
m_.time_stamp = 0;
m_.pose << 0, 0;
data_tree_.AddMeasurement(sys_, m_);
++size;

ASSERT_EQ(data_tree_.Size(), size);
ASSERT_EQ(data_tree_.data_.front().data_.front().front().pose,m_.pose);

// Create three clusters of 10 measurements that are disjoint by a spacial distance of three;
Meas m1 = m_;
m1.pose << 0,0;

Meas m2 = m_;
m2.pose << 3,0;
data_tree_.AddMeasurement(sys_, m2);
++size;

Meas m3 = m_;
m3.pose << sqrt(9- pow(1.5,2)), sqrt(9- pow(1.5,2));
data_tree_.AddMeasurement(sys_, m3);
++size;

for (int ii = 1; ii < 10; ++ii) {

    m1.time_stamp = ii;
    m2.time_stamp = ii;
    m3.time_stamp = ii;
    data_tree_.AddMeasurement(sys_, m1);
    data_tree_.AddMeasurement(sys_, m2);
    data_tree_.AddMeasurement(sys_, m3);
    size += 3;
}


// There should be 30 measurements and three clusters with 10 measurements each
ASSERT_EQ(data_tree_.Size(), size);
ASSERT_EQ(data_tree_.data_.size(), 3);

for (auto iter1 = data_tree_.data_.begin()->data_.begin(); iter1 != data_tree_.data_.begin()->data_.end(); ++iter1 ) {
    ASSERT_EQ(iter1->front().pose, m1.pose);
}

for (auto iter2 = std::next(data_tree_.data_.begin())->data_.begin(); iter2 != std::next(data_tree_.data_.begin())->data_.end(); ++iter2) {
    ASSERT_EQ(iter2->front().pose, m2.pose);
}

for (auto iter2 = std::next(data_tree_.data_.begin(),2)->data_.begin(); iter2 != std::next(data_tree_.data_.begin(),2)->data_.end(); ++iter2) {
    ASSERT_EQ(iter2->front().pose, m3.pose);
}

// Add another measurement that will cause the three clusters to merge into one
m_.pose << 1.5, 1.5;
m_.time_stamp = 9;
data_tree_.AddMeasurement(sys_, m_);
++size;
ASSERT_EQ(data_tree_.Size(), size);
ASSERT_EQ(data_tree_.data_.size(), 1);
// Make sure that all of the measurements are there
for (auto outer_iter = data_tree_.data_.begin()->data_.begin(); outer_iter !=  std::prev(data_tree_.data_.begin()->data_.end()); ++outer_iter) {
    auto inner_iter = outer_iter->begin();
    ASSERT_EQ(inner_iter->pose, m1.pose);
    ++inner_iter;
    ASSERT_EQ(inner_iter->pose, m2.pose);
    ++inner_iter;
    ASSERT_EQ(inner_iter->pose, m3.pose);

}

// Construct a vector of measurements to add.
unsigned int num_additional_meas = 100;
std::vector<Meas> vec_meas(num_additional_meas);
for (auto& m : vec_meas) {
    m.time_stamp = 10;
    m.pose = Eigen::Matrix<double,2,1>::Random();
}
size += num_additional_meas;
data_tree_.AddMeasurements(sys_,vec_meas);

ASSERT_EQ(data_tree_.Size(), size);

// Construct a list of measurements to add.
std::list<Meas> list_meas(num_additional_meas);
Meas list_m;
for (auto iter = list_meas.begin(); iter != list_meas.end(); ++iter) {
    list_m.time_stamp = 10;
    list_m.pose = Eigen::Matrix<double,2,1>::Random();
    (*iter) = list_m;
}
size += num_additional_meas;
data_tree_.AddMeasurements(sys_,list_meas);

ASSERT_EQ(data_tree_.Size(), size);


// Add a bunch of random measurements
Eigen::Matrix<double,1,1> time;

for (int ii = 0; ii < 1000; ++ii) {
    
    time.setRandom();
    m_.time_stamp = time(0,0)*20;
    m_.pose = MatData::Random();
    ASSERT_NO_THROW(data_tree_.AddMeasurement(sys_, m_) );
}


}