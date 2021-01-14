#include <gtest/gtest.h>
// #include "data_containers/data_tree/data_tree_cluster.h"
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
typedef TransformHomography<State> Transform;
typedef ModelRN<State, TransformHomography> Model;
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
params.RANSAC_minimum_subset_ = 5;

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
Meas m_;

};


// This tests the AddMeasurement, AddMeasurrement, MergeClusters and TransformMeasurements functions
TEST_F(DataTreeClustersTestObject, AddMeasurementsTest) {

double size = 0;
ASSERT_EQ(sys_.data_tree_.Size(), size);  // Make sure size is initialized to zero

// Add one measurement
m_.time_stamp = 0;
m_.pose << 0, 0;
sys_.data_tree_.AddMeasurement(sys_, m_);
++size;

ASSERT_EQ(sys_.data_tree_.Size(), size);
ASSERT_EQ(sys_.data_tree_.data_.front().data_.front().front().pose,m_.pose);

// Create three clusters of 10 measurements that are disjoint by a spacial distance of three;
Meas m1 = m_;
m1.pose << 0,0;

Meas m2 = m_;
m2.pose << 3,0;
sys_.data_tree_.AddMeasurement(sys_, m2);
++size;

Meas m3 = m_;
m3.pose << sqrt(9- pow(1.5,2)), sqrt(9- pow(1.5,2));
sys_.data_tree_.AddMeasurement(sys_, m3);
++size;

for (int ii = 1; ii < 10; ++ii) {

    m1.time_stamp = ii;
    m2.time_stamp = ii;
    m3.time_stamp = ii;
    sys_.data_tree_.AddMeasurement(sys_, m1);
    sys_.data_tree_.AddMeasurement(sys_, m2);
    sys_.data_tree_.AddMeasurement(sys_, m3);
    size += 3;
}


// There should be 30 measurements and three clusters with 10 measurements each
ASSERT_EQ(sys_.data_tree_.Size(), size);
ASSERT_EQ(sys_.data_tree_.data_.size(), 3);

for (auto iter1 = sys_.data_tree_.data_.begin()->data_.begin(); iter1 != sys_.data_tree_.data_.begin()->data_.end(); ++iter1 ) {
    ASSERT_EQ(iter1->front().pose, m1.pose);
}

for (auto iter2 = std::next(sys_.data_tree_.data_.begin())->data_.begin(); iter2 != std::next(sys_.data_tree_.data_.begin())->data_.end(); ++iter2) {
    ASSERT_EQ(iter2->front().pose, m2.pose);
}

for (auto iter2 = std::next(sys_.data_tree_.data_.begin(),2)->data_.begin(); iter2 != std::next(sys_.data_tree_.data_.begin(),2)->data_.end(); ++iter2) {
    ASSERT_EQ(iter2->front().pose, m3.pose);
}

// Add another measurement that will cause the three clusters to merge into one
m_.pose << 1.5, 1.5;
m_.time_stamp = 9;
sys_.data_tree_.AddMeasurement(sys_, m_);
++size;
ASSERT_EQ(sys_.data_tree_.Size(), size);
ASSERT_EQ(sys_.data_tree_.data_.size(), 1);
// Make sure that all of the measurements are there
for (auto outer_iter = sys_.data_tree_.data_.begin()->data_.begin(); outer_iter !=  std::prev(sys_.data_tree_.data_.begin()->data_.end()); ++outer_iter) {
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
sys_.data_tree_.AddMeasurements(sys_,vec_meas);

ASSERT_EQ(sys_.data_tree_.Size(), size);

// Construct a list of measurements to add.
std::list<Meas> list_meas(num_additional_meas);
Meas list_m;
for (auto iter = list_meas.begin(); iter != list_meas.end(); ++iter) {
    list_m.time_stamp = 10;
    list_m.pose = Eigen::Matrix<double,2,1>::Random();
    (*iter) = list_m;
}
size += num_additional_meas;
sys_.data_tree_.AddMeasurements(sys_,list_meas);

ASSERT_EQ(sys_.data_tree_.Size(), size);


// Add a bunch of random measurements
Eigen::Matrix<double,1,1> time;

for (int ii = 0; ii < 1000; ++ii) {
    
    time.setRandom();
    m_.time_stamp = time(0,0)*20;
    m_.pose = MatData::Random();
    ASSERT_NO_THROW(sys_.data_tree_.AddMeasurement(sys_, m_) );
}


// The transformation will send all of the measurement pose data to zero
sys_.data_tree_.TransformMeasurements(sys_.transformaion_);

for( auto cluster_iter = sys_.data_tree_.data_.begin(); cluster_iter !=sys_.data_tree_.data_.end(); ++cluster_iter ) {
    for(auto outer_iter = cluster_iter->data_.begin(); outer_iter != cluster_iter->data_.end(); ++outer_iter) {
        for(auto inner_iter = outer_iter->begin(); inner_iter != outer_iter->end(); ++inner_iter)
            ASSERT_DOUBLE_EQ(inner_iter->pose.norm(), 0);
    }
}

}

//-----------------------------------------------------------------------------------------------

TEST_F(DataTreeClustersTestObject, RemoveMeasurementsTest) {


Eigen::Matrix<double,1,1> time;
unsigned int num_measurements = 2000;
// Add a bunch of measurements and give them a unique source index to identify them
for (unsigned int ii; ii < num_measurements; ++ii) {
    time.setRandom();
    m_.time_stamp = time(0,0)*10;
    m_.source_index = ii;
    m_.pose = MatData::Random()*30;
    sys_.data_tree_.AddMeasurement(sys_,m_);
}

// Get random measurements to remove until all of them are removed
DataTreeClusters::MeasurementLocationInfo measurement_location;
std::vector<DataTreeClusters::MeasurementLocationInfo> measurements_location;

while (sys_.data_tree_.Size() != 0) {

    for(auto cluster_iter = sys_.data_tree_.data_.begin(); cluster_iter != sys_.data_tree_.data_.end(); ++cluster_iter) {

        for(auto outer_iter = cluster_iter->data_.begin(); outer_iter != cluster_iter->data_.end(); ++outer_iter ) {
            
            // Get an iterator to a random element
            auto inner_iter = outer_iter->begin();
            time.setRandom();
            inner_iter = std::next(inner_iter, std::round(fabs(time(0,0)*10)));  
            if(inner_iter == outer_iter->end())                      // Make sure it isn't the last iter
                ++inner_iter;
            measurement_location.cluster_iter = cluster_iter;
            measurement_location.iter_pair.outer_it = outer_iter;
            measurement_location.iter_pair.inner_it = inner_iter;
            measurements_location.push_back(measurement_location);
        }

    }

    // Make sure we are not removing zero elements
    ASSERT_GT(measurements_location.size(), 0);
    unsigned int size = sys_.data_tree_.Size();

    // Remove the elements
    sys_.data_tree_.RemoveMeasurements(measurements_location);
    size -= measurements_location.size();
    ASSERT_EQ(sys_.data_tree_.Size(), size);
}

}

//---------------------------------------------------------------------------------

TEST_F(DataTreeClustersTestObject, PruneMeasurementsTest) {

m_.time_stamp = 0;
m_.pose << 0,0;
sys_.data_tree_.AddMeasurement(sys_,m_);
m_.pose << 5,0;
sys_.data_tree_.AddMeasurement(sys_,m_);

m_.time_stamp = 1;
sys_.data_tree_.AddMeasurement(sys_,m_);
m_.pose << 0,0;
sys_.data_tree_.AddMeasurement(sys_,m_);

m_.time_stamp = 2;
sys_.data_tree_.AddMeasurement(sys_,m_);
m_.time_stamp = 3;
sys_.data_tree_.AddMeasurement(sys_,m_);

sys_.current_time_ = 3;

sys_.data_tree_.PruneDataTree(sys_, 0.1);

ASSERT_EQ(sys_.data_tree_.Size(), 3);
ASSERT_EQ(sys_.data_tree_.data_.size(),1);

// Remove everything
sys_.data_tree_.PruneDataTree(sys_, 10);
ASSERT_EQ(sys_.data_tree_.Size(), 0);
ASSERT_EQ(sys_.data_tree_.data_.size(),0);

// Add a bunch of measurements 
Eigen::Matrix<double,1,1> time;
unsigned int num_measurements = 2000;
for (unsigned int ii=0; ii < num_measurements; ++ii) {
    time.setRandom();
    m_.time_stamp = time(0,0)*10;
    m_.pose = MatData::Random()*10;
    sys_.data_tree_.AddMeasurement(sys_,m_);
}

// Get an expiration time from a measurement time stamp
double expiration_time;
auto iter = sys_.data_tree_.data_.begin();
time.setRandom()*iter->data_.size();
iter = std::next(iter,std::round(fabs(time(0,0)*sys_.data_tree_.data_.size())));  // Get a random cluster
if (iter == sys_.data_tree_.data_.end())
    --iter;

auto outer_iter = iter->data_.begin();
outer_iter = std::next(outer_iter,std::round(fabs(time(0,0)*iter->data_.size())));  // Get a random time step list
if (outer_iter == iter->data_.end())
    --outer_iter;

expiration_time = outer_iter->front().time_stamp;

// Add a measurement you know wont be pruned
m_.time_stamp = expiration_time + 1; 
sys_.current_time_ = m_.time_stamp; 
sys_.data_tree_.AddMeasurement(sys_,m_);

sys_.data_tree_.PruneDataTree(sys_, expiration_time);

// Make sure there are no measurements past the expiration time
for( auto cluster_iter = sys_.data_tree_.data_.begin(); cluster_iter !=sys_.data_tree_.data_.end(); ++cluster_iter ) {
    for(auto outer_iter = cluster_iter->data_.begin(); outer_iter != cluster_iter->data_.end(); ++outer_iter) {
        for(auto inner_iter = outer_iter->begin(); inner_iter != outer_iter->end(); ++inner_iter)
            ASSERT_GT(inner_iter->time_stamp, expiration_time);
    }
}

// The data tree should not be empty
ASSERT_NE(sys_.data_tree_.Size(), 0);

}



//-----------------------------------------------------------------------------------------------------------------------------------------------------------


TEST_F(DataTreeClustersTestObject, ConstructClusters) {

// Add a bunch of measurements 
unsigned int time_steps = sys_.params_.RANSAC_minimum_subset_*2;
for (unsigned int ii = 0; ii < time_steps; ++ii) {
    for (unsigned int jj = 0; jj < 200; ++jj) {
        m_.time_stamp = ii;

        // Make sure that one cluster is formed by controlling the pose at each time step
        if (jj == 0)
            m_.pose << 0,100;  
        else if (jj == 1) {
            m_.pose << -100,0;
        } else if (jj == 2) {
            m_.pose << 100,0;
        } else {
            m_.pose = MatData::Random()*30;
        }
        sys_.data_tree_.AddMeasurement(sys_,m_);
    }
    
}

sys_.data_tree_.ConstructClusters(sys_);

ASSERT_GE(sys_.clusters_.size(),3);

for (auto iter: sys_.clusters_) {
    ASSERT_GE(iter->data_.size(), sys_.params_.RANSAC_minimum_subset_);
}



}