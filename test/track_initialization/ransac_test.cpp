#include <gtest/gtest.h>
#include <stdlib.h>
#include <time.h>

#include "system.h"
#include "data_containers/cluster.h"
#include "track_initialization/ransac.h"

using namespace rransac;



TEST(RANSAC_TEST, GenerateMinimumSubsetTest) {

struct Empty {};

Cluster cluster;
Ransac<Empty> ransac;

srand(time(NULL));

Meas m;
m.pose = Eigen::Matrix<double,2,1>::Random();
m.type = MeasurementTypes::RN_POS;
unsigned int max_times = 20;
unsigned int max_meas_per_time = 100;

unsigned int num_times = 0;
while (num_times < 1) {
    num_times = rand() % max_times;
}

for (int ii = 0; ii < num_times; ++ii) {
    m.time_stamp = ii;
    unsigned int num_meas = 0;
    while (num_meas < 1) {
        num_meas = rand() % max_meas_per_time;
    }

    for (int jj = 0; jj < num_meas; ++jj) {
        m.pose << ii, jj;
        cluster.AddMeasurement(m);
    }
}



unsigned int num_min_subset = 0;

while (num_min_subset < 1) {
    num_min_subset = rand() % num_times;
}


std::vector<Cluster::ConstIteratorPair> meas_index = ransac.GenerateMinimumSubset(num_min_subset, cluster);

// Make sure there is a measurement from the current time step and that they are all from different time steps
int times[meas_index.size()];
bool meas_from_curr_time_step_found = false;
int index = 0;
for (auto iter_pair_iter = meas_index.begin(); iter_pair_iter != meas_index.end(); ++iter_pair_iter) {

    if(iter_pair_iter->inner_it->time_stamp == num_times -1) {
        meas_from_curr_time_step_found = true;
    }
    times[index] = iter_pair_iter->inner_it->time_stamp;
    // std::cout << "time " << index << ": " << iter_pair_iter->inner_it->time_stamp << std::endl;
    // std::cout << "meas " << index << ": " << iter_pair_iter->inner_it->pose(1,0) << std::endl;
    index++;
    
}

ASSERT_TRUE(meas_from_curr_time_step_found);

// Make sure that they are all from different time steps
for (int ii = 0; ii < meas_index.size(); ++ii) {
    for (int jj = ii+1; jj < meas_index.size(); ++jj) {
        ASSERT_FALSE(times[ii]==times[jj]);
    }
}



}