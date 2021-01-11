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
m.pose = Eigen::Matrix2d::Random();
m.type = MeasurementTypes::RN_POS;
unsigned int max_times = 20;
unsigned int max_meas_per_time = 100;

unsigned int num_times = 0;
while (num_times < 1) {
    num_times = rand() % max_times;
}

for (int ii = 0; ii < num_times; ++ii) {
    m.time_stamp = 0;
    unsigned int num_meas = 0;
    while (num_meas < 0) {
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


std::vector<Cluster::IteratorPair> meas_index = ransac.GenerateMinimumSubset(num_min_subset, cluster);


}