#include <gtest/gtest.h>
#include "Eigen/Core"
#include <vector>
#include "common/measurement/measurement_base.h"
#include "data_structures/consensus_set.h"
#include <stdlib.h>
#include <time.h> 
#include <chrono>

namespace rransac {

template<class M>
void Print(const ConsensusSet<M>& cs) {

for (auto it = cs.consensus_set_.begin(); it != cs.consensus_set_.end(); ++it) {
    std::cout << (*it).front().time_stamp << std::endl;
}

}

TEST(CONSENSUS_TEST, ADD_MEASUREMENT) {

MeasBase<int,int> m1, m2,m3,m4,m5,m6;
m1.time_stamp = 0;
m2.time_stamp = 0.1;
m3.time_stamp = 0.5;
m4.time_stamp = -0.1;
m5.time_stamp = 0.2;
m5.source_id = 1;
m6.time_stamp = 0.2;
m6.source_id = 2;
ConsensusSet<MeasBase<int,int>> cs;

cs.AddMeasToConsensusSet(m1);
ASSERT_EQ(cs.consensus_set_.size(), 1);



// m2 should come after m1
cs.AddMeasToConsensusSet(m2);
ASSERT_EQ(cs.consensus_set_.size(), 2);
ASSERT_EQ(cs.consensus_set_.back().back().time_stamp, m2.time_stamp);
ASSERT_EQ(cs.consensus_set_.front().front().time_stamp, m1.time_stamp);



// m3 should come after m1 and m2
cs.AddMeasToConsensusSet(m3);
std::list<std::vector<MeasBase<int,int>>>::iterator iter = cs.consensus_set_.begin();
ASSERT_EQ(cs.consensus_set_.size(), 3);
ASSERT_EQ(cs.consensus_set_.front().front().time_stamp, m1.time_stamp);
ASSERT_EQ((*(++iter)).back().time_stamp, m2.time_stamp);
ASSERT_EQ(cs.consensus_set_.back().front().time_stamp, m3.time_stamp);



// add a measurement to the front
cs.AddMeasToConsensusSet(m4);
iter = cs.consensus_set_.begin();
ASSERT_EQ(cs.consensus_set_.size(), 4);
ASSERT_EQ((*iter).front().time_stamp, m4.time_stamp);
ASSERT_EQ((*(++iter)).front().time_stamp, m1.time_stamp);
ASSERT_EQ((*(++iter)).back().time_stamp, m2.time_stamp);
ASSERT_EQ(cs.consensus_set_.back().front().time_stamp, m3.time_stamp);


// add a measurement to the middle
cs.AddMeasToConsensusSet(m5);

iter = cs.consensus_set_.begin();
ASSERT_EQ(cs.consensus_set_.size(), 5);
ASSERT_EQ((*iter).front().time_stamp, m4.time_stamp);
ASSERT_EQ((*(++iter)).front().time_stamp, m1.time_stamp);
ASSERT_EQ((*(++iter)).back().time_stamp, m2.time_stamp);
ASSERT_EQ((*(++iter)).back().time_stamp, m5.time_stamp);
ASSERT_EQ(cs.consensus_set_.back().front().time_stamp, m3.time_stamp);
// Print<Meas>(cs);

// std::cerr << "here 2" << std::endl;


// add a measurement with same time stamp as m5
cs.AddMeasToConsensusSet(m6);
// Print<Meas>(cs);
iter = cs.consensus_set_.begin();
ASSERT_EQ(cs.consensus_set_.size(), 5);
ASSERT_EQ((*iter).front().time_stamp, m4.time_stamp);
ASSERT_EQ((*(++iter)).front().time_stamp, m1.time_stamp);
ASSERT_EQ((*(++iter)).back().time_stamp, m2.time_stamp);
ASSERT_EQ((*(++iter)).front().time_stamp, m5.time_stamp);
ASSERT_EQ((*iter).front().source_id, m5.source_id);
ASSERT_EQ((*iter).back().source_id, m6.source_id);
ASSERT_EQ(cs.consensus_set_.back().front().time_stamp, m3.time_stamp);
// std::cerr << "here 3" << std::endl;


}

//----------------------------------------------------------------------

TEST(CONSENSUS_TEST, ADD_MEASUREMENTS) {


int num_meas = 1000;
ConsensusSet<MeasBase<int,int>> cs;
srand (time(NULL));
std::vector<MeasBase<int,int>> measurements;

for (int i = 0; i < num_meas; i++) {
    MeasBase<int,int> m1;
    m1.time_stamp = rand() % 200 - 100;
    measurements.push_back(m1);
}

// auto start = std::chrono::high_resolution_clock::now();
cs.AddMeasurementsToConsensusSet(measurements);
// auto finish = std::chrono::high_resolution_clock::now();
// std::chrono::duration<double> elapsed = finish - start;
// std::cout << "Elapsed time: " << elapsed.count() << " s\n";
// Print<Meas>(cs);

for (auto iter = cs.consensus_set_.begin(); iter != cs.consensus_set_.end(); ++iter) {

    // Make sure the time stamps are in ascending order
    if (iter != cs.consensus_set_.begin()) {
        // iter->front().time_stamp;
        ASSERT_LE((std::prev(iter,1))->front().time_stamp, iter->front().time_stamp);
    }

    for (auto it = (*iter).begin(); it != (*iter).end(); ++it) {
        if (it != (*iter).begin()) {
            // it -> time_stamp;
        ASSERT_EQ((*(std::prev(it,1))).time_stamp, it->time_stamp);
        }
    }
    

}
}


//----------------------------------------------------------------------

TEST(CONSENSUS_TEST, ADD_MEASUREMENTS_SAME_TIME_STAMP) {


ConsensusSet<MeasBase<int,int>> cs;
MeasBase<int,int> m1,m2,m3,m4,m5;
std::vector<MeasBase<int,int>> mv1, mv2, mv3, mv4, mv5, mv6;
m1.time_stamp = 0.1;
m2.time_stamp = -0.1;
m3.time_stamp = 0.4;
m4.time_stamp = 0.2;
m5.time_stamp = 0.2;
mv1.push_back(m1);
mv2.push_back(m2);
mv3.push_back(m3);
mv4.push_back(m4);
mv5.push_back(m5);
mv6.push_back(m5);
mv6.push_back(m5);

cs.AddMeasurementsToConsensusSetSameTimeStamp(mv1);
cs.AddMeasurementsToConsensusSetSameTimeStamp(mv2);
cs.AddMeasurementsToConsensusSetSameTimeStamp(mv3);
cs.AddMeasurementsToConsensusSetSameTimeStamp(mv4);
cs.AddMeasurementsToConsensusSetSameTimeStamp(mv5);
cs.AddMeasurementsToConsensusSetSameTimeStamp(mv6);

auto iter = cs.consensus_set_.begin();
ASSERT_EQ(cs.consensus_set_.size(), 4);
ASSERT_EQ(iter->front().time_stamp, m2.time_stamp);
ASSERT_EQ((*(++iter)).front().time_stamp, m1.time_stamp);
ASSERT_EQ((*(++iter)).front().time_stamp, m4.time_stamp);
ASSERT_EQ(iter->size(), 4);

// Add a bunch of random measurements
int num_meas = 1000;

srand (time(NULL));
std::vector<MeasBase<int,int>> measurements;

for (int i = 0; i < num_meas; i++) {
    MeasBase<int,int> m;
    m.time_stamp = rand() % 200 - 100;
    measurements.clear();
    measurements.push_back(m);
    cs.AddMeasurementsToConsensusSetSameTimeStamp(measurements);
}

// auto start = std::chrono::high_resolution_clock::now();
// cs.AddMeasurementsToConsensusSet(measurements);
// auto finish = std::chrono::high_resolution_clock::now();
// std::chrono::duration<double> elapsed = finish - start;
// std::cout << "Elapsed time: " << elapsed.count() << " s\n";
// Print<Meas>(cs);

for (auto iter = cs.consensus_set_.begin(); iter != cs.consensus_set_.end(); ++iter) {

    // Make sure the time stamps are in ascending order
    if (iter != cs.consensus_set_.begin()) {
        // iter->front().time_stamp;
        ASSERT_LE((std::prev(iter,1))->front().time_stamp, iter->front().time_stamp);
    }

    for (auto it = (*iter).begin(); it != (*iter).end(); ++it) {
        if (it != (*iter).begin()) {
            // it -> time_stamp;
        ASSERT_EQ((*(std::prev(it,1))).time_stamp, it->time_stamp);
        }
    }
    

}
}

//----------------------------------------------------------------------

TEST(CONSENSUS_TEST, PRUNE_CONSENSUS_SET) {

MeasBase<int,int> m1, m2,m3,m4,m5,m6;

m1.time_stamp = 0;
m2.time_stamp = 0.1;
m3.time_stamp = 0.5;
m4.time_stamp = -0.1;
m5.time_stamp = 0.2;
m6.time_stamp = 0.2;
std::vector<MeasBase<int,int>> meas {m1,m2,m3,m4,m5,m6};

ConsensusSet<MeasBase<int,int>> cs;
cs.AddMeasurementsToConsensusSet(meas);

// Remove all measurements
cs.PruneConsensusSet(1);
ASSERT_EQ(cs.consensus_set_.size(), 0);

cs.AddMeasurementsToConsensusSet(meas);

// Dont prune any measurements
cs.PruneConsensusSet(-0.5);
ASSERT_EQ(cs.consensus_set_.size(), 5);
ASSERT_EQ(cs.consensus_set_.front().front().time_stamp, m4.time_stamp);


// Prune the first measurement
cs.PruneConsensusSet(-0.01);
ASSERT_EQ(cs.consensus_set_.size(), 4);
ASSERT_EQ(cs.consensus_set_.front().front().time_stamp, m1.time_stamp);

// Prune the next two measurements
cs.PruneConsensusSet(0.15);
ASSERT_EQ(cs.consensus_set_.size(), 2);
ASSERT_EQ(cs.consensus_set_.front().front().time_stamp, m5.time_stamp);

}

}