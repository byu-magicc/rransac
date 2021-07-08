#include <gtest/gtest.h>
#include <Eigen/Core>
#include <vector>
#include <stdlib.h>
#include <time.h> 
#include <chrono>

#include "lie_groups/state.h"
#include "rransac/common/measurement/measurement_base.h"
#include "rransac/data_containers/consensus_set.h"
#include "rransac/common/transformations/transformation_base.h"


namespace rransac {

template<class M>
void Print(const ConsensusSet<M>& cs) {

    for (auto it = cs.consensus_set_.begin(); it != cs.consensus_set_.end(); ++it) {
        std::cout << (*it).front().time_stamp << std::endl;
    }
}

// --------------------------------------------------------------------------------------------
template<typename _State= lie_groups::R2_r2>
class TransformScalar : public TransformBase<TransformDerivedTraits<_State,double,Eigen::Matrix4d,false>, TransformScalar> {
public:

typedef TransformBase<TransformDerivedTraits<_State,double,Eigen::Matrix4d,false>, TransformScalar> Base;
typedef typename Base::State State;                                      /**< The State type being used. */
typedef typename Base::DataType DataType;                                /**< The scalar data type. */
typedef typename Base::TransformDataType TransformDataType;              /**< The transform data type being used. It is either an element of SE2 for R2 or SE3 for R3. */
typedef typename Base::MatCov MatCov;                                    /**< The covariance type of the track, and the transform jacobian type. */
typedef typename Base::Measurement Measurement;                          /**< The measurement type. */

void DerivedSetData(double data){ this->data_ = data;}

void DerivedTransformMeasurement(Measurement& meas) const{
    meas.pose = meas.pose*this->data_;
    meas.twist = meas.twist*this->data_;
}

void TransformTrack(lie_groups::R2_r2& state, Eigen::Matrix4d&cov) const {}

static TransformDataType DerivedGetRandomTransform(const DataType scalar = static_cast<DataType>(1.0)){
    return TransformDataType::Zero();
    }

};


// --------------------------------------------------------------------------------------------

TEST(CONSENSUS_TEST, ADD_MEASUREMENT) {


typedef Meas<double,Eigen::MatrixXd> Measurement;

Measurement m1, m2,m3,m4,m5,m6;
m1.time_stamp = 0;
m2.time_stamp = 0.1;
m3.time_stamp = 0.5;
m4.time_stamp = -0.1;
m5.time_stamp = 0.2;
m5.source_index = 1;
m6.time_stamp = 0.2;
m6.source_index = 2;
ConsensusSet<Measurement> cs;

cs.AddMeasToConsensusSet(m1);
ASSERT_EQ(cs.consensus_set_.size(), 1);



// m2 should come after m1
cs.AddMeasToConsensusSet(m2);
ASSERT_EQ(cs.consensus_set_.size(), 2);
ASSERT_EQ(cs.consensus_set_.back().back().time_stamp, m2.time_stamp);
ASSERT_EQ(cs.consensus_set_.front().front().time_stamp, m1.time_stamp);



// m3 should come after m1 and m2
cs.AddMeasToConsensusSet(m3);
std::list<std::vector<Measurement>>::iterator iter = cs.consensus_set_.begin();
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
ASSERT_EQ((*iter).front().source_index, m5.source_index);
ASSERT_EQ((*iter).back().source_index, m6.source_index);
ASSERT_EQ(cs.consensus_set_.back().front().time_stamp, m3.time_stamp);
// std::cerr << "here 3" << std::endl;


}

//----------------------------------------------------------------------

TEST(CONSENSUS_TEST, ADD_MEASUREMENTS) {

typedef Meas<double,Eigen::MatrixXd> Measurement;

int num_meas = 1000;
ConsensusSet<Measurement> cs;
srand (time(NULL));
std::vector<Measurement> measurements;

for (int i = 0; i < num_meas; i++) {
    Measurement m1;
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

typedef Meas<double,Eigen::MatrixXd> Measurement;

ConsensusSet<Measurement> cs;
Measurement m1,m2,m3,m4,m5;
std::vector<Measurement> mv1, mv2, mv3, mv4, mv5, mv6;
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
std::vector<Measurement> measurements;

for (int i = 0; i < num_meas; i++) {
    Measurement m;
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

typedef Meas<double,Eigen::MatrixXd> Measurement;
Measurement m1, m2,m3,m4,m5,m6;

m1.time_stamp = 0;
m2.time_stamp = 0.1;
m3.time_stamp = 0.5;
m4.time_stamp = -0.1;
m5.time_stamp = 0.2;
m6.time_stamp = 0.2;
std::vector<Measurement> meas {m1,m2,m3,m4,m5,m6};

ConsensusSet<Measurement> cs;
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

// -----------------------------------------------------------------

TEST(CONSENSUS_TEST, TRANSFORM_CONSENSUS_SET) {
typedef Meas<double,TransformScalar<>::TransformDataType> Measurement;
Measurement m1, m2,m3,m4,m5,m6;

m1.time_stamp = 0;
m2.time_stamp = 0.1;
m3.time_stamp = 0.5;
m4.time_stamp = -0.1;
m5.time_stamp = 0.2;
m6.time_stamp = 0.2;

m1.pose = Eigen::Vector2d::Random();
m2.pose = Eigen::Vector2d::Random();
m3.pose = Eigen::Vector2d::Random();
m4.pose = Eigen::Vector2d::Random();
m5.pose = Eigen::Vector2d::Random();
m6.pose = Eigen::Vector2d::Random();

m1.twist = Eigen::Vector2d::Random();
m2.twist = Eigen::Vector2d::Random();
m3.twist = Eigen::Vector2d::Random();
m4.twist = Eigen::Vector2d::Random();
m5.twist = Eigen::Vector2d::Random();
m6.twist = Eigen::Vector2d::Random();

std::vector<Measurement> meas {m1,m2,m3,m4,m5,m6};

ConsensusSet<Measurement> cs, cs_copy;
cs.AddMeasurementsToConsensusSet(meas);


TransformScalar<> trans;
double data = 2.0;
trans.SetData(data);

cs.TransformConsensusSet(trans);

auto iter = cs.consensus_set_.begin();

ASSERT_EQ( (*iter)[0].pose, m4.pose*data  );
++iter;
ASSERT_EQ( (*iter)[0].pose, m1.pose*data  );
++iter;
ASSERT_EQ( (*iter)[0].pose, m2.pose*data  );
++iter;
ASSERT_EQ( (*iter)[0].pose, m5.pose*data  );
ASSERT_EQ( (*iter)[1].pose, m6.pose*data  );
++iter;
ASSERT_EQ( (*iter)[0].pose, m3.pose*data  );


}

// -----------------------------------------------------------------

TEST(CONSENSUS_TEST, MERGE_CONSENSUS_SETS) {

typedef Meas<double,Eigen::MatrixXd> Measurement;
Measurement m1,m2,m3,m4,m5,m6,m7,m8,m9,m10;

m1.time_stamp = 0;
m2.time_stamp = 0.1;
m3.time_stamp = 0.5;
m4.time_stamp = -0.1;
m5.time_stamp = 0.2;
m6.time_stamp = 0.2;
m7.time_stamp = -0.2;
m8.time_stamp = 0.1;
m9.time_stamp = 0.3;
m10.time_stamp = 0.6;


std::vector<Measurement> meas1 {m1,m2,m3,m4,m5,m6};
std::vector<Measurement> meas2 {m7,m8,m9,m10};

ConsensusSet<Measurement> cs1;
ConsensusSet<Measurement> cs2;
ConsensusSet<Measurement> merged_cs;

cs1.AddMeasurementsToConsensusSet(meas1);
cs2.AddMeasurementsToConsensusSet(meas2);

merged_cs = ConsensusSet<Measurement>::MergeConsensusSets(cs1,cs2);

ASSERT_EQ(merged_cs.consensus_set_.size(), 8);

ASSERT_EQ(merged_cs.consensus_set_.front().front().time_stamp, m7.time_stamp);
ASSERT_EQ(merged_cs.consensus_set_.back().front().time_stamp, m10.time_stamp);

auto iter = merged_cs.consensus_set_.begin();
++iter;
ASSERT_EQ(iter->front().time_stamp, m4.time_stamp);
++iter;
ASSERT_EQ(iter->front().time_stamp, m1.time_stamp);
++iter;
ASSERT_EQ(iter->front().time_stamp, m2.time_stamp);
ASSERT_EQ(iter->size(), 2);
++iter;
ASSERT_EQ(iter->front().time_stamp, m5.time_stamp);
ASSERT_EQ(iter->size(), 2);
++iter;
ASSERT_EQ(iter->front().time_stamp, m9.time_stamp);
++iter;
ASSERT_EQ(iter->front().time_stamp, m3.time_stamp);


}

} // namespace rransac