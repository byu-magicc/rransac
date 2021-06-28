#include <gtest/gtest.h>
#include "lie_groups/state.h"
#include "rransac/common/transformations/transformation_null.h"
#include "rransac/common/sources/source_RN.h"
#include "rransac/common/sources/source_SEN_pose_twist.h"
#include "rransac/common/sources/source_container.h"
#include "rransac/parameters.h"

namespace rransac
{

using namespace lie_groups;


template <typename State>
struct CallbackClass0 {
    bool func(const State& state) {
        if (state.g_.data_.norm()  == 0)
            return true;
        else 
            return false;
}};

template <typename State>
struct CallbackClass1 {
    bool func(const State& state) {
        if (state.g_.data_.norm()  == 1)
            return true;
        else 
            return false;
}};

template <typename State>
struct CallbackClass2 {
    bool func(const State& state) {
        if (state.g_.data_.norm()  == 2)
            return true;
        else 
            return false;
}};

template <typename State>
struct CallbackClass3 {
    bool func(const State& state) {
        if (state.g_.data_.norm()  == 3)
            return true;
        else 
            return false;
}};

template <typename State>
struct CallbackClass4 {
    bool func(const State& state) {
        
            return false;
}};


typedef SourceRN<lie_groups::R2_r2,MeasurementTypes::RN_POS,TransformNULL> SourceR2Pos;
typedef SourceRN<lie_groups::R2_r2,MeasurementTypes::RN_POS_VEL,TransformNULL> SourceR2PosVel;

typedef SourceContainer<SourceR2Pos> SC1;
typedef SourceContainer<SourceR2Pos,SourceR2Pos,SourceR2Pos,SourceR2PosVel,SourceR2PosVel> SC;


class SourceContainerTest : public ::testing::Test
{

    protected:
    void SetUp() override {
    // Setup source parameters

    source_params.resize(num_sources);
    meas_types.resize(num_sources);
    meas_covs.resize(num_sources);
    indices.resize(num_sources);
    states.resize(num_sources);

    meas_types[0] = MeasurementTypes::RN_POS;
    meas_types[1] = MeasurementTypes::RN_POS;
    meas_types[2] = MeasurementTypes::RN_POS;
    meas_types[3] = MeasurementTypes::RN_POS_VEL;
    meas_types[4] = MeasurementTypes::RN_POS_VEL;

    meas_covs[0] = Eigen::Matrix2d::Identity();
    meas_covs[1] = Eigen::Matrix2d::Identity()*2;
    meas_covs[2] = Eigen::Matrix2d::Identity()*3;
    meas_covs[3] = Eigen::Matrix4d::Identity()*6;
    meas_covs[4] = Eigen::Matrix4d::Identity()*5;


    for (int ii = 0; ii < num_sources; ++ii) {
        source_params[ii].type_ = meas_types[ii];
        source_params[ii].meas_cov_ = meas_covs[ii];
        source_params[ii].source_index_ = ii;
        indices[ii] = ii;
        states[ii].g_.data_ << ii,0;
        states[ii].u_.data_ << ii,0;
    }

  
    
    }
private:
    /* data */
public:
    SourceContainerTest()=default;
    ~SourceContainerTest()=default;

    std::vector<SourceParameters> source_params; 
    std::vector<MeasurementTypes> meas_types;
    std::vector<Eigen::MatrixXd> meas_covs;
    std::vector<lie_groups::R2_r2> states;
 

    std::vector<int> indices;
    int num_sources = 5;


};





TEST_F(SourceContainerTest, ADD_SOURCE) {

    SC source_container_full;

    // Shouldn't have any issues adding these sources. Verify Parameters after adding them 
    for (int ii = 0; ii < source_params.size(); ++ii) {
        ASSERT_TRUE(source_container_full.AddSource(source_params[ii]));
        ASSERT_EQ(source_container_full.GetParams(ii).type_, meas_types[ii]);
        ASSERT_EQ(source_container_full.GetParams(ii).meas_cov_, meas_covs[ii]);
        ASSERT_EQ(source_container_full.GetParams(ii).source_index_, indices[ii]);
    }

    // Cant add sources again
    for (int ii = 0; ii < source_params.size(); ++ii) {
        ASSERT_ANY_THROW(source_container_full.AddSource(source_params[ii]));
    }

    // Adjust the parameters
    for (int ii =0; ii < num_sources; ++ii) {
        meas_covs[ii]*=2;
        source_params[ii].meas_cov_*=2;
        ASSERT_TRUE(source_container_full.ChangeSourceParameters(ii,source_params[ii]));
        ASSERT_EQ(source_container_full.GetParams(ii).meas_cov_, meas_covs[ii]);
    }

    

}

//------------------------------------------------------------------------------------------------------

TEST_F(SourceContainerTest, ADD_SOURCE_SURVEILLANCE_REGION) {

    SC source_container_full;
    bool transform_state = false;
    Eigen::MatrixXd EmptyMat;
    CallbackClass0<typename SourceR2Pos::State> call0;
    CallbackClass1<typename SourceR2Pos::State> call1;
    CallbackClass2<typename SourceR2Pos::State> call2;
    CallbackClass3<typename SourceR2PosVel::State> call3;
    CallbackClass4<typename SourceR2PosVel::State> call4;

    // Shouldn't have any issues adding these sources. Verify Parameters after adding them 
    ASSERT_TRUE(source_container_full.AddSource(source_params[0],std::bind(&CallbackClass0<typename SourceR2Pos::State>::func, call0, std::placeholders::_1)  ));
    ASSERT_TRUE(source_container_full.AddSource(source_params[1],std::bind(&CallbackClass1<typename SourceR2Pos::State>::func, call1, std::placeholders::_1)  ));
    ASSERT_TRUE(source_container_full.AddSource(source_params[2],std::bind(&CallbackClass2<typename SourceR2Pos::State>::func, call2, std::placeholders::_1)  ));
    ASSERT_TRUE(source_container_full.AddSource(source_params[3],std::bind(&CallbackClass3<typename SourceR2PosVel::State>::func, call3, std::placeholders::_1)  ));
    ASSERT_TRUE(source_container_full.AddSource(source_params[4],std::bind(&CallbackClass4<typename SourceR2PosVel::State>::func, call4, std::placeholders::_1)  ));

    for (int ii =0; ii < num_sources-1; ++ii) {
        ASSERT_TRUE(source_container_full.StateInsideSurveillanceRegion(ii,states[ii],transform_state,EmptyMat));
    }
    ASSERT_FALSE(source_container_full.StateInsideSurveillanceRegion(num_sources-1,states[num_sources-1],transform_state,EmptyMat));


    

}

//------------------------------------------------------------------------------------------------------

TEST_F(SourceContainerTest, JacobiansMeasOMinusRandomMeasurement) {

    SC source_container_full;
    bool transform_state = false;
    Eigen::MatrixXd EmptyMat;

    // Construct the Estimated Measurements
    typename SC::Source4::State state;
    state = SC::Source4::State::Random();
    typename SC::Source4 source;
    source.Init(source_params[4]);

    std::vector<Meas<double>> m1;
    std::vector<Meas<double>> m2;
    m1.resize(num_sources);
    m2.resize(num_sources);
    Eigen::Matrix<double,SC::Source4::meas_space_dim_*2,SC::Source4::meas_space_dim_*2> std;
    std.setIdentity();
    std *=0.1;

    for (int ii = 0; ii < num_sources; ++ii) {
        m1[ii] = source.GenerateRandomMeasurement(std,state,transform_state,EmptyMat);
        m2[ii] = source.GenerateRandomMeasurement(std,state,transform_state,EmptyMat);
    }
    
    


    for (int ii = 0; ii < source_params.size(); ++ii) {
        ASSERT_TRUE(source_container_full.AddSource(source_params[ii]));
        
    }

    ASSERT_EQ(source_container_full.GetLinObsMatState(0,states[0],transform_state,EmptyMat), SC::Source0::GetLinObsMatState(states[0],transform_state, EmptyMat));
    ASSERT_EQ(source_container_full.GetLinObsMatState(1,states[1],transform_state,EmptyMat), SC::Source1::GetLinObsMatState(states[1],transform_state, EmptyMat));
    ASSERT_EQ(source_container_full.GetLinObsMatState(2,states[2],transform_state,EmptyMat), SC::Source2::GetLinObsMatState(states[2],transform_state, EmptyMat));
    ASSERT_EQ(source_container_full.GetLinObsMatState(3,states[3],transform_state,EmptyMat), SC::Source3::GetLinObsMatState(states[3],transform_state, EmptyMat));
    ASSERT_EQ(source_container_full.GetLinObsMatState(4,states[4],transform_state,EmptyMat), SC::Source4::GetLinObsMatState(states[4],transform_state, EmptyMat));

    ASSERT_EQ(source_container_full.GetLinObsMatSensorNoise(0,states[0],transform_state,EmptyMat), SC::Source0::GetLinObsMatSensorNoise(states[0],transform_state, EmptyMat));
    ASSERT_EQ(source_container_full.GetLinObsMatSensorNoise(1,states[1],transform_state,EmptyMat), SC::Source1::GetLinObsMatSensorNoise(states[1],transform_state, EmptyMat));
    ASSERT_EQ(source_container_full.GetLinObsMatSensorNoise(2,states[2],transform_state,EmptyMat), SC::Source2::GetLinObsMatSensorNoise(states[2],transform_state, EmptyMat));
    ASSERT_EQ(source_container_full.GetLinObsMatSensorNoise(3,states[3],transform_state,EmptyMat), SC::Source3::GetLinObsMatSensorNoise(states[3],transform_state, EmptyMat));
    ASSERT_EQ(source_container_full.GetLinObsMatSensorNoise(4,states[4],transform_state,EmptyMat), SC::Source4::GetLinObsMatSensorNoise(states[4],transform_state, EmptyMat));
    
    ASSERT_EQ(source_container_full.GetEstMeas(0,states[0],transform_state,EmptyMat).pose, SC::Source0::GetEstMeas(states[0],transform_state, EmptyMat).pose);
    ASSERT_EQ(source_container_full.GetEstMeas(1,states[1],transform_state,EmptyMat).pose, SC::Source1::GetEstMeas(states[1],transform_state, EmptyMat).pose);
    ASSERT_EQ(source_container_full.GetEstMeas(2,states[2],transform_state,EmptyMat).pose, SC::Source2::GetEstMeas(states[2],transform_state, EmptyMat).pose);
    ASSERT_EQ(source_container_full.GetEstMeas(3,states[3],transform_state,EmptyMat).pose, SC::Source3::GetEstMeas(states[3],transform_state, EmptyMat).pose);
    ASSERT_EQ(source_container_full.GetEstMeas(4,states[4],transform_state,EmptyMat).pose, SC::Source4::GetEstMeas(states[4],transform_state, EmptyMat).pose);

    ASSERT_EQ(source_container_full.OMinus(0,m1[0],m2[0]), SC::Source0::OMinus(m1[0],m2[0]));
    ASSERT_EQ(source_container_full.OMinus(1,m1[1],m2[1]), SC::Source1::OMinus(m1[1],m2[1]));
    ASSERT_EQ(source_container_full.OMinus(2,m1[2],m2[2]), SC::Source2::OMinus(m1[2],m2[2]));
    ASSERT_EQ(source_container_full.OMinus(3,m1[3],m2[3]), SC::Source3::OMinus(m1[3],m2[3]));
    ASSERT_EQ(source_container_full.OMinus(4,m1[4],m2[4]), SC::Source4::OMinus(m1[4],m2[4]));
    

}


//------------------------------------------------------------------------------------------------------


TEST_F(SourceContainerTest, DistanceTests) {

    SC source_container_full;
    bool transform_state = false;
    Eigen::MatrixXd EmptyMat;

    // Construct the measurements
    Meas<double> m1, m2;
    typename SC::Source4::State state;
    state = SC::Source4::State::Random();
    typename SC::Source4 source;
    source.Init(source_params[4]);
    Eigen::Matrix<double,SC::Source4::meas_space_dim_*2,SC::Source4::meas_space_dim_*2> std;
    std.setIdentity();
    std *=0.1;
    m1 = source.GenerateRandomMeasurement(std,state,transform_state,EmptyMat);
    m2 = source.GenerateRandomMeasurement(std,state,transform_state,EmptyMat);
    m1.time_stamp = 0;
    m2.time_stamp = 1;
    m1.type = SC::Source4::meas_type_;
    m2.type = SC::Source4::meas_type_;
    m1.source_index = 4;
    m2.source_index = 4;

    // Setup parameters
    Parameters params;

    for (int ii = 0; ii < num_sources; ++ii) {

        ASSERT_EQ(source_container_full.GetTemporalDistance(m1,m2,params), source.GetTemporalDistance(m1,m2,params));
        ASSERT_EQ(source_container_full.GetSpatialDistance(m1,m2,params), source.GetSpatialDistance(m1,m2,params));
        ASSERT_EQ(source_container_full.GetVelocityDistance(m1,m2,params), source.GetVelocityDistance(m1,m2,params));

    }


}

    
} // namespace rransac
