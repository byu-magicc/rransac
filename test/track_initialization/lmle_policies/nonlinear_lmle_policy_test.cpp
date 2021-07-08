#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <math.h> 
#include <stdlib.h>
#include <time.h>


#include "lie_groups/state.h"
#include "rransac/system.h"
#include "rransac/common/transformations/transformation_null.h"
#include "rransac/common/models/model_SEN_pose_twist.h"
#include "rransac/common/models/model_SEN_pos_vel.h"
#include "rransac/common/models/model_RN.h"
#include "rransac/common/sources/source_R2_R3_radar.h"
#include "rransac/common/sources/source_base.h"
#include "rransac/data_containers/cluster.h"
#include "rransac/track_initialization/lmle_policies/nonlinear_lmle_policy.h"
#include "rransac/track_initialization/seed_policies/null_seed_policy.h"
#include "rransac/track_initialization/seed_policies/SE2_pos_seed_policy.h"
#include "rransac/track_initialization/seed_policies/radar_R2_R3_seed_policy.h"
#include "rransac/common/transformations/trans_radar_R2_R3_with_SE2_SE3.h"

using namespace rransac;
using namespace lie_groups;



// Define the state
typedef State<Rn,double,2,2> StateR2_2;
typedef State<Rn,double,2,3> StateR2_3;
typedef State<Rn,double,3,2> StateR3_2;
typedef State<Rn,double,3,3> StateR3_3;

// Define the sources
typedef SourceRadarR2_R3<StateR2_2,MeasurementTypes::R2_R3_RADAR,TransRadarR2R3WithSE2SE3> SourceR2_2Radar;
typedef SourceRadarR2_R3<StateR2_2,MeasurementTypes::R2_R3_RADAR_DEPTH_DERIV,TransRadarR2R3WithSE2SE3> SourceR2_2Radar_Depth;

typedef SourceRadarR2_R3<StateR2_3,MeasurementTypes::R2_R3_RADAR,TransRadarR2R3WithSE2SE3> SourceR2_3Radar;
typedef SourceRadarR2_R3<StateR2_3,MeasurementTypes::R2_R3_RADAR_DEPTH_DERIV,TransRadarR2R3WithSE2SE3> SourceR2_3Radar_Depth;

typedef SourceRadarR2_R3<StateR3_2,MeasurementTypes::R2_R3_RADAR,TransRadarR2R3WithSE2SE3> SourceR3_2Radar;
typedef SourceRadarR2_R3<StateR3_2,MeasurementTypes::R2_R3_RADAR_DEPTH_DERIV,TransRadarR2R3WithSE2SE3> SourceR3_2Radar_Depth;

typedef SourceRadarR2_R3<StateR3_3,MeasurementTypes::R2_R3_RADAR,TransRadarR2R3WithSE2SE3> SourceR3_3Radar;
typedef SourceRadarR2_R3<StateR3_3,MeasurementTypes::R2_R3_RADAR_DEPTH_DERIV,TransRadarR2R3WithSE2SE3> SourceR3_3Radar_Depth;

typedef SourceSENPoseTwist<SE3_se3,MeasurementTypes::SEN_POSE,TransformNULL> SourceSE3Pose;
typedef SourceSENPoseTwist<SE3_se3,MeasurementTypes::SEN_POSE_TWIST,TransformNULL> SourceSE3PoseTwist;

typedef SourceSENPosVel<SE2_se2,MeasurementTypes::SEN_POS,TransformNULL> SourceSE2Pos;
typedef SourceSENPosVel<SE2_se2,MeasurementTypes::SEN_POS_VEL,TransformNULL> SourceSE2PosVel;

// Define the source containers

typedef SourceContainer<SourceR2_2Radar,SourceR2_2Radar_Depth> SourceContainerR2_2Radar;
typedef SourceContainer<SourceR2_3Radar,SourceR2_3Radar_Depth> SourceContainerR2_3Radar;
typedef SourceContainer<SourceR3_2Radar,SourceR3_2Radar_Depth> SourceContainerR3_2Radar;
typedef SourceContainer<SourceR3_3Radar,SourceR3_3Radar_Depth> SourceContainerR3_3Radar;
typedef SourceContainer<SourceSE3Pose,SourceSE3PoseTwist> SourceContainerSE3;
typedef SourceContainer<SourceSE2Pos,SourceSE2PosVel> SourceContainerSE2;


// Setup the different test types
//---------------------------------------------------------------------------------------------------------------------------

struct TestR2_2Radar {

    typedef SourceContainerR2_2Radar SC;
    typedef ModelRN<SC> Model;
    typedef typename Model::Base::State State;
    typedef NonLinearLMLEPolicy<Model,RadarR2R3SeedPolicy> LMLEPolicy;
    typedef typename Model::Base::TransformDataType TransformDataType;
    typedef typename Model::Base::Measurement Measurement;
    typedef typename Model::Base::Transformation Transformation;
    bool transform_state_ = true;
    double noise_ = 1e-3;
    static constexpr double error_tolerance_ = 0.1;


    TransformDataType transform_data_;

    TestR2_2Radar(){
        transform_data_.setIdentity();
        if(transform_state_) {
            srand((unsigned int) time(0));
            transform_data_ = Transformation::GetRandomTransform(9);
        }
        std::cout << "transform data: " << std::endl << transform_data_ << std::endl;
        
    }

    
    
};

//---------------------------------------------------------------------------------------------------------------------------

struct TestR2_3Radar {

    typedef SourceContainerR2_3Radar SC;
    typedef ModelRN<SC> Model;
    typedef typename Model::Base::State State;
    typedef NonLinearLMLEPolicy<Model,RadarR2R3SeedPolicy> LMLEPolicy;
    typedef typename Model::Base::TransformDataType TransformDataType;
    typedef typename Model::Base::Measurement Measurement;
    typedef typename Model::Base::Transformation Transformation;
    bool transform_state_ = true;
    double noise_ = 1e-3;
    static constexpr double error_tolerance_ = 0.1;


    TransformDataType transform_data_;

    TestR2_3Radar(){
        transform_data_.setIdentity();
        if(transform_state_) {
            srand((unsigned int) time(0));
            
            transform_data_ = Transformation::GetRandomTransform(10);
        }
        
    }
    
};

//---------------------------------------------------------------------------------------------------------------------------

struct TestR3_2Radar {

    typedef SourceContainerR3_2Radar SC;
    typedef ModelRN<SC> Model;
    typedef typename Model::Base::State State;
    typedef NonLinearLMLEPolicy<Model,RadarR2R3SeedPolicy> LMLEPolicy;
    typedef typename Model::Base::TransformDataType TransformDataType;
    typedef typename Model::Base::Measurement Measurement;
    typedef typename Model::Base::Transformation Transformation;
    bool transform_state_ = true;
    double noise_ = 1e-3;
    static constexpr double error_tolerance_ = 0.1;


    TransformDataType transform_data_;

    TestR3_2Radar(){
        transform_data_.setIdentity();
        if(transform_state_) {
            srand((unsigned int) time(0));

            transform_data_ = Transformation::GetRandomTransform(10);
        }
        
    }
    
};

//---------------------------------------------------------------------------------------------------------------------------

struct TestR3_3Radar {

    typedef SourceContainerR3_3Radar SC;
    typedef ModelRN<SC> Model;
    typedef typename Model::Base::State State;
    typedef NonLinearLMLEPolicy<Model,RadarR2R3SeedPolicy> LMLEPolicy;
    typedef typename Model::Base::TransformDataType TransformDataType;
    typedef typename Model::Base::Measurement Measurement;
    typedef typename Model::Base::Transformation Transformation;
    bool transform_state_ = true;
    double noise_ = 1e-3;
    static constexpr double error_tolerance_ = 0.1;


    TransformDataType transform_data_;

    TestR3_3Radar(){
        transform_data_.setIdentity();
        if(transform_state_) {
            srand((unsigned int) time(0));

            transform_data_ = Transformation::GetRandomTransform(10);
        }
        
    }
    
};

//---------------------------------------------------------------------------------------------------------------------------

struct TestSE3PoseTwist {

    typedef SourceContainerSE3 SC;
    typedef ModelSENPoseTwist<SC> Model;
    typedef typename Model::Base::State State;
    typedef NonLinearLMLEPolicy<Model,NULLSeedPolicy> LMLEPolicy;
    typedef typename Model::Base::TransformDataType TransformDataType;
    typedef typename Model::Base::Measurement Measurement;
    typedef typename Model::Base::Transformation Transformation;
    bool transform_state_ = false;
    double noise_ = 1e-3;
    static constexpr double error_tolerance_ = 0.1;


    TransformDataType transform_data_;

    TestSE3PoseTwist(){

        
    }
    
};

//---------------------------------------------------------------------------------------------------------------------------

struct TestSE2PosVel {

    typedef SourceContainerSE2 SC;
    typedef ModelSENPosVel<SC> Model;
    typedef typename Model::Base::State State;
    typedef NonLinearLMLEPolicy<Model,SE2PosSeedPolicy> LMLEPolicy;
    typedef typename Model::Base::TransformDataType TransformDataType;
    typedef typename Model::Base::Measurement Measurement;
    typedef typename Model::Base::Transformation Transformation;
    bool transform_state_ = false;
    double noise_ = 1e-3;
    static constexpr double error_tolerance_ = 0.1;


    TransformDataType transform_data_;

    TestSE2PosVel(){ 

        
    }
    
};

//----------------------------------------------------------------------------------------------------------------------------------

template<typename T> 
class NonLinearLMLEPolicyTest : public testing::Test {

public:
typedef typename T::Model Model;
typedef System<typename T::Model> Sys;
typedef typename Sys::ClusterT ClusterT;
typedef typename T::TransformDataType TransformDataType;
typedef typename T::Measurement Measurement;

SourceParameters source_params0_, source_params1_;
Measurement m0_,m1_;
Parameters params_;
Sys sys_;
Model track_;
double start_time_ = 0;
double end_time_ = 0.5;
double dt_ = 0.1;
T test_traits_;
std::vector<typename ClusterT::IteratorPair> meas_subset_;
std::list<std::list<Measurement>> measurements_;

// protected:

void SetUp() override {
    
    // Setup source parameters

    source_params0_.meas_cov_ = T::SC::Source0::MatMeasCov::Identity()*test_traits_.noise_;
    source_params0_.source_index_ = 0;
    source_params0_.type_ = T::SC::Source0::measurement_type_;
    source_params1_.meas_cov_ =  T::SC::Source1::MatMeasCov::Identity()*test_traits_.noise_;
    source_params1_.source_index_ = 1;
    source_params1_.type_ = T::SC::Source1::measurement_type_;

    // Setup the measurement
    m0_.source_index = 0;
    m0_.type = T::SC::Source0::measurement_type_;
    m1_.source_index = 1;
    m1_.type = T::SC::Source1::measurement_type_;

    // Setup the system parameters
    params_.process_noise_covariance_ = Model::Base::MatModelCov::Identity()*test_traits_.noise_;
    params_.nonlinear_innov_cov_id_ = true;
    sys_.params_ = params_;
    sys_.source_container_.AddSource(source_params0_);
    sys_.source_container_.AddSource(source_params1_);

    // Generate minimum subset
    track_.state_ = Model::GetRandomState(10.0);
    // std::cout << "track pose" << std::endl << track_.state_.g_.data_ << std::endl;
    // std::cout << "track twist" << std::endl << track_.state_.u_.data_ << std::endl;
    // track_.state_.u_.data_*=10;

    Measurement tmp0, tmp1;
    
    std::list<Measurement> meas_time;
    
    typename ClusterT::IteratorPair iter_pair;



    for (double ii = start_time_; ii < end_time_; ii += dt_) {

        if (ii != start_time_) {
            track_.PropagateModel(dt_);
        }
        meas_time.clear();

        tmp0 = sys_.source_container_.GenerateRandomMeasurement(m0_.source_index,T::SC::Source0::MatMeasCov::Identity()*sqrt(test_traits_.noise_*0),track_.state_,test_traits_.transform_state_,test_traits_.transform_data_);
        tmp1 = sys_.source_container_.GenerateRandomMeasurement(m1_.source_index,T::SC::Source1::MatMeasCov::Identity()*sqrt(test_traits_.noise_*0),track_.state_,test_traits_.transform_state_,test_traits_.transform_data_);
        m0_.pose = tmp0.pose;
        m0_.time_stamp = ii;
        m0_.transform_state = test_traits_.transform_state_;
        m0_.transform_data_t_m = test_traits_.transform_data_;
        m1_.pose = tmp1.pose;
        m1_.twist = tmp1.twist;
        m1_.time_stamp = ii;
        m1_.transform_state = test_traits_.transform_state_;
        m1_.transform_data_t_m = test_traits_.transform_data_;
        meas_time.push_back(m0_);
        meas_time.push_back(m1_);
        measurements_.push_back(meas_time);
        sys_.current_time_ = ii;

    }

    for (auto outer_iter = measurements_.begin(); outer_iter != measurements_.end(); ++outer_iter) {
        for (auto inner_iter = outer_iter->begin(); inner_iter != outer_iter->end(); ++ inner_iter) {
            iter_pair.inner_it = inner_iter;
            iter_pair.outer_it = outer_iter;
            meas_subset_.push_back(iter_pair);
        }
    }

    std::random_shuffle(meas_subset_.begin(), meas_subset_.end());


}



};



using MyTypes = ::testing::Types<TestR2_2Radar,TestR2_3Radar,TestR3_2Radar,TestR3_3Radar,TestSE3PoseTwist,TestSE2PosVel>;
// using MyTypes = ::testing::Types<TestR2_2Radar>;
TYPED_TEST_SUITE(NonLinearLMLEPolicyTest, MyTypes);

TYPED_TEST(NonLinearLMLEPolicyTest, BigTest){
    
    typename TypeParam::Model::State state;
    bool success = false;
    typename TypeParam::LMLEPolicy policy;
    state = policy.GenerateHypotheticalStateEstimatePolicy(this->meas_subset_,this->sys_,success);
    double error_tolerance = TypeParam::error_tolerance_;

    // std::cout << "pose error " << std::endl << this->track_.state_.g_.data_-state.g_.data_ << std::endl;
    // std::cout << "twist error" << std::endl << this->track_.state_.u_.data_-state.u_.data_ << std::endl;

    ASSERT_LT( (this->track_.state_.g_.data_-state.g_.data_).norm(), error_tolerance  );
    ASSERT_LT( (this->track_.state_.u_.data_-state.u_.data_).norm(), error_tolerance  );



}



