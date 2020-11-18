#include <gtest/gtest.h>
#include "common/sources/source_base.h"
#include <stdlib.h>     /* srand, rand */
#include "memory.h"


namespace rransac {

template<class Derived, class S>
class Base
{
public:
    void Init(S& state){
        // std::static_cast<Derived*>(this)->Init(params_);
        static_cast<Derived*>(this)->Init(params_);
    }

    SourceParameters params_;
};

// class test : public

// Dummy function needed to initialize SourceBase
template <class S>
class Dummy : public SourceBase<S,Dummy<S>> 
// class Dummy : public Base<Dummy<S>,S> 
{
public:

    typedef S state_type_;

    /** Initializes the measurement source. This function must set the parameters.  */
    void Init(const SourceParameters& params) {
        this->params_ = params;
        this->H_ = Eigen::Matrix2d::Identity();
        this->V_ = Eigen::Matrix2d::Identity();
    }


    /** Returns the jacobian of the observation function w.r.t. the states */
    Eigen::MatrixXd GetLinObsMatState(const S& state){
        return this->H_;
    }                              

    /** Returns the jacobian of the observation function w.r.t. the sensor noise */
    Eigen::MatrixXd GetLinObsMatSensorNoise(const S& state){
        return this->V_;
    }                         

    /** Computes the estimated measurement given a state */
    Meas GetEstMeas(const S& state){
        Meas&& tmp = Meas();
        tmp.pose = state.g_.data_;
        return tmp;
    } /** Returns an estimated measurement according to the state. */


};

///////////////////////////////////////////////////////////////////////
//                             All functions execpt spatial test
///////////////////////////////////////////////////////////////////////

TEST(SOURCE_BASE, TemporalDistance) {




/* initialize random seed: */
srand (time(NULL));
SourceParameters source_params_;
Parameters params_;
source_params_.type_ = MeasurementTypes::RN_POS;
Dummy<lie_groups::R2_r2> source;

lie_groups::R2_r2 state;
state.g_.data_.setRandom();

Meas m1, m2;
m1.time_stamp = rand() % 100 -50;
m2.time_stamp = rand() % 100 -50;

// Test init function
source.Init(source_params_);
ASSERT_EQ(source.params_.type_, source_params_.type_);

// Test GetLinObsMatState
ASSERT_EQ(source.GetLinObsMatState(state), Eigen::Matrix2d::Identity());

// Test GetLinObsMatSensorNoise
ASSERT_EQ(source.GetLinObsMatSensorNoise(state), Eigen::Matrix2d::Identity());

// Test GetEstMeas
ASSERT_EQ(source.GetEstMeas(state).pose, state.g_.data_);

// lie_groups::R2_r2 state;
// SourceParameters params;
// source->Init(params);

ASSERT_EQ(source.GetTemporalDistance(m1,m2,params_),fabs(m1.time_stamp-m2.time_stamp));


// Random number generator;
int num_randn = 50000;
std::vector<Eigen::Matrix<double,4,1>> randn_nums(num_randn);
Eigen::Matrix<double,4,1> mean;
Eigen::Matrix<double,4,4> cov;
mean.setZero();
cov.setZero();

// Get random numbers and calculate mean
for(Eigen::Matrix<double,4,1>& randn_num : randn_nums) {
    randn_num = source.GaussianRandomGenerator(4);
    mean += randn_num;
    // std::cout << randn_num << std::endl << std::endl;
}

mean /= num_randn;

// Calculate std
for(Eigen::Matrix<double,4,1> randn_num : randn_nums) {
    cov += (mean - randn_num)*(mean-randn_num).transpose();
    
}

cov /= num_randn;

ASSERT_LE( (mean - Eigen::Matrix<double,4,1>::Zero()).norm(), 0.1);
ASSERT_LE( (cov - Eigen::Matrix<double,4,4>::Identity()).norm(), 0.1);


}

///////////////////////////////////////////////////////////////////////
//                             SPATIAL DISTANCE RN
///////////////////////////////////////////////////////////////////////

TEST(SOURCE_BASE, SpatialDistance_RN) {

/* initialize random seed: */
srand (time(NULL));
Parameters params_;

constexpr unsigned int max_iter = 20;
for (unsigned int i; i <max_iter; ++i) {

Dummy<lie_groups::R2_r2> source1;
Meas m_R2_Pos_1, m_R2_Pos_2, m_R2_Pos_Vel_1, m_R2_Pos_Vel_2;
m_R2_Pos_1.pose = Eigen::Matrix<double,max_iter,1>::Random().block(0,0,i,1);
m_R2_Pos_2.pose = Eigen::Matrix<double,max_iter,1>::Random().block(0,0,i,1);
m_R2_Pos_1.type = MeasurementTypes::RN_POS;
m_R2_Pos_2.type = MeasurementTypes::RN_POS;


m_R2_Pos_Vel_1.pose  = Eigen::Matrix<double,max_iter,1>::Random().block(0,0,i,1);
m_R2_Pos_Vel_1.twist = Eigen::Matrix<double,max_iter,1>::Random().block(0,0,i,1);
m_R2_Pos_Vel_2.pose  = Eigen::Matrix<double,max_iter,1>::Random().block(0,0,i,1);
m_R2_Pos_Vel_2.twist = Eigen::Matrix<double,max_iter,1>::Random().block(0,0,i,1);
m_R2_Pos_Vel_1.type  = MeasurementTypes::RN_POS_VEL;
m_R2_Pos_Vel_2.type  = MeasurementTypes::RN_POS_VEL;

ASSERT_DOUBLE_EQ(source1.GetSpatialDistance( m_R2_Pos_1,     m_R2_Pos_2,params_),         (m_R2_Pos_1.pose-    m_R2_Pos_2.pose).norm());
ASSERT_DOUBLE_EQ(source1.GetSpatialDistance( m_R2_Pos_1,     m_R2_Pos_Vel_1,params_),     (m_R2_Pos_1.pose-    m_R2_Pos_Vel_1.pose).norm());
ASSERT_DOUBLE_EQ(source1.GetSpatialDistance( m_R2_Pos_Vel_1, m_R2_Pos_Vel_2,params_),     (m_R2_Pos_Vel_1.pose-m_R2_Pos_Vel_2.pose).norm());
ASSERT_DOUBLE_EQ(source1.GetSpatialDistance( m_R2_Pos_Vel_1, m_R2_Pos_1,params_),         (m_R2_Pos_Vel_1.pose-m_R2_Pos_1.pose).norm());

}



}

///////////////////////////////////////////////////////////////////////
//             SPATIAL DISTANCE SEN POS with VEL
///////////////////////////////////////////////////////////////////////
TEST(SOURCE_BASE, SpatialDistance_SEN_POS) {

/* initialize random seed: */
srand (time(NULL));
Parameters params_;

constexpr unsigned int max_iter = 20;
for (unsigned int i; i <max_iter; ++i) {

Dummy<lie_groups::R2_r2> source1;
Meas m_R2_Pos_1, m_R2_Pos_2, m_R2_Pos_Vel_1, m_R2_Pos_Vel_2;
m_R2_Pos_1.pose = Eigen::Matrix<double,max_iter,1>::Random().block(0,0,i,1);
m_R2_Pos_2.pose = Eigen::Matrix<double,max_iter,1>::Random().block(0,0,i,1);
m_R2_Pos_1.type = MeasurementTypes::SEN_POS;
m_R2_Pos_2.type = MeasurementTypes::SEN_POS;


m_R2_Pos_Vel_1.pose  = Eigen::Matrix<double,max_iter,1>::Random().block(0,0,i,1);
m_R2_Pos_Vel_1.twist = Eigen::Matrix<double,max_iter,1>::Random().block(0,0,i,1);
m_R2_Pos_Vel_2.pose  = Eigen::Matrix<double,max_iter,1>::Random().block(0,0,i,1);
m_R2_Pos_Vel_2.twist = Eigen::Matrix<double,max_iter,1>::Random().block(0,0,i,1);
m_R2_Pos_Vel_1.type  = MeasurementTypes::SEN_POS_VEL;
m_R2_Pos_Vel_2.type  = MeasurementTypes::SEN_POS_VEL;

ASSERT_DOUBLE_EQ(source1.GetSpatialDistance( m_R2_Pos_1,     m_R2_Pos_2,params_),         (m_R2_Pos_1.pose-    m_R2_Pos_2.pose).norm());
ASSERT_DOUBLE_EQ(source1.GetSpatialDistance( m_R2_Pos_1,     m_R2_Pos_Vel_1,params_),     (m_R2_Pos_1.pose-    m_R2_Pos_Vel_1.pose).norm());
ASSERT_DOUBLE_EQ(source1.GetSpatialDistance( m_R2_Pos_Vel_1, m_R2_Pos_Vel_2,params_),     (m_R2_Pos_Vel_1.pose-m_R2_Pos_Vel_2.pose).norm());
ASSERT_DOUBLE_EQ(source1.GetSpatialDistance( m_R2_Pos_Vel_1, m_R2_Pos_1,params_),         (m_R2_Pos_Vel_1.pose-m_R2_Pos_1.pose).norm());

}



}

///////////////////////////////////////////////////////////////////////
//             SPATIAL DISTANCE SE2 and SE3 POSE with TWIST
///////////////////////////////////////////////////////////////////////

TEST(SOURCE_BASE, SpatialDistance_SEN_POSE) {

/* initialize random seed: */
srand (time(NULL));
Parameters params_;

// R2
Dummy<lie_groups::SE2_se2> source1;
Meas m_SE2_Pose_1, m_SE2_Pose_2, m_SE2_Pose_Twist_1, m_SE2_Pose_Twist_2;
m_SE2_Pose_1.pose  = lie_groups::se2::Exp(Eigen::Matrix<double,3,1>::Random());
m_SE2_Pose_2.pose  = lie_groups::se2::Exp(Eigen::Matrix<double,3,1>::Random());
m_SE2_Pose_1.type = MeasurementTypes::SEN_POSE;
m_SE2_Pose_2.type = MeasurementTypes::SEN_POSE;


m_SE2_Pose_Twist_1.pose  = lie_groups::se2::Exp(Eigen::Matrix<double,3,1>::Random());
m_SE2_Pose_Twist_1.twist = Eigen::Matrix<double,3,1>::Random();
m_SE2_Pose_Twist_2.pose  = lie_groups::se2::Exp(Eigen::Matrix<double,3,1>::Random());
m_SE2_Pose_Twist_2.twist = Eigen::Matrix<double,3,1>::Random();
m_SE2_Pose_Twist_1.type = MeasurementTypes::SEN_POSE_TWIST;
m_SE2_Pose_Twist_2.type = MeasurementTypes::SEN_POSE_TWIST;

ASSERT_DOUBLE_EQ(source1.GetSpatialDistance( m_SE2_Pose_1,       m_SE2_Pose_2,params_),           (lie_groups::SE2::OMinus(m_SE2_Pose_1.pose,      m_SE2_Pose_2.pose)).norm());
ASSERT_DOUBLE_EQ(source1.GetSpatialDistance( m_SE2_Pose_1,       m_SE2_Pose_Twist_1,params_),     (lie_groups::SE2::OMinus(m_SE2_Pose_1.pose,      m_SE2_Pose_Twist_1.pose)).norm());
ASSERT_DOUBLE_EQ(source1.GetSpatialDistance( m_SE2_Pose_Twist_1, m_SE2_Pose_Twist_2,params_),     (lie_groups::SE2::OMinus(m_SE2_Pose_Twist_1.pose,m_SE2_Pose_Twist_2.pose)).norm());
ASSERT_DOUBLE_EQ(source1.GetSpatialDistance( m_SE2_Pose_Twist_1, m_SE2_Pose_1,params_),           (lie_groups::SE2::OMinus(m_SE2_Pose_Twist_1.pose,m_SE2_Pose_1.pose)).norm());

// R3
Dummy<lie_groups::SE3_se3> source2;
Meas m_SE3_Pose_1, m_SE3_Pose_2, m_SE3_Pose_Twist_1, m_SE3_Pose_Twist_2;
m_SE3_Pose_1.pose  = lie_groups::se3::Exp(Eigen::Matrix<double,6,1>::Random());
m_SE3_Pose_2.pose  = lie_groups::se3::Exp(Eigen::Matrix<double,6,1>::Random());
m_SE3_Pose_1.type = MeasurementTypes::SEN_POSE;
m_SE3_Pose_2.type = MeasurementTypes::SEN_POSE;


m_SE3_Pose_Twist_1.pose  = lie_groups::se3::Exp(Eigen::Matrix<double,6,1>::Random());
m_SE3_Pose_Twist_1.twist = Eigen::Matrix<double,6,1>::Random();
m_SE3_Pose_Twist_2.pose  = lie_groups::se3::Exp(Eigen::Matrix<double,6,1>::Random());
m_SE3_Pose_Twist_2.twist = Eigen::Matrix<double,6,1>::Random();
m_SE3_Pose_Twist_1.type = MeasurementTypes::SEN_POSE_TWIST;
m_SE3_Pose_Twist_2.type = MeasurementTypes::SEN_POSE_TWIST;

ASSERT_DOUBLE_EQ(source2.GetSpatialDistance( m_SE3_Pose_1,       m_SE3_Pose_2,params_),           (lie_groups::SE3::OMinus(m_SE3_Pose_1.pose,      m_SE3_Pose_2.pose)).norm());
ASSERT_DOUBLE_EQ(source2.GetSpatialDistance( m_SE3_Pose_1,       m_SE3_Pose_Twist_1,params_),     (lie_groups::SE3::OMinus(m_SE3_Pose_1.pose,      m_SE3_Pose_Twist_1.pose)).norm());
ASSERT_DOUBLE_EQ(source2.GetSpatialDistance( m_SE3_Pose_Twist_1, m_SE3_Pose_Twist_2,params_),     (lie_groups::SE3::OMinus(m_SE3_Pose_Twist_1.pose,m_SE3_Pose_Twist_2.pose)).norm());
ASSERT_DOUBLE_EQ(source2.GetSpatialDistance( m_SE3_Pose_Twist_1, m_SE3_Pose_1,params_),           (lie_groups::SE3::OMinus(m_SE3_Pose_Twist_1.pose,m_SE3_Pose_1.pose)).norm());


}



} // namespace rransac