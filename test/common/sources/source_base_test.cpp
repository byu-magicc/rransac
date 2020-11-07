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

    typedef S type_;

    /** Initializes the measurement source. This function must set the parameters.  */
    void Init(const SourceParameters& params) {
        this->params_ = params;
        // this->H_ = Eigen::Matrix2d::Identity();
        // this->V_ = Eigen::Matrix2d::Identity();
    }


    // /** Returns the jacobian of the observation function w.r.t. the states */
    // Eigen::MatrixXd GetLinObsMatState(const S& state){
    //     return this->H_;
    // }                              

    // /** Returns the jacobian of the observation function w.r.t. the sensor noise */
    // Eigen::MatrixXd GetLinObsMatSensorNoise(const S& state){
    //     return this->V_;
    // }                         

    // /** Computes the estimated measurement given a state */
    // Eigen::MatrixXd GetEstMeas(const S& state){
    //     return state.g_.data_;
    // } /** Returns an estimated measurement according to the state. */


};

///////////////////////////////////////////////////////////////////////
//                             Distance_Test
///////////////////////////////////////////////////////////////////////
typedef Dummy<lie_groups::R2_r2> DummyType;
typedef SourceBase<lie_groups::R2_r2,DummyType> SourceType;

TEST(SOURCE_BASE, TemporalDistance) {



/* initialize random seed: */
srand (time(NULL));
Parameters params_;
// SourceType* source = new Dummy<DummyType>;
// Base<DummyType,lie_groups::R2_r2>* source = new DummyType;
SourceBase<lie_groups::R2_r2,DummyType>* source = new Dummy<lie_groups::R2_r2>;
// std::unique_ptr<SourceType> source { new Dummy<DummyType>() };
Meas m1, m2;
m1.time_stamp = rand() % 100 -50;
m2.time_stamp = rand() % 100 -50;

// lie_groups::R2_r2 state;
// SourceParameters params;
// source->Init(params);

ASSERT_EQ(source->GetTemporalDistance(m1,m2,params_),fabs(m1.time_stamp-m2.time_stamp));

}

///////////////////////////////////////////////////////////////////////
//                             SPATIAL DISTANCE R2 and R3
///////////////////////////////////////////////////////////////////////

// TEST(SOURCE_BASE, SpatialDistance_RN) {

// /* initialize random seed: */
// srand (time(NULL));
// Parameters params_;

// // R2
// SourceBase<lie_groups::R2_r2> source1;
// Meas m_R2_Pos_1, m_R2_Pos_2, m_R2_Pos_Vel_1, m_R2_Pos_Vel_2;
// m_R2_Pos_1.pose = Eigen::Matrix<double,2,1>::Random();
// m_R2_Pos_2.pose = Eigen::Matrix<double,2,1>::Random();
// m_R2_Pos_1.type = MeasurementTypes::RN_POS;
// m_R2_Pos_2.type = MeasurementTypes::RN_POS;


// m_R2_Pos_Vel_1.pose = Eigen::Matrix<double,2,1>::Random();
// m_R2_Pos_Vel_1.twist = Eigen::Matrix<double,2,1>::Random();
// m_R2_Pos_Vel_2.pose = Eigen::Matrix<double,2,1>::Random();
// m_R2_Pos_Vel_2.twist = Eigen::Matrix<double,2,1>::Random();
// m_R2_Pos_Vel_1.type = MeasurementTypes::RN_POS_VEL;
// m_R2_Pos_Vel_2.type = MeasurementTypes::RN_POS_VEL;

// ASSERT_DOUBLE_EQ(source1.GetSpatialDistance( m_R2_Pos_1,     m_R2_Pos_2,params_),         (m_R2_Pos_1.pose-    m_R2_Pos_2.pose).norm());
// ASSERT_DOUBLE_EQ(source1.GetSpatialDistance( m_R2_Pos_1,     m_R2_Pos_Vel_1,params_),     (m_R2_Pos_1.pose-    m_R2_Pos_Vel_1.pose).norm());
// ASSERT_DOUBLE_EQ(source1.GetSpatialDistance( m_R2_Pos_Vel_1, m_R2_Pos_Vel_2,params_),     (m_R2_Pos_Vel_1.pose-m_R2_Pos_Vel_2.pose).norm());
// ASSERT_DOUBLE_EQ(source1.GetSpatialDistance( m_R2_Pos_Vel_1, m_R2_Pos_1,params_),         (m_R2_Pos_Vel_1.pose-m_R2_Pos_1.pose).norm());

// // R3
// SourceBase<lie_groups::R3_r3> source2;
// Meas m_R3_Pos_1, m_R3_Pos_2, m_R3_Pos_Vel_1, m_R3_Pos_Vel_2;
// m_R3_Pos_1.pose  = Eigen::Matrix<double,3,1>::Random();
// m_R3_Pos_2.pose  = Eigen::Matrix<double,3,1>::Random();
// m_R3_Pos_1.type = MeasurementTypes::RN_POS;
// m_R3_Pos_2.type = MeasurementTypes::RN_POS;


// m_R3_Pos_Vel_1.pose  = Eigen::Matrix<double,3,1>::Random();
// m_R3_Pos_Vel_1.twist = Eigen::Matrix<double,3,1>::Random();
// m_R3_Pos_Vel_2.pose  = Eigen::Matrix<double,3,1>::Random();
// m_R3_Pos_Vel_2.twist = Eigen::Matrix<double,3,1>::Random();
// m_R3_Pos_Vel_1.type = MeasurementTypes::RN_POS_VEL;
// m_R3_Pos_Vel_2.type = MeasurementTypes::RN_POS_VEL;

// ASSERT_DOUBLE_EQ(source2.GetSpatialDistance( m_R3_Pos_1,     m_R3_Pos_2,params_),         (m_R3_Pos_1.pose-    m_R3_Pos_2.pose).norm());
// ASSERT_DOUBLE_EQ(source2.GetSpatialDistance( m_R3_Pos_1,     m_R3_Pos_Vel_1,params_),     (m_R3_Pos_1.pose-    m_R3_Pos_Vel_1.pose).norm());
// ASSERT_DOUBLE_EQ(source2.GetSpatialDistance( m_R3_Pos_Vel_1, m_R3_Pos_Vel_2,params_),     (m_R3_Pos_Vel_1.pose-m_R3_Pos_Vel_2.pose).norm());
// ASSERT_DOUBLE_EQ(source2.GetSpatialDistance( m_R3_Pos_Vel_1, m_R3_Pos_1,params_),         (m_R3_Pos_Vel_1.pose-m_R3_Pos_1.pose).norm());


// }

// ///////////////////////////////////////////////////////////////////////
// //             SPATIAL DISTANCE SE2 and SE3 POS with VEL
// ///////////////////////////////////////////////////////////////////////
// TEST(SOURCE_BASE, SpatialDistance_SEN_POS) {

// /* initialize random seed: */
// srand (time(NULL));
// Parameters params_;

// // R2
// SourceBase<lie_groups::SE2_se2> source1;
// Meas m_SE2_Pos_1, m_SE2_Pos_2, m_SE2_Pos_Vel_1, m_SE2_Pos_Vel_2;
// m_SE2_Pos_1.pose  = Eigen::Matrix<double,2,1>::Random();
// m_SE2_Pos_2.pose  = Eigen::Matrix<double,2,1>::Random();
// m_SE2_Pos_1.type = MeasurementTypes::SEN_POS;
// m_SE2_Pos_2.type = MeasurementTypes::SEN_POS;


// m_SE2_Pos_Vel_1.pose  = Eigen::Matrix<double,2,1>::Random();
// m_SE2_Pos_Vel_1.twist = Eigen::Matrix<double,2,1>::Random();
// m_SE2_Pos_Vel_2.pose  = Eigen::Matrix<double,2,1>::Random();
// m_SE2_Pos_Vel_2.twist = Eigen::Matrix<double,2,1>::Random();
// m_SE2_Pos_Vel_1.type = MeasurementTypes::SEN_POS_VEL;
// m_SE2_Pos_Vel_2.type = MeasurementTypes::SEN_POS_VEL;

// ASSERT_DOUBLE_EQ(source1.GetSpatialDistance( m_SE2_Pos_1,     m_SE2_Pos_2,params_),         (m_SE2_Pos_1.pose-    m_SE2_Pos_2.pose).norm());
// ASSERT_DOUBLE_EQ(source1.GetSpatialDistance( m_SE2_Pos_1,     m_SE2_Pos_Vel_1,params_),     (m_SE2_Pos_1.pose-    m_SE2_Pos_Vel_1.pose).norm());
// ASSERT_DOUBLE_EQ(source1.GetSpatialDistance( m_SE2_Pos_Vel_1, m_SE2_Pos_Vel_2,params_),     (m_SE2_Pos_Vel_1.pose-m_SE2_Pos_Vel_2.pose).norm());
// ASSERT_DOUBLE_EQ(source1.GetSpatialDistance( m_SE2_Pos_Vel_1, m_SE2_Pos_1,params_),         (m_SE2_Pos_Vel_1.pose-m_SE2_Pos_1.pose).norm());

// // R3
// SourceBase<lie_groups::SE3_se3> source2;
// Meas m_SE3_Pos_1, m_SE3_Pos_2, m_SE3_Pos_Vel_1, m_SE3_Pos_Vel_2;
// m_SE3_Pos_1.pose  = Eigen::Matrix<double,3,1>::Random();
// m_SE3_Pos_2.pose  = Eigen::Matrix<double,3,1>::Random();
// m_SE3_Pos_1.type = MeasurementTypes::SEN_POS;
// m_SE3_Pos_2.type = MeasurementTypes::SEN_POS;


// m_SE3_Pos_Vel_1.pose  = Eigen::Matrix<double,3,1>::Random();
// m_SE3_Pos_Vel_1.twist = Eigen::Matrix<double,3,1>::Random();
// m_SE3_Pos_Vel_2.pose  = Eigen::Matrix<double,3,1>::Random();
// m_SE3_Pos_Vel_2.twist = Eigen::Matrix<double,3,1>::Random();
// m_SE3_Pos_Vel_1.type = MeasurementTypes::SEN_POS_VEL;
// m_SE3_Pos_Vel_2.type = MeasurementTypes::SEN_POS_VEL;

// ASSERT_DOUBLE_EQ(source2.GetSpatialDistance( m_SE3_Pos_1,     m_SE3_Pos_2,params_),         (m_SE3_Pos_1.pose-    m_SE3_Pos_2.pose).norm());
// ASSERT_DOUBLE_EQ(source2.GetSpatialDistance( m_SE3_Pos_1,     m_SE3_Pos_Vel_1,params_),     (m_SE3_Pos_1.pose-    m_SE3_Pos_Vel_1.pose).norm());
// ASSERT_DOUBLE_EQ(source2.GetSpatialDistance( m_SE3_Pos_Vel_1, m_SE3_Pos_Vel_2,params_),     (m_SE3_Pos_Vel_1.pose-m_SE3_Pos_Vel_2.pose).norm());
// ASSERT_DOUBLE_EQ(source2.GetSpatialDistance( m_SE3_Pos_Vel_1, m_SE3_Pos_1,params_),         (m_SE3_Pos_Vel_1.pose-m_SE3_Pos_1.pose).norm());


// }

// ///////////////////////////////////////////////////////////////////////
// //             SPATIAL DISTANCE SE2 and SE3 POSE with TWIST
// ///////////////////////////////////////////////////////////////////////

// TEST(SOURCE_BASE, SpatialDistance_SEN_POSE) {

// /* initialize random seed: */
// srand (time(NULL));
// Parameters params_;

// // R2
// SourceBase<lie_groups::SE2_se2> source1;
// Meas m_SE2_Pose_1, m_SE2_Pose_2, m_SE2_Pose_Twist_1, m_SE2_Pose_Twist_2;
// m_SE2_Pose_1.pose  = lie_groups::se2::Exp(Eigen::Matrix<double,3,1>::Random());
// m_SE2_Pose_2.pose  = lie_groups::se2::Exp(Eigen::Matrix<double,3,1>::Random());
// m_SE2_Pose_1.type = MeasurementTypes::SEN_POSE;
// m_SE2_Pose_2.type = MeasurementTypes::SEN_POSE;


// m_SE2_Pose_Twist_1.pose  = lie_groups::se2::Exp(Eigen::Matrix<double,3,1>::Random());
// m_SE2_Pose_Twist_1.twist = Eigen::Matrix<double,3,1>::Random();
// m_SE2_Pose_Twist_2.pose  = lie_groups::se2::Exp(Eigen::Matrix<double,3,1>::Random());
// m_SE2_Pose_Twist_2.twist = Eigen::Matrix<double,3,1>::Random();
// m_SE2_Pose_Twist_1.type = MeasurementTypes::SEN_POSE_TWIST;
// m_SE2_Pose_Twist_2.type = MeasurementTypes::SEN_POSE_TWIST;

// ASSERT_DOUBLE_EQ(source1.GetSpatialDistance( m_SE2_Pose_1,       m_SE2_Pose_2,params_),           (lie_groups::SE2::OMinus(m_SE2_Pose_1.pose,      m_SE2_Pose_2.pose)).norm());
// ASSERT_DOUBLE_EQ(source1.GetSpatialDistance( m_SE2_Pose_1,       m_SE2_Pose_Twist_1,params_),     (lie_groups::SE2::OMinus(m_SE2_Pose_1.pose,      m_SE2_Pose_Twist_1.pose)).norm());
// ASSERT_DOUBLE_EQ(source1.GetSpatialDistance( m_SE2_Pose_Twist_1, m_SE2_Pose_Twist_2,params_),     (lie_groups::SE2::OMinus(m_SE2_Pose_Twist_1.pose,m_SE2_Pose_Twist_2.pose)).norm());
// ASSERT_DOUBLE_EQ(source1.GetSpatialDistance( m_SE2_Pose_Twist_1, m_SE2_Pose_1,params_),           (lie_groups::SE2::OMinus(m_SE2_Pose_Twist_1.pose,m_SE2_Pose_1.pose)).norm());

// // R3
// SourceBase<lie_groups::SE3_se3> source2;
// Meas m_SE3_Pose_1, m_SE3_Pose_2, m_SE3_Pose_Twist_1, m_SE3_Pose_Twist_2;
// m_SE3_Pose_1.pose  = lie_groups::se3::Exp(Eigen::Matrix<double,6,1>::Random());
// m_SE3_Pose_2.pose  = lie_groups::se3::Exp(Eigen::Matrix<double,6,1>::Random());
// m_SE3_Pose_1.type = MeasurementTypes::SEN_POSE;
// m_SE3_Pose_2.type = MeasurementTypes::SEN_POSE;


// m_SE3_Pose_Twist_1.pose  = lie_groups::se3::Exp(Eigen::Matrix<double,6,1>::Random());
// m_SE3_Pose_Twist_1.twist = Eigen::Matrix<double,6,1>::Random();
// m_SE3_Pose_Twist_2.pose  = lie_groups::se3::Exp(Eigen::Matrix<double,6,1>::Random());
// m_SE3_Pose_Twist_2.twist = Eigen::Matrix<double,6,1>::Random();
// m_SE3_Pose_Twist_1.type = MeasurementTypes::SEN_POSE_TWIST;
// m_SE3_Pose_Twist_2.type = MeasurementTypes::SEN_POSE_TWIST;

// ASSERT_DOUBLE_EQ(source2.GetSpatialDistance( m_SE3_Pose_1,       m_SE3_Pose_2,params_),           (lie_groups::SE3::OMinus(m_SE3_Pose_1.pose,      m_SE3_Pose_2.pose)).norm());
// ASSERT_DOUBLE_EQ(source2.GetSpatialDistance( m_SE3_Pose_1,       m_SE3_Pose_Twist_1,params_),     (lie_groups::SE3::OMinus(m_SE3_Pose_1.pose,      m_SE3_Pose_Twist_1.pose)).norm());
// ASSERT_DOUBLE_EQ(source2.GetSpatialDistance( m_SE3_Pose_Twist_1, m_SE3_Pose_Twist_2,params_),     (lie_groups::SE3::OMinus(m_SE3_Pose_Twist_1.pose,m_SE3_Pose_Twist_2.pose)).norm());
// ASSERT_DOUBLE_EQ(source2.GetSpatialDistance( m_SE3_Pose_Twist_1, m_SE3_Pose_1,params_),           (lie_groups::SE3::OMinus(m_SE3_Pose_Twist_1.pose,m_SE3_Pose_1.pose)).norm());


// }

// TEST(SOURCE_BASE, R2_POS) {

// // Construct a state
// lie_groups::R2_r2 state;
// Eigen::Matrix<double,2,1> g;
// g << 1,2;
// state.g_.data_ = g;

// // Get parameters
// SourceParameters params;
// params.source_index_ = 100;
// params.type_ = MeasurementTypes::R2_POSE;
// params.expected_num_false_meas_ = 0.1;
// Eigen::Matrix2d R;
// R.setRandom();
// unsigned int id = 100;
// params.meas_cov_ = R;
// params.meas_cov_fixed_ = false;

// // Initialize source
// SourceBase source;
// const MeasurementTypes type = params.type_;
// source.Init<type>(params);

// // Construct Jacobians
// Eigen::Matrix<double,2,4> H;
// H << 1,0,0,0,0,1,0,0;

// Eigen::Matrix2d V;
// V.setIdentity();

// // Test to see that parameters were set properly
// ASSERT_EQ(params.expected_num_false_meas_, source.params_.expected_num_false_meas_);
// ASSERT_EQ(params.meas_cov_, source.params_.meas_cov_);
// ASSERT_EQ(params.meas_cov_fixed_, source.params_.meas_cov_fixed_);
// ASSERT_EQ(params.source_index_, source.params_.source_index_);
// ASSERT_EQ(params.type_, source.params_.type_);

// Test Jacobia.
// ASSERT_EQ(H, source.GetLinObsMatState<source.params_.type_>(state));
// ASSERT_EQ(V, source.GetLinObsMatSensorNoise<source.params_.type_>(state));

// source->GetEstMeas<rransac::SourceTypes::R2_POS>(state)

// Test Estimated measurement
// ASSERT_EQ(g, source.GetEstMeas<rransac::SourceTypes::R2_POS>(state));


// }


///////////////////////////////////////////////////////////////////////
//                             R2_POS_VEL
///////////////////////////////////////////////////////////////////////


// TEST(SOURCE_BASE, R2_POS_VEL) {

// // Construct a state
// lie_groups::R2_r2 state;
// Eigen::Matrix<double,2,1> g;
// Eigen::Matrix<double,2,1> u;
// g.setRandom();
// u.setRandom();
// state.g_.data_ = g;
// state.u_.data_ = u;

// // Get parameters
// rransac::SourceParameters params;
// params.source_id_ = 100;
// params.type_ = SourceTypes::R2_POS_VEL;
// params.expected_num_false_meas_ = 0.1;
// Eigen::Matrix4d R;
// R.setRandom();
// unsigned int id = 100;
// params.meas_cov_ = R;
// params.meas_cov_fixed_ = false;

// // Initialize source
// SourceBase source;
// source.Init<SourceTypes::R2_POS_VEL>(params);

// // Construct Jacobians
// Eigen::Matrix<double,4,4> H;
// H.setIdentity();

// Eigen::Matrix4d V;
// V.setIdentity();

// // Construct expected measurement
// Eigen::Matrix<double,4,1> meas;
// meas << g(0), g(1), u(0), u(1);


// // Test to see that parameters were set properly
// ASSERT_EQ(params.expected_num_false_meas_, source.params_.expected_num_false_meas_);
// ASSERT_EQ(params.meas_cov_, source.params_.meas_cov_);
// ASSERT_EQ(params.meas_cov_fixed_, source.params_.meas_cov_fixed_);
// ASSERT_EQ(id, source.params_.source_id_);
// ASSERT_EQ(SourceTypes::R2_POS_VEL, source.params_.type_);



// // Test Jacobians
// ASSERT_EQ(H, source.GetLinObsMatState<SourceTypes::R2_POS_VEL>(state));
// ASSERT_EQ(V, source.GetLinObsMatSensorNoise<SourceTypes::R2_POS_VEL>(state));

// // Test Estimated measurement
// ASSERT_EQ(meas, source.GetEstMeas<SourceTypes::R2_POS_VEL>(state));

// // #if CMAKE_BUILD_TYPE==Debug
// // std::cout << "debug" << std::endl;
// // #endif


// }


}