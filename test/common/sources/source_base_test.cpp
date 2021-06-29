#include <gtest/gtest.h>
#include <functional>
#include <stdlib.h>     /* srand, rand */
#include "memory.h"

#include "rransac/common/sources/source_base.h"
#include "rransac/common/transformations/transformation_base.h"
#include "rransac/common/measurement/measurement_base.h"

namespace rransac {

template <typename State>
struct CallbackClass {

bool func(const State& state) {
    if (state.g_.data_.norm() < 5)
        return true;
    else 
        return false;
}

};

//////////////////////////////////////////////////
// Dummy Source needed to initialize SourceBase
/////////////////////////////////////////////////

template<class tState>
class TransformDummy : public TransformBase<Eigen::Matrix<typename tState::DataType,Eigen::Dynamic, Eigen::Dynamic>, tState, Eigen::Matrix<typename tState::DataType,tState::dim_,tState::dim_>,false,TransformDummy<tState>> {
public:


typedef Eigen::Matrix<typename tState::DataType,Eigen::Dynamic, Eigen::Dynamic> MatXd;
typedef MatXd Data;
typedef typename tState::DataType DataType;  /**< The scalar object for the data. Ex. float, double, etc. */
typedef tState State;
typedef Eigen::Matrix<typename tState::DataType, tState::dim_, tState::dim_> MatCov;

void DerivedInit(){}

void DerivedSetData(const Data& data){}

void GetData() const {}

void DerivedTransformMeasurement(Meas<DataType>& meas) const {}

void DerivedTransformTrack(State& state, MatCov& cov) const {}

static void DerivedTransformTrack(State& state, MatCov& cov, const Data& transform_data) {
   state.g_.data_ = transform_data * state.g_.data_;
   state.u_.data_ *=2;
   cov *=2;
}

static State DerivedTransformState(const State& state, const Data& transform_data) {
   State state_transformed = state;
   state_transformed.g_.data_ = transform_data * state.g_.data_;
   state_transformed.u_.data_ *=2;
   return state_transformed;
}

static MatXd DerivedGetTransformationJacobian(const State& state, const Data& transform_data) {
   Eigen::Matrix<DataType, tState::dim_,tState::dim_> J;
   J.setIdentity();
   J *= transform_data(0,0);
   return J;
}

};



//////////////////////////////////////////////////
// Dummy Source needed to initialize SourceBase
/////////////////////////////////////////////////
template <typename tState, MeasurementTypes tMeasurementType, template <typename > typename tTransformation>
class DummySource : public SourceBase<tState, tMeasurementType, tTransformation,tState::Group::dim_, DummySource> 
{
public:

    using State = tState;
    typedef tTransformation<State> Transformation;
    typedef typename tState::DataType DataType;
    typedef Eigen::Matrix<DataType,Eigen::Dynamic,Eigen::Dynamic> MatXd;
    static constexpr unsigned int meas_space_dim_ = tState::Group::dim_;
    static constexpr MeasurementTypes meas_type_ = tMeasurementType;
    typedef SourceBase<tState, tMeasurementType, tTransformation, tState::Group::dim_, DummySource> Base;

    static constexpr int dim_mult_ = MeasHasVelocity<tMeasurementType>::value ? 2 : 1;

    // static MatXd H_;
    // static MatXd V_;

    /** Initializes the measurement source. This function must set the parameters.  */
    void DerivedInit(const SourceParameters& params) { 
        Base::H_ = Eigen::Matrix<double,dim_mult_*meas_space_dim_,tState::Group::dim_*2>::Identity();
        Base::V_ = Eigen::Matrix<double,dim_mult_*meas_space_dim_,dim_mult_*meas_space_dim_>::Identity();
    }


    /** Returns the jacobian of the observation function w.r.t. the states */
    static Eigen::MatrixXd DerivedGetLinObsMatState(const State& state)  {
        return Base::H_*state.g_.data_(0,0);
    }                              

    /** Returns the jacobian of the observation function w.r.t. the sensor noise */
    static Eigen::MatrixXd DerivedGetLinObsMatSensorNoise(const State& state)  {
        return Base::V_*state.g_.data_(0,0);
    }                         

    /** Computes the estimated measurement given a state */
    static Meas<DataType> DerivedGetEstMeas(const State& state) {
        Meas<DataType> tmp;
        tmp.pose = state.g_.data_;
        return tmp;
    } /** Returns an estimated measurement according to the state. */


    /**
     * Performs the OMinus operation between two measurement (m1 ominus m2) of the same type. In other words, this
     * method computes the geodesic distance between two measurements of the same type.
     */ 
    static MatXd DerivedOMinus(const Meas<DataType>& m1, const Meas<DataType>& m2) {
        Eigen::Matrix<double,meas_space_dim_,1> err = m1.pose - m2.pose;
        return err;
    } 

   /**
     * Generates a random measurement from a Gaussian distribution with mean defined by the state and standard deviation defined by meas_std. This
     * method is used primarily in simulations and tests.
     */ 
    Meas<DataType> DerivedGenerateRandomMeasurement(const MatXd& meas_std, const State& state) const {
        Meas<DataType> m;
        m.type = tMeasurementType;
        m.source_index = this->params_.source_index_;
        m.pose = Eigen::Matrix<double,meas_space_dim_,1>::Ones()*state.g_.data_(0,0);
        return m;

    }

};

typedef DummySource<lie_groups::R2_r2, MeasurementTypes::RN_POS, TransformDummy> DummySource1;
typedef DummySource<lie_groups::R2_r2, MeasurementTypes::RN_POS_VEL, TransformDummy> DummySource2;
typedef DummySource<lie_groups::R3_r3, MeasurementTypes::RN_POS_VEL, TransformDummy> DummySource3;


///////////////////////////////////////////////////////////////////////
//                             Test the Init functions and state withing surveillance region
///////////////////////////////////////////////////////////////////////
TEST(SOURCE_BASE, InitFunction) {

SourceParameters source_params;
DummySource1 source;
CallbackClass<lie_groups::R2_r2> call;


//
// Invalid source parameters
//

// Empty measurement covariance 
source_params.spacial_density_of_false_meas_ = 0.1;
source_params.type_ = MeasurementTypes::RN_POS;
source_params.gate_probability_ = 0.8;
source_params.probability_of_detection_ = 0.9;
ASSERT_ANY_THROW(source.Init(source_params, std::bind(&CallbackClass<lie_groups::R2_r2>::func, call, std::placeholders::_1)));

// Measurement covariance of wrong dimension
source_params.meas_cov_ = Eigen::Matrix3d::Identity();
ASSERT_ANY_THROW(source.Init(source_params, std::bind(&CallbackClass<lie_groups::R2_r2>::func, call, std::placeholders::_1)));

// not symmetric measurement covariance
source_params.meas_cov_ = Eigen::Matrix2d::Identity();
source_params.meas_cov_ << 0, 1, 2, 0;
ASSERT_ANY_THROW(source.Init(source_params, std::bind(&CallbackClass<lie_groups::R2_r2>::func, call, std::placeholders::_1)));

// not positive definite measurement covariance
source_params.meas_cov_ << -1, 0, 0, -1;
ASSERT_ANY_THROW(source.Init(source_params, std::bind(&CallbackClass<lie_groups::R2_r2>::func, call, std::placeholders::_1)));

// Invalid spacial density
source_params.meas_cov_ = Eigen::Matrix2d::Identity();
source_params.spacial_density_of_false_meas_ = -0.1;
ASSERT_ANY_THROW(source.Init(source_params, std::bind(&CallbackClass<lie_groups::R2_r2>::func, call, std::placeholders::_1)));
source_params.spacial_density_of_false_meas_ = 0.9;

// Invalid measurement type
source_params.type_ = MeasurementTypes::NUM_TYPES;
ASSERT_ANY_THROW(source.Init(source_params, std::bind(&CallbackClass<lie_groups::R2_r2>::func, call, std::placeholders::_1)));
source_params.type_ = MeasurementTypes::RN_POS;

// Invalid gate_probability_
source_params.gate_probability_ = -0.1;
ASSERT_ANY_THROW(source.Init(source_params, std::bind(&CallbackClass<lie_groups::R2_r2>::func, call, std::placeholders::_1)));
source_params.gate_probability_ = 1.1;
ASSERT_ANY_THROW(source.Init(source_params, std::bind(&CallbackClass<lie_groups::R2_r2>::func, call, std::placeholders::_1)));
source_params.gate_probability_ = 0.8;

// Invalid probability of detection
source_params.probability_of_detection_ = -0.1;
ASSERT_ANY_THROW(source.Init(source_params, std::bind(&CallbackClass<lie_groups::R2_r2>::func, call, std::placeholders::_1)));
source_params.probability_of_detection_ = 1.1;
ASSERT_ANY_THROW(source.Init(source_params, std::bind(&CallbackClass<lie_groups::R2_r2>::func, call, std::placeholders::_1)));
source_params.probability_of_detection_ = 0.9;



// Valid source parameters
source_params.meas_cov_ = Eigen::Matrix2d::Identity();
source_params.spacial_density_of_false_meas_ = 0.1;
source_params.type_ = MeasurementTypes::RN_POS;
source_params.gate_probability_ = 0.393469340287367;
source_params.probability_of_detection_ = 0.9;

ASSERT_NO_THROW(source.Init(source_params, std::bind(&CallbackClass<lie_groups::R2_r2>::func, call, std::placeholders::_1)));

// Check that the dimesions were set.
int meas_space_dim = DummySource1::meas_space_dim_;
int meas_space_dim_mult = DummySource1::dim_mult_;
ASSERT_EQ(source.GetParams().meas_space_dim_,meas_space_dim);
ASSERT_EQ(source.GetParams().meas_space_dim_mult_,meas_space_dim_mult);

// Check the gate threshold and unit hyptersphere.
ASSERT_LE( fabs(1- source.GetParams().gate_threshold_), 1e-6);
ASSERT_LE( fabs(1- source.GetParams().gate_threshold_sqrt_), 1e-4);
ASSERT_LE( fabs(M_PI- source.GetParams().vol_unit_hypershpere_  ), 1e-6);

DummySource2 source2;

source_params.type_ = MeasurementTypes::RN_POS_VEL;
source_params.meas_cov_ = Eigen::Matrix4d::Identity();
source_params.gate_probability_ = 0.593994150290162;
ASSERT_NO_THROW(source2.Init(source_params, std::bind(&CallbackClass<lie_groups::R2_r2>::func, call, std::placeholders::_1)));
ASSERT_LE( fabs(4- source2.GetParams().gate_threshold_), 1e-3);
ASSERT_LE( fabs(2- source2.GetParams().gate_threshold_sqrt_), 1e-3);
ASSERT_LE( fabs(4.934802200544679- source2.GetParams().vol_unit_hypershpere_  ), 1e-3);

meas_space_dim = DummySource2::meas_space_dim_;
meas_space_dim_mult = DummySource2::dim_mult_;
ASSERT_EQ(source2.GetParams().meas_space_dim_,meas_space_dim);
ASSERT_EQ(source2.GetParams().meas_space_dim_mult_,meas_space_dim_mult);

// State is within surveillance region.
Eigen::MatrixXd EmptyMat;
lie_groups::R2_r2 state;
state.g_.data_ << 1,1;
ASSERT_TRUE(source.StateInsideSurveillanceRegion(state, false, EmptyMat));

// Going to transform the data so that it is outside of the surveillance region.
state.g_.data_ << 3,3;
Eigen::MatrixXd transform_data = Eigen::Matrix2d::Identity()*2;
ASSERT_FALSE(source.StateInsideSurveillanceRegion(state,true, transform_data));


// State outside of surveillance region
state.g_.data_ << 5,5;
ASSERT_FALSE(source.StateInsideSurveillanceRegion(state,false, EmptyMat));

// Test default StateInsideSurveillanceRegionDefaultCallback
source_params.type_ = MeasurementTypes::RN_POS;
source_params.meas_cov_ = Eigen::Matrix2d::Identity();
source.Init(source_params);
ASSERT_TRUE(source.StateInsideSurveillanceRegion(state, false, EmptyMat));


}

      

///////////////////////////////////////////////////////////////////////
//                            Distances, Jacobians, Est Meas, O Minus, Generate Random Measurements
///////////////////////////////////////////////////////////////////////

TEST(SOURCE_BASE, TemporalDistance) {




/* initialize random seed: */
srand (time(NULL));
SourceParameters source_params_;

source_params_.meas_cov_ = Eigen::Matrix2d::Identity();
source_params_.spacial_density_of_false_meas_ = 0.1;
source_params_.type_ = MeasurementTypes::RN_POS;
source_params_.gate_probability_ = 0.8;
source_params_.probability_of_detection_ = 0.9;
Parameters params_;
DummySource1 source;

lie_groups::R2_r2 state;
state.g_.data_.setRandom();

Meas<double> m1, m2;
m1.time_stamp = rand() % 100 -50;
m2.time_stamp = rand() % 100 -50;

// Test init function
source.Init(source_params_);
ASSERT_EQ(source.GetParams().type_, source_params_.type_);

Eigen::MatrixXd EmptyMat;
// Test GetLinObsMatState
ASSERT_EQ(source.GetLinObsMatState(state, false, EmptyMat), (DummySource1::H_*state.g_.data_(0,0)));

// Test GetLinObsMatSensorNoise
ASSERT_EQ(source.GetLinObsMatSensorNoise(state, false, EmptyMat), DummySource1::V_*state.g_.data_(0,0));

// Test GetEstMeas
ASSERT_EQ(source.GetEstMeas(state, false, EmptyMat).pose, state.g_.data_);

// Test the transformation
Eigen::MatrixXd transform_data = Eigen::Matrix2d::Identity()*2;
typename DummySource1::State state_transformed = state;
state_transformed.g_.data_ = transform_data * state_transformed.g_.data_;
ASSERT_EQ(source.GetLinObsMatState(state, true, transform_data), DummySource1::H_*state_transformed.g_.data_(0,0)*transform_data(0,0));
ASSERT_EQ(source.GetLinObsMatSensorNoise(state, true, transform_data), DummySource1::V_*state_transformed.g_.data_(0,0));
ASSERT_EQ(source.GetEstMeas(state, true, transform_data).pose, state_transformed.g_.data_);



m1.pose = state.g_.data_;
m2.pose = state.g_.data_ *2;

ASSERT_EQ(source.OMinus(m1,m2), m1.pose - m2.pose);

ASSERT_EQ(source.GenerateRandomMeasurement(EmptyMat, state, false, EmptyMat).pose, (Eigen::Matrix<double,2,1>::Ones()*state.g_.data_(0,0))  );
ASSERT_EQ(source.GenerateRandomMeasurement(EmptyMat, state, true, transform_data).pose, (Eigen::Matrix<double,2,1>::Ones()*state_transformed.g_.data_(0,0)));





// lie_groups::R2_r2 state;
// SourceParameters params;
// source->Init(params);

ASSERT_EQ(source.GetTemporalDistance(m1,m2,params_),fabs(m1.time_stamp-m2.time_stamp));


// // Random number generator;
// int num_randn = 50000;
// std::vector<Eigen::Matrix<double,4,1>> randn_nums(num_randn);
// Eigen::Matrix<double,4,1> mean;
// Eigen::Matrix<double,4,4> cov;
// mean.setZero();
// cov.setZero();

// // Get random numbers and calculate mean
// for(Eigen::Matrix<double,4,1>& randn_num : randn_nums) {
//     randn_num = source.GaussianRandomGenerator(4);
//     mean += randn_num;
//     // std::cout << randn_num << std::endl << std::endl;
// }

// mean /= num_randn;

// // Calculate std
// for(Eigen::Matrix<double,4,1> randn_num : randn_nums) {
//     cov += (mean - randn_num)*(mean-randn_num).transpose();
    
// }

// cov /= num_randn;

// ASSERT_LE( (mean - Eigen::Matrix<double,4,1>::Zero()).norm(), 0.1);
// ASSERT_LE( (cov - Eigen::Matrix<double,4,4>::Identity()).norm(), 0.1);


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

DummySource2 source1;
Meas<double> m_R2_Pos_1, m_R2_Pos_2, m_R2_Pos_Vel_1, m_R2_Pos_Vel_2;
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

DummySource2 source1;
Meas<double> m_R2_Pos_1, m_R2_Pos_2, m_R2_Pos_Vel_1, m_R2_Pos_Vel_2;
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


// se2
typedef lie_groups::SE2_se2 State_SE2;
DummySource<lie_groups::SE2_se2, MeasurementTypes::SEN_POSE, TransformDummy> source1;
Meas<double> m_SE2_Pose_1, m_SE2_Pose_2, m_SE2_Pose_Twist_1, m_SE2_Pose_Twist_2;
m_SE2_Pose_1.pose  = State_SE2::Algebra::Exp(Eigen::Matrix<double,3,1>::Random());
m_SE2_Pose_2.pose  = State_SE2::Algebra::Exp(Eigen::Matrix<double,3,1>::Random());
m_SE2_Pose_1.type = MeasurementTypes::SEN_POSE;
m_SE2_Pose_2.type = MeasurementTypes::SEN_POSE;


m_SE2_Pose_Twist_1.pose  = State_SE2::Algebra::Exp(Eigen::Matrix<double,3,1>::Random());
m_SE2_Pose_Twist_1.twist = Eigen::Matrix<double,3,1>::Random();
m_SE2_Pose_Twist_2.pose  = State_SE2::Algebra::Exp(Eigen::Matrix<double,3,1>::Random());
m_SE2_Pose_Twist_2.twist = Eigen::Matrix<double,3,1>::Random();
m_SE2_Pose_Twist_1.type = MeasurementTypes::SEN_POSE_TWIST;
m_SE2_Pose_Twist_2.type = MeasurementTypes::SEN_POSE_TWIST;

ASSERT_DOUBLE_EQ(source1.GetSpatialDistance( m_SE2_Pose_1,       m_SE2_Pose_2,params_),           (State_SE2::Group::OMinus(m_SE2_Pose_1.pose,      m_SE2_Pose_2.pose)).norm());
ASSERT_DOUBLE_EQ(source1.GetSpatialDistance( m_SE2_Pose_1,       m_SE2_Pose_Twist_1,params_),     (State_SE2::Group::OMinus(m_SE2_Pose_1.pose,      m_SE2_Pose_Twist_1.pose)).norm());
ASSERT_DOUBLE_EQ(source1.GetSpatialDistance( m_SE2_Pose_Twist_1, m_SE2_Pose_Twist_2,params_),     (State_SE2::Group::OMinus(m_SE2_Pose_Twist_1.pose,m_SE2_Pose_Twist_2.pose)).norm());
ASSERT_DOUBLE_EQ(source1.GetSpatialDistance( m_SE2_Pose_Twist_1, m_SE2_Pose_1,params_),           (State_SE2::Group::OMinus(m_SE2_Pose_Twist_1.pose,m_SE2_Pose_1.pose)).norm());


typedef lie_groups::SE3_se3 State_SE3;
DummySource<lie_groups::SE3_se3, MeasurementTypes::SEN_POSE, TransformDummy> source2;
Meas<double> m_SE3_Pose_1, m_SE3_Pose_2, m_SE3_Pose_Twist_1, m_SE3_Pose_Twist_2;
m_SE3_Pose_1.pose  = State_SE3::Algebra::Exp(Eigen::Matrix<double,6,1>::Random());
m_SE3_Pose_2.pose  = State_SE3::Algebra::Exp(Eigen::Matrix<double,6,1>::Random());
m_SE3_Pose_1.type = MeasurementTypes::SEN_POSE;
m_SE3_Pose_2.type = MeasurementTypes::SEN_POSE;


m_SE3_Pose_Twist_1.pose  = State_SE3::Algebra::Exp(Eigen::Matrix<double,6,1>::Random());
m_SE3_Pose_Twist_1.twist = Eigen::Matrix<double,6,1>::Random();
m_SE3_Pose_Twist_2.pose  = State_SE3::Algebra::Exp(Eigen::Matrix<double,6,1>::Random());
m_SE3_Pose_Twist_2.twist = Eigen::Matrix<double,6,1>::Random();
m_SE3_Pose_Twist_1.type = MeasurementTypes::SEN_POSE_TWIST;
m_SE3_Pose_Twist_2.type = MeasurementTypes::SEN_POSE_TWIST;

ASSERT_DOUBLE_EQ(source2.GetSpatialDistance( m_SE3_Pose_1,       m_SE3_Pose_2,params_),           (State_SE3::Group::OMinus(m_SE3_Pose_1.pose,      m_SE3_Pose_2.pose)).norm());
ASSERT_DOUBLE_EQ(source2.GetSpatialDistance( m_SE3_Pose_1,       m_SE3_Pose_Twist_1,params_),     (State_SE3::Group::OMinus(m_SE3_Pose_1.pose,      m_SE3_Pose_Twist_1.pose)).norm());
ASSERT_DOUBLE_EQ(source2.GetSpatialDistance( m_SE3_Pose_Twist_1, m_SE3_Pose_Twist_2,params_),     (State_SE3::Group::OMinus(m_SE3_Pose_Twist_1.pose,m_SE3_Pose_Twist_2.pose)).norm());
ASSERT_DOUBLE_EQ(source2.GetSpatialDistance( m_SE3_Pose_Twist_1, m_SE3_Pose_1,params_),           (State_SE3::Group::OMinus(m_SE3_Pose_Twist_1.pose,m_SE3_Pose_1.pose)).norm());


}



} // namespace rransac