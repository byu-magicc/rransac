#include <gtest/gtest.h>
#include <Eigen/Core>
#include "rransac/common/measurement/measurement_base.h"

#include "rransac/common/models/model_RN.h"
#include "rransac/common/models/model_SEN_pos_vel.h"
#include "rransac/common/models/model_SEN_pose_twist.h"
#include "rransac/common/sources/source_RN.h"
#include "rransac/common/sources/source_SEN_pos_vel.h"
#include "rransac/common/sources/source_SEN_pose_twist.h"
#include "rransac/common/sources/source_container.h"
#include "rransac/parameters.h"
#include "rransac/common/transformations/transformation_null.h"
#include "rransac/system.h"

#include "rransac/common/data_association/validation_region_policies/validation_region_fixed_policy.h"
#include "rransac/common/data_association/validation_region_policies/validation_region_fixed_pos_policy.h"
#include "rransac/common/data_association/validation_region_policies/validation_region_innov_policy.h"

using namespace lie_groups;
using namespace rransac;

typedef SourceRN<R2_r2,MeasurementTypes::RN_POS,TransformNULL> SourceR2Pos;
typedef SourceRN<R2_r2,MeasurementTypes::RN_POS_VEL,TransformNULL> SourceR2PosVel;
typedef SourceRN<R3_r3,MeasurementTypes::RN_POS,TransformNULL> SourceR3Pos;
typedef SourceRN<R3_r3,MeasurementTypes::RN_POS_VEL,TransformNULL> SourceR3PosVel;

typedef SourceSENPosVel<SE2_se2,MeasurementTypes::SEN_POS,TransformNULL> SourceSE2Pos;
typedef SourceSENPosVel<SE2_se2,MeasurementTypes::SEN_POS_VEL,TransformNULL> SourceSE2PosVel;
typedef SourceSENPosVel<SE3_se3,MeasurementTypes::SEN_POS,TransformNULL> SourceSE3Pos;
typedef SourceSENPosVel<SE3_se3,MeasurementTypes::SEN_POS_VEL,TransformNULL> SourceSE3PosVel;

typedef SourceContainer<SourceR2Pos,SourceR2PosVel> SourceContainerR2;
typedef SourceContainer<SourceR3Pos,SourceR3PosVel> SourceContainerR3;
typedef SourceContainer<SourceSE2Pos,SourceSE2PosVel> SourceContainerSE2;
typedef SourceContainer<SourceSE3Pos,SourceSE3PosVel> SourceContainerSE3;



typedef ModelRN<SourceContainerR2> Model1;
typedef ModelRN<SourceContainerR3> Model2;
typedef ModelSENPosVel<SourceContainerSE2> Model3;
typedef ModelSENPosVel<SourceContainerSE3> Model4;

// N is the dimension of the position
template<class M, MeasurementTypes MT1, MeasurementTypes MT2, int N>
struct ModelHelper {
    typedef M Model;
    static constexpr MeasurementTypes MeasType1 = MT1;
    static constexpr MeasurementTypes MeasType2 = MT2;
    typedef typename Model::Measurement Measurement;

    static Eigen::MatrixXd GetPosError(const M& track, const Measurement& meas) {
        return meas.pose.block(0,0,N,1) - track.state_.g_.data_.block(0,track.state_.g_.data_.cols()-1,N,1);
    }

    static Eigen::MatrixXd GetError(const M& track, const Measurement& meas, const typename M::SourceContainer& source_container) {
        bool transform_state = false;
        Eigen::MatrixXd EmptyMat;
        return source_container.OMinus(meas.source_index, meas, source_container.GetEstMeas(meas.source_index, track.state_, transform_state, EmptyMat));
    }
};





typedef ModelHelper<Model1, MeasurementTypes::RN_POS, MeasurementTypes::RN_POS_VEL,2> ModelHelper1;
typedef ModelHelper<Model2, MeasurementTypes::RN_POS, MeasurementTypes::RN_POS_VEL,3> ModelHelper2;
typedef ModelHelper<Model3, MeasurementTypes::SEN_POS, MeasurementTypes::SEN_POS_VEL,2> ModelHelper3;
typedef ModelHelper<Model4, MeasurementTypes::SEN_POS, MeasurementTypes::SEN_POS_VEL,3> ModelHelper4;


using MyTypes = ::testing::Types<ModelHelper1, ModelHelper2, ModelHelper3, ModelHelper4>;
// using MyTypes = ::testing::Types<ModelHelper3>;

template<class ModelHelper>
class ValidationRegionTest : public ::testing::Test {

protected:



static constexpr unsigned int meas_dim = ModelHelper::Model::SourceContainer::Source0::meas_pose_dim_;
static constexpr unsigned int state_dim = ModelHelper::Model::State::g_type_::dim_*2;
static constexpr unsigned int cov_dim = ModelHelper::Model::cov_dim_;
static constexpr unsigned int a_vel_dim = ModelHelper::Model::cov_dim_ - ModelHelper::Model::State::g_type_::dim_-1;
static constexpr unsigned int t_vel_dim = state_dim - ModelHelper::Model::State::g_type_::dim_ - a_vel_dim;


typedef typename ModelHelper::Model::Base::Measurement Measurement;

void SetUp() override {






double meas_cov_scale = 0.1;
double system_cov_scale = 0.1;
meas_cov1.setIdentity();
meas_cov1 *= meas_cov_scale;
meas_std1.setIdentity();
meas_std1 *= sqrt(meas_cov_scale);

source_params1.meas_cov_ = meas_cov1;
source_params1.type_ = ModelHelper::MeasType1;
source_params1.source_index_ = 0;
source_params1.spacial_density_of_false_meas_ = 0.8;
source_params1.probability_of_detection_ = 0.8;
source_params1.gate_probability_ = 0.8;
source_params1.gate_threshold_ = 1;

meas_cov2.setIdentity();
meas_cov2 *= meas_cov_scale;
// source_params2.meas_cov_fixed_ = true;
// source_params2.meas_cov_ = meas_cov2;
source_params2.type_ = ModelHelper::MeasType2;
source_params2.source_index_ = 1;
source_params2.spacial_density_of_false_meas_ = 0.8;
source_params2.probability_of_detection_ = 0.8;
source_params2.gate_probability_ = 0.8;
source_params2.meas_cov_ = meas_cov2;
source_params1.gate_threshold_ = 1;
meas_std2.setIdentity();
meas_std2 *= sqrt(meas_cov_scale);


// std::cerr << "here00" << std::endl;



sys.params_.process_noise_covariance_ = Eigen::Matrix<double,cov_dim,cov_dim>::Identity()*system_cov_scale;
sys.source_container_.AddSource(source_params1);
sys.source_container_.AddSource(source_params2);
// std::cerr << "here 1" << std::endl;

// Generate Trajectory and measurements
Measurement m1;
Measurement m2;

track.Init(sys.params_);
track.state_ = ModelHelper::Model::State::Random();




bool transform_state = false;
Eigen::MatrixXd EmptyMat;




for (unsigned int ii = 0; ii < num_meas; ++ii) {


    m1 = sys.source_container_.GenerateRandomMeasurement(0,meas_std1,track.state_,transform_state,EmptyMat);
    m2 = sys.source_container_.GenerateRandomMeasurement(1,meas_std2,track.state_,transform_state,EmptyMat);
    m1.time_stamp = 0;
    m1.source_index = 0;
    m1.type = ModelHelper::MeasType1;
    m2.time_stamp = 0;
    m2.source_index = 1;
    m2.type = ModelHelper::MeasType2;
        
    sys.new_meas_.push_back(m1);
    sys.new_meas_.push_back(m2);

}



}

//-----------------------------------------------------------------------------------------------
typedef typename ModelHelper::Model model;
SourceParameters source_params1;
SourceParameters source_params2;
Eigen::Matrix<double,meas_dim, meas_dim> meas_cov1;
Eigen::Matrix<double,meas_dim*2, meas_dim*2> meas_cov2;
Eigen::Matrix<double,meas_dim, meas_dim> meas_std1;
Eigen::Matrix<double,meas_dim*2, meas_dim*2> meas_std2;
std::vector<std::vector<std::vector<Measurement>>> new_meas; // time, source, measurements
System<model> sys;

int num_meas = 1000;

typename ModelHelper::Model track;

};


TYPED_TEST_SUITE(ValidationRegionTest, MyTypes);

//-------------------------------------------------------------------------------------------------------------------------------

TYPED_TEST(ValidationRegionTest,ValidationRegionFixedPolicy ) {


    ValidationRegionFixedPolicy<typename TypeParam::ModelHelper::Model> validation_region;


    bool in_validation_region = false;
    Eigen::MatrixXd err;
    int num_true = 0;




    for (auto& meas: this->sys.new_meas_) {


        err = TypeParam::ModelHelper::GetError(this->track,meas,this->sys.source_container_);

        if(err.norm() <= this->sys.source_container_.GetParams(meas.source_index).gate_threshold_) {
            in_validation_region = true;
            num_true++;
        } else {
            in_validation_region = false;
            // std::cout << "err: " << std::endl <<  err << std::endl;
            // std::cout << "err norm: "  <<  err.norm() << std::endl;
            // std::cout << "meas: " << std::endl <<  meas.pose << std::endl;
            // std::cout << "track: " << std::endl <<  this->track.state_.g_.data_ << std::endl;
        }

        ASSERT_EQ(in_validation_region, validation_region.PolicyInValidationRegion(this->sys,meas,this->track));
        EXPECT_GT(num_true,0) << "Try running the test case again" << std::endl;


    }
}

//-------------------------------------------------------------------------------------------------------------------------------

TYPED_TEST(ValidationRegionTest,ValidationRegionFixedPosPolicy ) {


    ValidationRegionFixedPosPolicy<typename TypeParam::ModelHelper::Model> validation_region;


    bool in_validation_region = false;
    Eigen::MatrixXd err;

    int num_true = 0;


    for (auto& meas: this->sys.new_meas_) {


        err = TypeParam::ModelHelper::GetPosError(this->track,meas);

        if(err.norm() <= this->sys.source_container_.GetParams(meas.source_index).gate_threshold_) {
            in_validation_region = true;
            num_true++;
        } else {
            in_validation_region = false;
            // std::cout << "err: " << std::endl <<  err << std::endl;
            // std::cout << "err norm: "  <<  err.norm() << std::endl;
            // std::cout << "meas: " << std::endl <<  meas.pose << std::endl;
            // std::cout << "track: " << std::endl <<  this->track.state_.g_.data_ << std::endl;
        }

        ASSERT_EQ(in_validation_region, validation_region.PolicyInValidationRegion(this->sys,meas,this->track));
        EXPECT_GT(num_true,0) << "Try running the test case again" << std::endl;

    }
}

//-------------------------------------------------------------------------------------------------------------------------------

TYPED_TEST(ValidationRegionTest,ValidationRegionInnovPolicy ) {


    ValidationRegionInnovPolicy<typename TypeParam::ModelHelper::Model> validation_region;


    bool in_validation_region = false;
    Eigen::MatrixXd err;
    Eigen::MatrixXd S;
    bool transform_state = false;
    Eigen::MatrixXd EmptyMat;

    int num_true = 0;





    for (auto& meas: this->sys.new_meas_) {


        err = TypeParam::ModelHelper::GetError(this->track,meas, this->sys.source_container_);
        S = this->track.GetInnovationCovariance(this->sys.source_container_, meas.source_index, transform_state, EmptyMat);

        // std::cout << "err: " << std::endl << err << std::endl;
        // std::cout << "S: " << std::endl << S << std::endl;

        if((err.transpose()*S.inverse()*err)(0,0) <= this->sys.source_container_.GetParams(meas.source_index).gate_threshold_) {
            in_validation_region = true;
            num_true++;
        } else {
            in_validation_region = false;
            // std::cout << "err: " << std::endl <<  err << std::endl;
            // std::cout << "err norm: "  <<  err.norm() << std::endl;
            // std::cout << "meas: " << std::endl <<  meas.pose << std::endl;
            // std::cout << "track: " << std::endl <<  this->track.state_.g_.data_ << std::endl;
        }

        ASSERT_EQ(in_validation_region, validation_region.PolicyInValidationRegion(this->sys,meas,this->track));
        EXPECT_GT(num_true,0) << "Try running the test case again" << std::endl;


    }
}





