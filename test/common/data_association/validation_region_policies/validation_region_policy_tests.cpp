#include <gtest/gtest.h>
#include <Eigen/Core>
#include "rransac/common/measurement/measurement_base.h"

#include "rransac/common/models/model_RN.h"
#include "rransac/common/models/model_SEN_pos_vel.h"
#include "rransac/common/models/model_SEN_pose_twist.h"
#include "rransac/common/sources/source_RN.h"
#include "rransac/common/sources/source_SEN_pos_vel.h"
#include "rransac/common/sources/source_SEN_pose_twist.h"
#include "rransac/parameters.h"
#include "rransac/common/transformations/transformation_null.h"
#include "rransac/system.h"

#include "rransac/common/data_association/validation_region_policies/validation_region_fixed_policy.h"
#include "rransac/common/data_association/validation_region_policies/validation_region_fixed_pos_policy.h"
#include "rransac/common/data_association/validation_region_policies/validation_region_innov_policy.h"

using namespace lie_groups;
using namespace rransac;

typedef ModelRN<R2_r2, TransformNULL, SourceRN> Model1;
typedef ModelRN<R3_r3, TransformNULL, SourceRN> Model2;
typedef ModelSENPosVel<SE2_se2, TransformNULL, SourceSENPosVel> Model3;
typedef ModelSENPosVel<SE3_se3, TransformNULL, SourceSENPosVel> Model4;

// N is the dimension of the position
template<class M, MeasurementTypes MT1, MeasurementTypes MT2, int N>
struct ModelHelper {
    typedef M Model;
    static constexpr MeasurementTypes MeasType1 = MT1;
    static constexpr MeasurementTypes MeasType2 = MT2;

    static Eigen::MatrixXd GetPosError(const M& track, const Meas<typename M::DataType>& meas) {
        return meas.pose.block(0,0,N,1) - track.state_.g_.data_.block(0,track.state_.g_.data_.cols()-1,N,1);
    }

    static Eigen::MatrixXd GetError(const M& track, const Meas<typename M::DataType>& meas, const std::vector<typename M::Source>& sources) {
        return sources[meas.source_index].OMinus(meas, sources[meas.source_index].GetEstMeas(track.state_));
    }
};





typedef ModelHelper<Model1, MeasurementTypes::RN_POS, MeasurementTypes::RN_POS_VEL,2> ModelHelper1;
typedef ModelHelper<Model2, MeasurementTypes::RN_POS, MeasurementTypes::RN_POS_VEL,3> ModelHelper2;
typedef ModelHelper<Model3, MeasurementTypes::SEN_POS, MeasurementTypes::SEN_POS_VEL,2> ModelHelper3;
typedef ModelHelper<Model4, MeasurementTypes::SEN_POS, MeasurementTypes::SEN_POS_VEL,3> ModelHelper4;


using MyTypes = ::testing::Types<ModelHelper1, ModelHelper2, ModelHelper3, ModelHelper4>;
// using MyTypes = ::testing::Types<ModelHelper1>;

template<class ModelHelper>
class ValidationRegionTest : public ::testing::Test {

protected:

typedef Meas<double> Measurement;

static constexpr unsigned int meas_dim = ModelHelper::Model::Source::meas_space_dim_;
static constexpr unsigned int state_dim = ModelHelper::Model::State::g_type_::dim_*2;
static constexpr unsigned int cov_dim = ModelHelper::Model::cov_dim_;
static constexpr unsigned int a_vel_dim = ModelHelper::Model::cov_dim_ - ModelHelper::Model::State::g_type_::dim_-1;
static constexpr unsigned int t_vel_dim = state_dim - ModelHelper::Model::State::g_type_::dim_ - a_vel_dim;

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

// std::cerr << "here0" << std::endl;
typename ModelHelper::Model::Source source1;
typename ModelHelper::Model::Source source2;
source1.Init(source_params1);
source2.Init(source_params2);
// std::cerr << "here00" << std::endl;



sys.params_.process_noise_covariance_ = Eigen::Matrix<double,cov_dim,cov_dim>::Identity()*system_cov_scale;
sys.sources_.push_back(source1);
sys.sources_.push_back(source2);
// std::cerr << "here 1" << std::endl;

// Generate Trajectory and measurements
Measurement m1;
Measurement m2;

track.Init(sys.params_,sys.sources_.size());
track.state_ = ModelHelper::Model::State::Random();









for (unsigned int ii = 0; ii < num_meas; ++ii) {


    m1 = source1.GenerateRandomMeasurement(track.state_, meas_std1);
    m2 = source2.GenerateRandomMeasurement(track.state_, meas_std2);
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




    for (auto& meas: this->sys.new_meas_) {


        err = TypeParam::ModelHelper::GetError(this->track,meas,this->sys.sources_);

        if(err.norm() <= this->sys.sources_[meas.source_index].params_.gate_threshold_) {
            in_validation_region = true;
        } else {
            in_validation_region = false;
            // std::cout << "err: " << std::endl <<  err << std::endl;
            // std::cout << "err norm: "  <<  err.norm() << std::endl;
            // std::cout << "meas: " << std::endl <<  meas.pose << std::endl;
            // std::cout << "track: " << std::endl <<  this->track.state_.g_.data_ << std::endl;
        }

        ASSERT_EQ(in_validation_region, validation_region.PolicyInValidationRegion(this->sys,meas,this->track));


    }
}

//-------------------------------------------------------------------------------------------------------------------------------

TYPED_TEST(ValidationRegionTest,ValidationRegionFixedPosPolicy ) {


    ValidationRegionFixedPosPolicy<typename TypeParam::ModelHelper::Model> validation_region;


    bool in_validation_region = false;
    Eigen::MatrixXd err;




    for (auto& meas: this->sys.new_meas_) {


        err = TypeParam::ModelHelper::GetPosError(this->track,meas);

        if(err.norm() <= this->sys.sources_[meas.source_index].params_.gate_threshold_) {
            in_validation_region = true;
        } else {
            in_validation_region = false;
            // std::cout << "err: " << std::endl <<  err << std::endl;
            // std::cout << "err norm: "  <<  err.norm() << std::endl;
            // std::cout << "meas: " << std::endl <<  meas.pose << std::endl;
            // std::cout << "track: " << std::endl <<  this->track.state_.g_.data_ << std::endl;
        }

        ASSERT_EQ(in_validation_region, validation_region.PolicyInValidationRegion(this->sys,meas,this->track));


    }
}

//-------------------------------------------------------------------------------------------------------------------------------

TYPED_TEST(ValidationRegionTest,ValidationRegionInnovPolicy ) {


    ValidationRegionInnovPolicy<typename TypeParam::ModelHelper::Model> validation_region;


    bool in_validation_region = false;
    Eigen::MatrixXd err;
    Eigen::MatrixXd S;




    for (auto& meas: this->sys.new_meas_) {


        err = TypeParam::ModelHelper::GetError(this->track,meas, this->sys.sources_);
        S = this->track.GetInnovationCovariance(this->sys.sources_, meas.source_index);

        // std::cout << "err: " << std::endl << err << std::endl;
        // std::cout << "S: " << std::endl << S << std::endl;

        if((err.transpose()*S.inverse()*err)(0,0) <= this->sys.sources_[meas.source_index].params_.gate_threshold_) {
            in_validation_region = true;
        } else {
            in_validation_region = false;
            // std::cout << "err: " << std::endl <<  err << std::endl;
            // std::cout << "err norm: "  <<  err.norm() << std::endl;
            // std::cout << "meas: " << std::endl <<  meas.pose << std::endl;
            // std::cout << "track: " << std::endl <<  this->track.state_.g_.data_ << std::endl;
        }

        ASSERT_EQ(in_validation_region, validation_region.PolicyInValidationRegion(this->sys,meas,this->track));


    }
}





