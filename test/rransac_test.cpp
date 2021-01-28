#include <string>
#include <gtest/gtest.h>
#include "rransac.h"
#include "common/models/model_RN.h"
#include "state.h"
#include "common/transformations/transformation_null.h"
#include "common/transformations/trans_homography.h"
#include "common/data_association/model_policies/model_pdf_policy.h"
#include "common/data_association/cluster_data_tree_policies/data_tree_cluster_association_policy.h"
#include "track_initialization/seed_policies/null_policy.h"
#include "track_initialization/lmle_policies/linear_lmle_policy.h"



using namespace rransac;
using namespace lie_groups;

template <typename State>
struct CallbackClass {

bool func(const State& state) {
    if (state.g_.data_.norm() < 5)
        return true;
    else 
        return false;
}
};

typedef ModelRN<R2_r2,TransformHomography> Model;
typedef RRANSACTemplateParameters<Model,ModelPDFPolicy,DataTreeClusterAssociationPolicy,NULLSeedPolicy,LinearLMLEPolicy> RANSACParams;

TEST(RRANSACTest, AddSource) {

RRANSAC<RANSACParams> rransac;

CallbackClass<RANSACParams::State_> call;

// Init sources
typename RANSACParams::Source_ source1, source2, source3, source4;
SourceParameters source_params1, source_params2, source_params3, source_params4;

source_params1.meas_cov_ = Eigen::Matrix2d::Identity();
source_params1.type_ = MeasurementTypes::RN_POS;
source_params2.meas_cov_ = Eigen::Matrix2d::Identity();
source_params2.type_ = MeasurementTypes::RN_POS;
source_params3.meas_cov_ = Eigen::Matrix4d::Identity();
source_params3.type_ = MeasurementTypes::RN_POS_VEL;
source_params4.meas_cov_ = Eigen::Matrix4d::Identity();
source_params4.type_ = MeasurementTypes::RN_POS_VEL;

// Test invalid source index
source_params1.source_index_ = -1;
ASSERT_ANY_THROW(rransac.AddSource(source_params1));
ASSERT_ANY_THROW(rransac.AddSource(source_params1,std::bind(&CallbackClass<RANSACParams::State_>::func, call, std::placeholders::_1)));
source_params1.source_index_ = 1;
ASSERT_ANY_THROW(rransac.AddSource(source_params1));
ASSERT_ANY_THROW(rransac.AddSource(source_params1,std::bind(&CallbackClass<RANSACParams::State_>::func, call, std::placeholders::_1));

// Add valid source index
source_params1 = 0;
ASSERT_NO_THROW(rransac.AddSource(source_params1));

// Test invalid source index
ASSERT_ANY_THROW(rransac.AddSource(source_params1,std::bind(&CallbackClass<RANSACParams::State_>::func, call, std::placeholders::_1));

// Test valid source index
source_params2.source_index_ = 1;
ASSERT_ANY_THROW(rransac.AddSource(source_params2,std::bind(&CallbackClass<RANSACParams::State_>::func, call, std::placeholders::_1));

// Test invalid source index
source_params3.source_index_ = 0;
ASSERT_ANY_THROW(rransac.AddSource(source_params3));

// Test valid source index
source_params3.source_index_ = 2;
ASSERT_NO_THROW(rransac.AddSource(source_params3));

// Test invalid source index
source_params4.source_index_ = 4;
ASSERT_ANY_THROW(rransac.AddSource(source_params4));

// Test valid source index
source_params4.source_index_ = 3;
ASSERT_NO_THROW(rransac.AddSource(source_params4));

}







