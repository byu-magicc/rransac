#include <gtest/gtest.h>
#include "common/models/model_base.h"
#include "common/sources/source_base.h"
#include "common/sources/source_RN.h"
#include "common/sources/source_SEN_pos_vel.h"
#include "common/sources/source_SEN_pose_twist.h"
#include "parameters.h"
#include "common/transformations/transformation_base.h"
#include "common/transformations/trans_homography.h"


namespace rransac
{

using namespace lie_groups;

typedef ModelBase<R2_r2, R2_r2::StateType, SourceRN<R2_r2>, TransformNULL<R2_r2>> Model2;
typedef ModelBase<R3_r3, R3_r3::StateType, SourceRN<R3_r3>, TransformNULL<R3_r3>> Model2;
typedef ModelBase<SE2_se2, SE2_se2::StateType, SourceSENPosVel<SE2_se2>, TransformNULL<SE2_se2>> Model3;
typedef ModelBase<SE2_se2, SE2_se2::StateType, SourceSENPoseTwist<SE2_se2>, TransformNULL<SE2_se2>> Model3;
typedef ModelBase<SE3_Se3, SE3_Se3::StateType, SourceSENPosVel<SE3_Se3>, TransformNULL<SE3_Se3>> Model3;
typedef ModelBase<SE3_Se3, SE3_Se3::StateType, SourceSENPoseTwist<SE3_Se3>, TransformNULL<SE3_Se3>> Model3;

template<class S>
class ModelTest : public ::testing::Test {

public:

SourceParameters source_params1;
SourceParameters source_params2;




};

using MyTypes = ::testing::Types< Model1, Model2, Model3, Model4 >;
TYPED_TEST_SUITE(ModelTest, MyTypes);


template <typename State, typename StateType, typename Source, typename Transformation> 
class ModelBase

TYPED_TEST(ModelTest, SetData) {



    
} // namespace rransac
