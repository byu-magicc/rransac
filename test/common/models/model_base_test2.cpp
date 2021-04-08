#include <gtest/gtest.h>
#include <typeinfo>
#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>


#include "rransac/common/models/model_base.h"
#include "rransac/common/models/model_RN.h"
#include "rransac/common/models/model_SEN_pos_vel.h"
#include "rransac/common/models/model_SEN_pose_twist.h"
#include "rransac/common/sources/source_base.h"
#include "rransac/common/sources/source_RN.h"
#include "rransac/parameters.h"
#include "rransac/common/transformations/transformation_base.h"
#include "rransac/common/transformations/trans_homography.h"
#include "rransac/common/transformations/transformation_null.h"

namespace rransac
{

using namespace lie_groups;


typedef ModelRN<R3_r3,TransformNULL,SourceRN> Model2;




// ----------------------------------------------------------------------------------

TEST(ModelBaseTest2_, OMinus){

ModelSENPosVel<SE2_se2,TransformNULL,SourceSENPosVel> model1, model2;

typename SE2_se2::Mat_SC cartesian;
cartesian << 1,2,0.1,3,4,0.2;

model1.state_ = model1.state_.Random();
model2.state_ = model1.state_.OPlus(cartesian);

Eigen::MatrixXd diff = model1.OMinus(model2,model1);

ASSERT_DOUBLE_EQ(diff(0,0), cartesian(0,0));
ASSERT_DOUBLE_EQ(diff(1,0), cartesian(1,0));
ASSERT_DOUBLE_EQ(diff(2,0), cartesian(2,0));
ASSERT_DOUBLE_EQ(diff(3,0), cartesian(3,0));
ASSERT_DOUBLE_EQ(diff(4,0), cartesian(5,0));

ModelSENPoseTwist<SE2_se2,TransformNULL,SourceSENPoseTwist> model3, model4;
model3.state_ = model3.state_.Random();
model4.state_ = model3.state_.OPlus(cartesian);

Eigen::MatrixXd diff2 = model3.OMinus(model4,model3);
ASSERT_LE( ( diff2 - cartesian).norm(), 1e-10);

ModelRN<R3_r3,TransformNULL,SourceRN> model5, model6;
model5.state_ = model5.state_.Random();
model6.state_ = model5.state_.OPlus(cartesian);

Eigen::MatrixXd diff3 = model5.OMinus(model6,model5);
ASSERT_LE( ( diff3 - cartesian).norm(), 1e-10);


}



} // namespace rransac
