#include <gtest/gtest.h>
#include <typeinfo>
#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>


#include "common/models/model_base.h"
#include "common/models/model_RN.h"
#include "common/models/model_SEN_pos_vel.h"
#include "common/models/model_SEN_pose_twist.h"
#include "common/sources/source_base.h"
#include "common/sources/source_RN.h"
#include "parameters.h"
#include "common/transformations/transformation_base.h"
#include "common/transformations/trans_homography.h"
#include "common/transformations/transformation_null.h"

namespace rransac
{


TEST(ModelBaseTest2, INIT) {

}

TEST(ModelBaseTest2, UpdateLikelihood) {
    
}

TEST(ModelBaseTest2, ConsensusSet) {
    
}

} // namespace rransac
