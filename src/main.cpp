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
#include "common/sources/source_SEN_pos_vel.h"
#include "common/sources/source_SEN_pose_twist.h"
#include "parameters.h"
#include "common/transformations/transformation_base.h"
#include "common/transformations/trans_homography.h"
#include "common/transformations/transformation_null.h"


using namespace rransac;
using namespace lie_groups;
int main(){



typedef ModelBase<SourceR2, TransformNULL<R2_r2>, 4, ModelRN<R2_r2, TransformNULL<R2_r2>>> Model1;

R2_r2 state;

Model1 model;

// model.GetLinTransFuncMatState(state,0.1);
// // boost::math::inverse_chi_squared dist_inv(3,1);
// // boost::math::chi_squared dist(3,1);
// boost::math::chi_squared dist(2);
// float test = boost::math::tgamma(1.0/2 +1);
// std::cout << pow(M_PI, 1.0/2.0)/test << std::endl;

// // std::cout << boost::math::quantile(dist,0.199) << std::endl;
// std::cout << boost::math::quantile(dist, 0.9997) << std::endl;
// // std::cout << boost::math::quantile(dist_inv, 0.198) << std::endl;

// SourceBase<lie_groups::R2_r2, SourceRN<lie_groups::R2_r2>> source;

return 0;
 
 

} 

