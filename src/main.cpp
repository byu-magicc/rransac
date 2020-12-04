
#include <typeinfo>
#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>
#include <iostream>


// #include "common/models/model_base.h"
// #include "common/models/model_RN.h"
// #include "common/models/model_SEN_pos_vel.h"
// #include "common/models/model_SEN_pose_twist.h"
// #include "common/sources/source_base.h"
// #include "common/sources/source_RN.h"
// #include "common/sources/source_SEN_pos_vel.h"
// #include "common/sources/source_SEN_pose_twist.h"
// #include "parameters.h"
// #include "common/transformations/transformation_base.h"
// #include "common/transformations/trans_homography.h"
// #include "common/transformations/transformation_null.h"


// using namespace rransac;
// using namespace lie_groups;
int main(){


Eigen::Matrix2d blah;
blah.setIdentity();

std::cout << blah.sqrt() << std::endl;
Eigen::Matrix2d blah_Sqrt = blah.sqrt();
return 0;
 
 

} 

