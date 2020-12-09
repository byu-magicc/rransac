
#include <typeinfo>
#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>
#include <iostream>
#include "system.h"
#include "state.h"
#include "common/models/model_RN.h"
#include "common/sources/source_RN.h"
#include "common/transformations/transformation_null.h"
#include "common/measurement/measurement_base.h"


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


using namespace rransac;
using namespace lie_groups;
int main(){

typedef ModelRN<R2_r2, TransformNULL> Model;

System<Model> sys;

Meas m;
m.time_stamp = 0;
m.pose = Eigen::Matrix2d::Identity();
sys.data_tree_.AddMeasurement(sys,m);

std::cout << sys.data_tree_.Size() << std::endl;

return 0;
 
 

} 

