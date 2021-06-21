
#include "lie_groups/state.h"
#include <Eigen/Dense>
#include "rransac/common/models/model_SEN_pos_vel.h"
#include "rransac/common/sources/source_SEN_pos_vel.h"
#include "rransac/common/transformations/transformation_null.h"
#include <iostream>


using namespace lie_groups;

int main() {


  rransac::ModelSENPosVel<lie_groups::SE3_se3,rransac::TransformNULL,rransac::SourceSENPosVel> model;

  model.state_ = model.GetRandomState();

  Eigen::MatrixXd F = model.GetLinTransFuncMatState(model.state_,0.1);

  std::cout << "F: " << std::endl << F << std::endl;




 return 0;   
}