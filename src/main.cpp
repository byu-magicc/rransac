
#include "common/models/model_RN.h"
#include "common/sources/source_RN.h"
#include "common/sources/source_SEN_pose_twist.h"
#include "state.h"
#include "common/transformations/transformation_null.h"

// struct Blah1{};
// struct Blah2{};

// template<typename T>
// class Test{
//     typedef double ModelCompatibility;
// }







int main() {


    rransac::ModelRN<lie_groups::R2_r2,rransac::SourceRN,rransac::TransformNULL> model;


    return 0;
}