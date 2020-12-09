
#include <typeinfo>
#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>
#include <iostream>
#include <chrono>
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

Eigen::MatrixXd md1, md2, md3;
Eigen::Matrix<double,10,10> ms1, ms2, ms3;
ms1.setRandom();
ms2.setRandom();
md1 = ms1;
md2 = ms2;
unsigned int num_iters = 100000;
auto start = std::chrono::high_resolution_clock::now();

for (unsigned int ii = 0; ii < num_iters; ++ii) {
    md3 = md1*md2;
}
 
auto finish = std::chrono::high_resolution_clock::now();
std::chrono::duration<double> elapsed = finish - start;
std::cout << "dynamic multiply: " << elapsed.count() << std::endl;

//-----------------------------------------------------------------------------

start = std::chrono::high_resolution_clock::now();
for (unsigned int ii = 0; ii < num_iters; ++ii) {
    ms3 = ms1*ms2;
}

finish = std::chrono::high_resolution_clock::now();
elapsed = finish - start;
std::cout << "static multiply: " << elapsed.count() << std::endl;

//-----------------------------------------------------------------------------

start = std::chrono::high_resolution_clock::now();
for (unsigned int ii = 0; ii < num_iters; ++ii) {
    md1 = md2;
}

finish = std::chrono::high_resolution_clock::now();
elapsed = finish - start;

std::cout << "dynamic assign: " << elapsed.count() << std::endl;

//-----------------------------------------------------------------------------

start = std::chrono::high_resolution_clock::now();
for (unsigned int ii = 0; ii < num_iters; ++ii) {
    ms1 = ms2;
}

finish = std::chrono::high_resolution_clock::now();
elapsed = finish - start;

std::cout << "static assign: " << elapsed.count() << std::endl;

//-----------------------------------------------------------------------------


start = std::chrono::high_resolution_clock::now();
for (unsigned int ii = 0; ii < num_iters; ++ii) {
    md3 = ms1*md2;
    ms3 = ms1*md2;

}

finish = std::chrono::high_resolution_clock::now();
elapsed = finish - start;

std::cout << "static dynamic multiply: " << elapsed.count()/2.0 << std::endl;

//-----------------------------------------------------------------------------


start = std::chrono::high_resolution_clock::now();
for (unsigned int ii = 0; ii < num_iters; ++ii) {
    md3 = md2;
    ms3 = ms1;

}

finish = std::chrono::high_resolution_clock::now();
elapsed = finish - start;

std::cout << "static dynamic assign: " << elapsed.count()/2.0 << std::endl;

return 0;

} 

