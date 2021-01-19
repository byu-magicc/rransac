
// #include <typeinfo>
// #include <Eigen/Core>
// #include <unsupported/Eigen/MatrixFunctions>
#include <iostream>
// #include <chrono>
// #include "system.h"
// #include "state.h"
// #include "common/models/model_RN.h"
// #include "common/sources/source_RN.h"
// #include "common/transformations/transformation_null.h"
// #include "common/measurement/measurement_base.h"


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

template<typename T>
struct Blah {
    static void print(){
        std::cout << "here" << std::endl;
    }
};

template<typename T>
struct Sup {
    static void print(){
        std::cout << "here" << std::endl;
    }
};

template<typename T>
struct Meh {
    void print() { T::print();}
};

template<>
template<typename Q>
void Meh<Blah<Q>>::print() {Blah<T>::print(); Blah<T>::print();}

// template<typename T>
// struct Meh<Blah<T>>{
//     void print() {Blah<T>::print(); Blah<T>::print();}
// };


// using namespace rransac;
// using namespace lie_groups;
int main(){

Meh<Sup<int>> S;
Meh<Blah<double>> B;

S.print();
B.print();

return 0;

} 

