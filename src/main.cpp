#include <iostream>
#include <vector>
#include "common/measurement/measurement_base.h"
#include "common/sources/source_base.h"
#include <typeinfo>
#include "state.h"
#include <unsupported/Eigen/MatrixFunctions>

#include <chrono> 

   
// using namespace std;
 
// class IVectorParam
// {
// public:
//     virtual ~IVectorParam() {}
//     float y;
// };

// template <class T>
// class CVectorParam : public IVectorParam
// {
//     public:
//     std::vector<T> m_vect;
//     ~CVectorParam()=default;
// };
 
//  enum class Blah { A, B, C};

//  template<Blah type>
//  void test(){
//      std::cout << "hi" << std::endl;
//  }

//  template<>
//  void test<Blah::A>(){
//      std::cout << "A" << std::endl;
//  }

using namespace rransac;

enum Color {Red, Blue, Green};



Meas Func1() {
    Meas&& m = Meas();
    m.pose = Eigen::Matrix2d::Identity();
    m.twist = Eigen::Matrix2d::Identity();
    return m;

}

Meas Func2() {
    Meas m;
    m.pose = Eigen::Matrix2d::Identity();
    m.twist = Eigen::Matrix2d::Identity();
    return m;

}

int main(){

Eigen::Matrix4d tmp;
tmp.setIdentity();
tmp = tmp*0.1;
std::cout << tmp.inverse() << std::endl << std::endl; 
std::cout << tmp.sqrt() << std::endl << std::endl; 

return 0;
 
 

} 

