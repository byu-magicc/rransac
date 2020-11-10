#include <iostream>
#include <vector>
#include "common/measurement/measurement_base.h"
#include "common/sources/source_base.h"
#include <typeinfo>
#include "state.h"

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
 
 unsigned long num_iter = 100000;
auto start = std::chrono::high_resolution_clock::now();

// Eigen::MatrixXd dynamic1 = Eigen::Matrix4d::Identity();
// Eigen::MatrixXd dynamic2 = Eigen::Matrix4d::Identity();

for (unsigned long i=0; i < num_iter; ++i){
   Meas m =  Func2();
}

auto finish = std::chrono::high_resolution_clock::now();

std::chrono::duration<double> elapsed = finish - start;


std::cout << "Elapsed time: " << elapsed.count() << " s\n";

start = std::chrono::high_resolution_clock::now();

Eigen::Matrix4d static1 = Eigen::Matrix4d::Identity();
Eigen::Matrix4d static2 = Eigen::Matrix4d::Identity();

for (unsigned long i=0; i < num_iter; ++i){
     Meas m =  Func1();
}




finish = std::chrono::high_resolution_clock::now();

elapsed = finish - start;

std::cout << "Elapsed time: " << elapsed.count() << " s\n";
return 0;
 
 

} 

