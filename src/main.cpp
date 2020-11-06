#include <iostream>
#include <vector>
#include "common/measurement/measurement_base.h"
#include "common/sources/source_base.h"
#include <typeinfo>
#include "state.h"

   
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

enum Color {Red, Blue, Green};

int main(){
 
// const Blah tmp = Blah::A;

// rransac::Meas m1, m2;

// m1.data = 5;
// m2.data = Eigen::Matrix2d::Identity();
// std::cout << m1.data << std::endl;
// std::cout << m2.data << std::endl;

rransac::SourceBase<lie_groups::R2_r2> source();

// rransac::Meas test
// std::cout << Red << std::endl;

// std::cout << rransac::MeasurementTypes::R2_POSE << std::endl;
return 0;
 
 

} 