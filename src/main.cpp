
#include "lie_groups/state.h"
#include <Eigen/Dense>
// #include "rransac/common/models/model_SEN_pos_vel.h"
// #include "rransac/common/sources/source_SEN_pos_vel.h"
// #include "rransac/common/transformations/transformation_null.h"
#include <iostream>
#include <tuple>
#include <variant>

using namespace lie_groups;


template <typename tDerived>
class Animals {

public: 
int Bark(){
  return static_cast<tDerived*>(this)->DerivedBark();
}

};

class Dog : public Animals<Dog> {

public:
int a = 1;
int DerivedBark() {
  return a;
}


};


class Cat : public Animals<Cat> {

public:
int a = 2;
int DerivedBark() {
  return a;
}

};

class Dummy : public Animals<Cat> {

public:
int a = 0;
int DerivedBark() {
  return a;
}

};


template<typename... Ts>
struct Count;

template<typename Animal, typename... Ts>
struct Count<Animal, Ts...>
{


static constexpr int num = std::is_same<Animal,Dummy>::value ? 0 :1 + Count<Ts...>::num;

};

template <>
struct Count<> {
  static constexpr int num = 0;
};




template<typename T1=Dummy, typename T2=Dummy, typename T3=Dummy >
class AnimalContainter {
public: 
std::tuple<T1,T2,T3> animals;

static constexpr int num = Count<T1,T2,T3>::num;

int Bark(int index){


  switch (index)
  {
  case 0:
    return std::get<0>(animals).Bark();
    break;
  case 1:
    return std::get<1>(animals).Bark();
    break;
  
  default:
    return std::get<2>(animals).Bark();
    break;
  }
  
}

};



int main() {


AnimalContainter<Cat,Cat,Dog> animals;




long int num_sim = 1e3;





// std::get<1>(animals.animals).a = 4;

std::cout << "bark 1 " << animals.Bark(0) << std::endl;
std::cout << "bark 2 " << animals.Bark(1) << std::endl;
std::cout << "num " << animals.num << std::endl;





 return 0;   
}