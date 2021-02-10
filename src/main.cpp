
#include "lie_groups/state.h"
#include<Eigen/Dense>
#include <iostream>


using namespace lie_groups;

int main() {

double dt = 1e-7;

  SE2<double> x;
  x = SE2<double>::Random();
  x.t_*=10;

  Eigen::Matrix<double,3,1> t;
  Eigen::Matrix<double,3,1> y;
  y.setRandom();
  y*=10;
  t.setRandom();
  t*=10;
  se2<double> l_g(se2<double>::Log(x.data_));
  se2<double> l_gh(se2<double>::Log(x.data_ *se2<double>::Exp(t*dt)));
  

  Eigen::Matrix<double,3,1> ans = (l_gh.JlInv()*l_g.Jl()*l_gh.JrInv()*y - l_g.JrInv()*y)/dt;

  std::cout << "ans: " << std::endl << y << std::endl;

  std::cout << "t: " << std::endl << se2<double>::Wedge(t)*y  << std::endl;




 return 0;   
}