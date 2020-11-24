#include <iostream>
#include <vector>
#include "common/measurement/measurement_base.h"
#include "common/sources/source_base.h"
#include "common/sources/source_RN.h"
#include <typeinfo>
#include "state.h"
#include <unsupported/Eigen/MatrixFunctions>

#include <chrono> 
#include <boost/math/distributions/inverse_chi_squared.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/special_functions/gamma.hpp>
using namespace rransac;
int main(){

// // boost::math::inverse_chi_squared dist_inv(3,1);
// // boost::math::chi_squared dist(3,1);
// boost::math::chi_squared dist(2);
// float test = boost::math::tgamma(1.0/2 +1);
// std::cout << pow(M_PI, 1.0/2.0)/test << std::endl;

// // std::cout << boost::math::quantile(dist,0.199) << std::endl;
// std::cout << boost::math::quantile(dist, 0.9997) << std::endl;
// // std::cout << boost::math::quantile(dist_inv, 0.198) << std::endl;

SourceBase<lie_groups::R2_r2, SourceRN<lie_groups::R2_r2>> source;

return 0;
 
 

} 

