#ifndef RRANSAC_COMMON_UTILITIES_H_
#define RRANSAC_COMMON_UTILITIES_H_


#include <random>
#include <Eigen/Core>

namespace rransac {
    namespace utilities {


Eigen::MatrixXd GaussianRandomGenerator(const int size){

    std::default_random_engine gen(std::random_device{}());

    std::normal_distribution<double> dist_(0,1);
    Eigen::MatrixXd randn_nums(size,1);
    for (unsigned int ii = 0; ii < size; ++ii){
        randn_nums(ii,0) = dist_(gen);
    }

    return randn_nums;
}

    } // namespace utilities
} // namespace rransac



#endif //RRANSAC_COMMON_UTILITIES_H_