#ifndef RRANSAC_COMMON_UTILITIES_H_
#define RRANSAC_COMMON_UTILITIES_H_
#pragma once


#include <random>
#include <Eigen/Core>

namespace rransac {
namespace utilities {

template <typename T>
struct AlwaysFalse : std::false_type
{
};


inline Eigen::MatrixXd GaussianRandomGenerator(const int size){

    std::default_random_engine gen(std::random_device{}());

    std::normal_distribution<double> dist_(0,1);
    Eigen::MatrixXd randn_nums(size,1);
    for (unsigned int ii = 0; ii < size; ++ii){
        randn_nums(ii,0) = dist_(gen);
    }

    return randn_nums;
}

   


//  Structs to determine source compatibility

struct CompatibleWithModelRN{};
struct CompatibleWithModelSENPosVel{};
struct CompatibleWithModelSENPoseTwist{};
struct CompatibleWithModelNull{};



} // namespace utilities
} // namespace rransac

#endif //RRANSAC_COMMON_UTILITIES_H_