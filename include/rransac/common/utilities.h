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


// Dynamic matrix
template<typename _DataType>
using MatXT = Eigen::Matrix<_DataType,Eigen::Dynamic,Eigen::Dynamic>;


//  Structs to determine source compatibility

struct CompatibleWithModelRN{};
struct CompatibleWithModelSENPosVel{};
struct CompatibleWithModelSENPoseTwist{};
struct CompatibleWithModelNull{};


inline unsigned int factorial(unsigned int n)
{
    if (n == 0)
        return 1;
    return n * factorial(n - 1);
}

template<typename T>
inline static void WrapAngle(T& angle) {
    while (angle < -M_PI) {
        angle += 2.0*M_PI;
    }
    while (angle > M_PI) {
        angle -= 2.0*M_PI;
    }


}



} // namespace utilities
} // namespace rransac

#endif //RRANSAC_COMMON_UTILITIES_H_