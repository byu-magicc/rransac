#ifndef RRANSAC_COMMON_MODELS__MODEL_RN_JACOBIAN_H_
#define RRANSAC_COMMON_MODELS__MODEL_RN_JACOBIAN_H_
#pragma once



#include "lie_groups/state.h"

namespace rransac {


/**
 * \class ConstructF
 * This class is used to construct the F Jacobian for RN with some cases be optimal.
 */ 
template<typename _DataType, int _GroupDim, int _NumTangentSpaces>
class ConstructF {

typedef _DataType DataType;
static constexpr unsigned int u_dim_ = _GroupDim*_NumTangentSpaces;             /**< The dimension of the total tangent space. */
static constexpr unsigned int cov_dim_ = _GroupDim+u_dim_;                      /**< The dimension of the error covariance. */
typedef Eigen::Matrix<_DataType,cov_dim_,cov_dim_> Mat;                         /**< The object type of the error covariance, Jacobians, and others. */

public:

    ConstructF()=default;

   static  Mat GetF(DataType dt) {

        Mat A0(Mat::Zero());
        A0.block(0,_GroupDim,u_dim_,u_dim_).setIdentity();
        Mat Ac = A0;
    
        Mat F = Mat::Identity();
        for (int ii = 0; ii < _NumTangentSpaces; ++ii) {
            F+= Ac*pow(dt,ii+1)/utilities::factorial(ii+1);
            Ac*=A0;
        }
        return F;
    }

private:

};

//////////////////////////////////////////////////////////
// Optimized for 2, 1
/////////////////////////////////////////////////////////

template<typename _DataType>
class ConstructF<_DataType,2,1> {
public:

static constexpr unsigned int g_dim_ = 2;
static constexpr unsigned int num_tangent_spaces_ = 1;
static constexpr unsigned int u_dim_ = g_dim_*num_tangent_spaces_;             /**< The dimension of the total tangent space. */
static constexpr unsigned int cov_dim_ = g_dim_+u_dim_;                      /**< The dimension of the error covariance. */
typedef Eigen::Matrix<_DataType,cov_dim_,cov_dim_> Mat;                         /**< The object type of the error covariance, Jacobians, and others. */




    ConstructF()=default;

   static  Mat GetF(_DataType dt) {
       _DataType one_(1.0);
       _DataType zero_(0.0);

       Mat F;
       F << one_, zero_, dt, zero_, zero_, one_, zero_, dt,zero_, zero_, one_, zero_,zero_, zero_, zero_, one_;
      return F;
    }

private:

};

//////////////////////////////////////////////////////////
// Optimized for 2, 2
/////////////////////////////////////////////////////////


template<typename _DataType>
class ConstructF<_DataType,2,2> {
static constexpr unsigned int g_dim_ = 2;
static constexpr unsigned int num_tangent_spaces_ = 2;
static constexpr unsigned int u_dim_ = g_dim_*num_tangent_spaces_;             /**< The dimension of the total tangent space. */
static constexpr unsigned int cov_dim_ = g_dim_+u_dim_;                      /**< The dimension of the error covariance. */
typedef Eigen::Matrix<_DataType,cov_dim_,cov_dim_> Mat;                         /**< The object type of the error covariance, Jacobians, and others. */



public:

    ConstructF()=default;

   static  Mat GetF(_DataType dt) {
       _DataType one_(1.0);
       _DataType zero_(0.0);
       _DataType dt2 = dt*dt/2.0;


       Mat F;
       F << one_, zero_, dt, zero_, dt2, zero_, zero_, one_, zero_, dt, zero_, dt2, zero_, zero_, one_, zero_, dt, zero_, zero_, zero_, zero_, one_, zero_, dt, zero_, zero_, zero_, zero_, one_, zero_, zero_, zero_, zero_, zero_, zero_, one_;
      return F;
    }

private:

};

//////////////////////////////////////////////////////////
// Optimized for 2, 3
/////////////////////////////////////////////////////////


template<typename _DataType>
class ConstructF<_DataType,2,3> {
static constexpr unsigned int g_dim_ = 2;
static constexpr unsigned int num_tangent_spaces_ = 3;
static constexpr unsigned int u_dim_ = g_dim_*num_tangent_spaces_;             /**< The dimension of the total tangent space. */
static constexpr unsigned int cov_dim_ = g_dim_+u_dim_;                      /**< The dimension of the error covariance. */
typedef Eigen::Matrix<_DataType,cov_dim_,cov_dim_> Mat;                         /**< The object type of the error covariance, Jacobians, and others. */


public:

    ConstructF()=default;

   static  Mat GetF(_DataType dt) {
       _DataType one_(1.0);
       _DataType zero_(0.0);
       _DataType dt2 = dt*dt/2.0;
       _DataType dt3 = dt2*dt/3.0;
       Mat F;
       F << one_, zero_, dt, zero_, dt2, zero_, dt3, zero_, zero_, one_, zero_, dt, zero_, dt2, zero_, dt3, zero_, zero_, one_, zero_, dt, zero_, dt2, zero_, zero_, zero_, zero_, one_, zero_, dt, zero_, dt2, zero_, zero_, zero_, zero_, one_, zero_, dt, zero_, zero_, zero_, zero_, zero_, zero_, one_, zero_, dt, zero_, zero_, zero_, zero_, zero_, zero_, one_, zero_, zero_, zero_, zero_, zero_, zero_, zero_, zero_, one_;
      return F;
    }

private:

};

//////////////////////////////////////////////////////////
// Optimized for 3, 1
/////////////////////////////////////////////////////////


template<typename _DataType>
class ConstructF<_DataType,3,1> {
static constexpr unsigned int g_dim_ = 3;
static constexpr unsigned int num_tangent_spaces_ = 1;
static constexpr unsigned int u_dim_ = g_dim_*num_tangent_spaces_;             /**< The dimension of the total tangent space. */
static constexpr unsigned int cov_dim_ = g_dim_+u_dim_;                      /**< The dimension of the error covariance. */
typedef Eigen::Matrix<_DataType,cov_dim_,cov_dim_> Mat;                         /**< The object type of the error covariance, Jacobians, and others. */


public:

    ConstructF()=default;

   static  Mat GetF(_DataType dt) {
       _DataType one_(1.0);
       _DataType zero_(0.0);

       Mat F;
       F << one_, zero_, zero_, dt, zero_, zero_, zero_, one_, zero_, zero_, dt, zero_, zero_, zero_, one_, zero_, zero_, dt, zero_, zero_, zero_, one_, zero_, zero_, zero_, zero_, zero_, zero_, one_, zero_, zero_, zero_, zero_, zero_, zero_, one_;
      return F;
    }

private:

};

//////////////////////////////////////////////////////////
// Optimized for 3, 2
/////////////////////////////////////////////////////////


template<typename _DataType>
class ConstructF<_DataType,3,2> {
static constexpr unsigned int g_dim_ = 3;
static constexpr unsigned int num_tangent_spaces_ = 2;
static constexpr unsigned int u_dim_ = g_dim_*num_tangent_spaces_;             /**< The dimension of the total tangent space. */
static constexpr unsigned int cov_dim_ = g_dim_+u_dim_;                      /**< The dimension of the error covariance. */
typedef Eigen::Matrix<_DataType,cov_dim_,cov_dim_> Mat;                         /**< The object type of the error covariance, Jacobians, and others. */


public:
     

    ConstructF()=default;

   static  Mat GetF(_DataType dt) {
       _DataType one_(1.0);
       _DataType zero_(0.0);
       _DataType dt2 = dt*dt/2.0;
       Mat F;
       F << one_, zero_, zero_, dt, zero_, zero_, dt2, zero_, zero_, zero_, one_, zero_, zero_, dt, zero_, zero_, dt2, zero_, zero_, zero_, one_, zero_, zero_, dt, zero_, zero_, dt2, zero_, zero_, zero_, one_, zero_, zero_, dt, zero_, zero_, zero_, zero_, zero_, zero_, one_, zero_, zero_, dt, zero_, zero_, zero_, zero_, zero_, zero_, one_, zero_, zero_, dt, zero_, zero_, zero_, zero_, zero_, zero_, one_, zero_, zero_, zero_, zero_, zero_, zero_, zero_, zero_, zero_, one_, zero_, zero_, zero_, zero_, zero_, zero_, zero_, zero_, zero_, one_;
      return F;
    }

private:

};

//////////////////////////////////////////////////////////
// Optimized for 3, 3
/////////////////////////////////////////////////////////


template<typename _DataType>
class ConstructF<_DataType,3,3> {
static constexpr unsigned int g_dim_ = 3;
static constexpr unsigned int num_tangent_spaces_ = 3;
static constexpr unsigned int u_dim_ = g_dim_*num_tangent_spaces_;             /**< The dimension of the total tangent space. */
static constexpr unsigned int cov_dim_ = g_dim_+u_dim_;                      /**< The dimension of the error covariance. */
typedef Eigen::Matrix<_DataType,cov_dim_,cov_dim_> Mat;                         /**< The object type of the error covariance, Jacobians, and others. */


public:

    ConstructF()=default;

   static  Mat GetF(_DataType dt) {
       _DataType one_(1.0);
       _DataType zero_(0.0);

       _DataType dt2 = dt*dt/2.0;
       _DataType dt3 = dt2*dt/3.0;
       Mat F;
       F << one_, zero_, zero_, dt, zero_, zero_, dt2, zero_, zero_, dt3, zero_, zero_, zero_, one_, zero_, zero_, dt, zero_, zero_, dt2, zero_, zero_, dt3, zero_, zero_, zero_, one_, zero_, zero_, dt, zero_, zero_, dt2, zero_, zero_, dt3, zero_, zero_, zero_, one_, zero_, zero_, dt, zero_, zero_, dt2, zero_, zero_, zero_, zero_, zero_, zero_, one_, zero_, zero_, dt, zero_, zero_, dt2, zero_, zero_, zero_, zero_, zero_, zero_, one_, zero_, zero_, dt, zero_, zero_, dt2, zero_, zero_, zero_, zero_, zero_, zero_, one_, zero_, zero_, dt, zero_, zero_, zero_, zero_, zero_, zero_, zero_, zero_, zero_, one_, zero_, zero_, dt, zero_, zero_, zero_, zero_, zero_, zero_, zero_, zero_, zero_, one_, zero_, zero_, dt, zero_, zero_, zero_, zero_, zero_, zero_, zero_, zero_, zero_, one_, zero_, zero_, zero_, zero_, zero_, zero_, zero_, zero_, zero_, zero_, zero_, zero_, one_, zero_, zero_, zero_, zero_, zero_, zero_, zero_, zero_, zero_, zero_, zero_, zero_, one_;
      return F;
    }

private:

};


} // rransac


#endif // RRANSAC_COMMON_MODELS__MODEL_RN_JACOBIAN_H_