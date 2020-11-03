#ifndef RRANSAC_COMMON_SOURCES_SOURCE_BASE_H_
#define RRANSAC_COMMON_SOURCES_SOURCE_BASE_H_

#include <Eigen/Core>
#include "state.h"

namespace rransac
{

// Lists the different types of sources available
enum class SourceTypes {
    R2_POS,                // The Target's manifold is R2 and only position is observed
    R2_POS_VEL,            // The Target's manifold is R2 and position and velocity is observed
    SE2_POS,               // The Target's manifold is SE2 and only position is observed
    SE2_POS_ATT            // The target's manifold is SE2 and the position and attitude is observed
};

/** \class SourceParameters
 * This struct contains the parameters needed by a measurement source.
 */ 
struct SourceParameters {
    bool meas_cov_fixed_;            /** < Flag used to indicate if the measurement covariance is the same for every measurement. 
                                      If it is, then the measurement covariance needs to be given to the source. If it isn't, 
                                      then the measurement covariance needs to be given with the measurement. @see Meas */
    Eigen::MatrixXd meas_cov_;       /** < The fixed measurement covariance.*/
    float expected_num_false_meas_;  /** < The expected number of false measurements. We assume that the expected number of false 
                                        measurements from a source per sensor scan can be modeled using a Poisson distribution. */
    unsigned int source_id_;    /** < An identifier that is unique for every source. When a measurement is given, the measurement
                                      indicates which source it came from using the source ID. @see Meas */

    SourceTypes type_;          /** < The source type @see SourceTypes */

};

/** \class SourceBase 
 * A source is an algorithm that takes data from a sensor and extracts measurements from it. Due to the nature of
 * the source, there is measurement noise and false measurements. The measurement noise is assumed to be sampled from
 * a white-noise, zero-mean, Gaussian distribution. If the measurement covariance is fixed, then it should be added to 
 * the source; otherwise, it is supplied with each measurement. We assume that the false measurements are uniformly
 * distributed in the local surveillance region, and that the spatial density of false measurements can 
 * be modeled using a Poisson distribution. 
 * 
 * Since there can be many different sources, we allow users to create their own source type using polymorphism. The 
 * base source must contain all of the necessary member functions and the child classes provide the implementation for
 * specific sources. 
 * 
 * When a measurement is received, the measurement indicates which source it came from using the source ID tag @see Meas.
 * We also need to know the dimension of the measurement.
 * 
 */ 


class SourceBase
{

public:

    SourceParameters params_;  /** < The source parameters @see SourceParameters */

    // unsigned int meas_dim_;     /** < The dimension of the measurement. You can also consider this as the number of generalized coordinates measured.
    //                               For example, if the measurement containted the position of a 2D object in the
    //                               xy-plane, then the dimension of the measurement would be two. This is used by RANSAC when creating model hypotheses. */

    template<SourceTypes type>
    inline void Init(const SourceParameters& params); /** Initializes the measurement source */
    
    template <SourceTypes type, class S>
    inline Eigen::MatrixXd GetLinObsMatState(const S& state);                               /** Returns the jacobian of the observation function w.r.t. the states */
    
    template <SourceTypes type, class S>
    inline Eigen::MatrixXd GetLinObsMatSensorNoise(const S& state);                         /** Returns the jacobian of the observation function w.r.t. the sensor noise */

    template <SourceTypes type, class S>
    inline Eigen::MatrixXd GetEstMeas(const S& state); /** Returns an estimated measurement according to the state. */
    
private:
    Eigen::MatrixXd H_;
    Eigen::MatrixXd V_;

};


///////////////////////////////////////////////////////////////////////
//                             R2_POS
///////////////////////////////////////////////////////////////////////
template<>
inline void SourceBase::Init<SourceTypes::R2_POS>(const SourceParameters& params) {
    params_ = params;
    Eigen::Matrix<double,2,4> H;
    Eigen::Matrix<double, 2,2> V;
    H.block(0,0,2,2).setIdentity();
    H.block(0,2,2,2).setZero();
    V.setIdentity();

    H_ = H;
    V_ = V;
}

//-------------------------------------------

template <>
inline Eigen::MatrixXd SourceBase::GetLinObsMatState<SourceTypes::R2_POS,lie_groups::R2_r2>(const lie_groups::R2_r2& state) {
    return H_;
}                             

//-------------------------------------------

template <>
inline Eigen::MatrixXd SourceBase::GetLinObsMatSensorNoise<SourceTypes::R2_POS,lie_groups::R2_r2>(const lie_groups::R2_r2& state) {
    return V_;
}

//-------------------------------------------

template <>
inline Eigen::MatrixXd SourceBase::GetEstMeas<SourceTypes::R2_POS,lie_groups::R2_r2>(const lie_groups::R2_r2& state) {
    return state.g_.data_;
}

///////////////////////////////////////////////////////////////////////
//                             R2_POS_VEL
///////////////////////////////////////////////////////////////////////
template<>
inline void SourceBase::Init<SourceTypes::R2_POS_VEL>(const SourceParameters& params) {
    params_ = params;
    Eigen::Matrix<double,4,4> H;
    Eigen::Matrix<double, 4,4> V;
    H.setIdentity();
    V.setIdentity();

    H_ = H;
    V_ = V;
}

//-------------------------------------------

template <>
inline Eigen::MatrixXd SourceBase::GetLinObsMatState<SourceTypes::R2_POS_VEL>(const lie_groups::R2_r2& state) {
    return H_;
}                             

//-------------------------------------------

template <>
inline Eigen::MatrixXd SourceBase::GetLinObsMatSensorNoise<SourceTypes::R2_POS_VEL>(const lie_groups::R2_r2& state) {
    return V_;
}

//-------------------------------------------

template <>
inline Eigen::MatrixXd SourceBase::GetEstMeas<SourceTypes::R2_POS_VEL>(const lie_groups::R2_r2& state) {
    Eigen::Matrix<double,4,1> tmp;
    tmp.block(0,0,2,1) = state.g_.data_;
    tmp.block(2,0,2,1) = state.u_.data_;
    return tmp;
}

} // namespace rransac

#endif // RRANSAC_COMMON_SOURCES_SOURCE_BASE_H_