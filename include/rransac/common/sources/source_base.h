#ifndef RRANSAC_COMMON_SOURCES_SOURCE_BASE_H_
#define RRANSAC_COMMON_SOURCES_SOURCE_BASE_H_
#pragma once


#include <Eigen/Core>
#include <cmath>
#include <typeinfo>
#include <random>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <functional>


#include "lie_groups/state.h"
#include "rransac/parameters.h"
#include "rransac/common/measurement/measurement_base.h"
#include "rransac/common/utilities.h"

namespace rransac
{

template <typename _State, typename _MatH, typename _MatV, template <typename> typename _Transformation, int _MeasPoseRows, int _MeasPoseCols, int _MeasTwistRows, int _MeasTwistCols, int _MeasPoseDim, int _MeasTwistDim, bool _HasVel, MeasurementTypes _MeasurementType, typename _ModelCompatibility>
struct SourceDerivedTraits {

typedef _State State;                                            /**< The state of the target. @see State. */
typedef typename _State::DataType DataType;                      /**< The scalar object for the data. Ex. float, double, etc. */
typedef _MatH MatH;                                              /**< The object type of the Jacobian H. */
typedef _MatV MatV;                                              /**< The object type of the Jacobians V. */
template<typename T>
using Transformation = _Transformation<T>;                       /**< The transformation used to transform the measurements and tracks. */
static constexpr unsigned int meas_pose_rows_  =_MeasPoseRows;   /**< The number of rows in the pose measurement. */
static constexpr unsigned int meas_pose_cols_  =_MeasPoseCols;   /**< The number of columns in the pose measurement. */
static constexpr unsigned int meas_twist_rows_ =_MeasTwistRows;  /**< The number of rows in the twist measurement. */
static constexpr unsigned int meas_twist_cols_ =_MeasTwistCols;  /**< The number of columns in the twist measurement. */
static constexpr unsigned int meas_pose_dim_   =_MeasPoseDim;    /**< The measurement pose dimension. */
static constexpr unsigned int meas_twist_dim_  =_MeasTwistDim;   /**< The measurement twist dimension. */
static constexpr unsigned int total_meas_dim_  =meas_pose_dim_+meas_twist_dim_;   /**< The total measurement dimension. */
static constexpr bool has_vel_ = _HasVel;                        /**< Indicates if the measurement contains velocity.  */
static constexpr MeasurementTypes measurement_type_ = _MeasurementType; /**< The measurement type of the source. */
typedef Eigen::Matrix<DataType,total_meas_dim_,total_meas_dim_> MatMeasCov; /**< The measurement covariance type. */
typedef Eigen::Matrix<DataType,total_meas_dim_,1> VecMeas;       /**< The error type of the difference between two measurements. */
typedef _ModelCompatibility ModelCompatibility;                  /**< Indicates which model the source type is compatible with. */


};

/** \class SourceParameters
 * This struct contains the parameters needed by a measurement source. Some of the parameters are defined by the user, while others are calculated by RRANSAC.
 * The parameters that the user must specify are SourceParameters::source_index_, SourceParameters::type_, SourceParameters::meas_cov_, SourceParameters::spacial_density_of_false_meas_,
 * SourceParameters::probability_of_detection_, SourceParameters::gate_probability_, and SourceParameters::RANSAC_inlier_probability_. Some of the parameters have specified conditions on them, 
 * make sure that you set them correctly. 
 */ 
struct SourceParameters {

    unsigned int source_index_;       /**< This is a unique identifiable index for a measurement source. The first source added to RRANSAC must have a source index value of 0, and 
                                           the second a value of 1, and so on. The source index is used to verify that a measurement corresponds to a proper source. This parameter must be set by the user. */
    MeasurementTypes type_;           /**< The measurement type. This parameter must be set by the user. @see MeasurementTypes */
                                      
    Eigen::MatrixXd meas_cov_;        /**< The measurement covariance. It must be positive definite. This parameter must be set by the user. */

    double spacial_density_of_false_meas_;  /**< We assume that measurement sources produce false measurements uniformly across their surveillance region. The number
                                                 of false measurements per sensor scan is modeled using a Poison distribution with expected value \f$ \lambda \f$. The spacial density
                                                 of false measurements is the expected value of the Poison distribution divided by the volume of the surveillance region. See Tracking and 
                                                 Data Fusion by Bar-Shalom 2011 for more detail. This parameter must be set by the user. */    

    double probability_of_detection_; /**< The probability that the target of interest is detected by the source during
                                          a sensor scan. This value must be between 0 and 1. This parameter must be set by the user. */

    double gate_probability_;         /**< The probability that a true measurement will be inside the validation region of a track. This value must
                                          be between 0 and 1. This parameter must be set by the user. See Tracking and Data Fusion by Bar-Shalom 2011 for more detail on the validation region. 
                                          It is also used during track initialization. 
                                          During track initialization, a hypothetical state is scored according to the number of measurement inliers. A measurement is an inlier if
                                          it is a certain distance from the hypothetical state. The threshold on the distance is determined by this parameter. Thus, it is 
                                          the probability that a measurement is an inlier according to a chi squared distribution. See Tracking and Data Fusion by Bar-Shalom 2011 for 
                                          more detail. This parameter must be set by the user. */

    double gate_threshold_;              /**< The gate threshold of the validation region. If this value is non positive, then the gate_threshold_ will be calculated using SourceParameters::gate_probability_ as 
                                              described in Tracking and Data Fusion by Bar-Shalom 2011 regarding the validation region. If this value is set to a positive value, then it will retain that 
                                              value. By default, this value is set to zero in the constructor. This gate threshold is used to associate measurements to tracks and to determine if a measurement
                                              is an inlier during the track initialization phase.*/
    

    // These parameters are not defined by the user, but are calculated depending on the user specified parameters.
    double gate_threshold_sqrt_;         /**< The square root of the gate threshold. */    
    double vol_unit_hypershpere_;        /**< The volume of the unit hypershpere in the measurement space. This is determined by the dimension of the measurement source. */
    unsigned int meas_pose_rows_;        /**< The number of rows in the pose measurement. */
    unsigned int meas_pose_cols_;        /**< The number of columns in the pose measurement. */
    unsigned int meas_twist_rows_;       /**< The number of rows in the twist measurement. */
    unsigned int meas_twist_cols_;       /**< The number of columns in the twist measurement. */
    unsigned int meas_pose_dim_;         /**< The measurement pose dimension. */
    unsigned int meas_twist_dim_;        /**< The measurement twist dimension. */
    unsigned int total_meas_dim_;        /**< The total measurement dimension. */
    bool has_vel_;                       /**< Indicates if the measurement contains velocity or twist data.  */


    // Sets some parameters to default
    SourceParameters() {
        spacial_density_of_false_meas_ = 0.1;
        probability_of_detection_ = 0.9;
        gate_probability_ = 0.8;
        gate_threshold_ = 0;
    }

};

/** \class SourceBase 
 * A source is an algorithm that takes data from a sensor and extracts measurements from it. Due to the nature of
 * the source, there is measurement noise and false measurements. The measurement noise is assumed to be sampled from
 * a white-noise, zero-mean, Gaussian distribution. We assume that the false measurements are uniformly
 * distributed in the local surveillance region, and that the spatial density of false measurements can 
 * be modeled using a Poisson distribution. 
 * 
 * The observation model for a specific source is 
 * \f[
 * y = h(x,v)
 * \f]
 * with \f$ y \f$ denoting the measurement, \f$ h \f$ the observation function, \f$ x \f$ the state of the target, and \f$ v \f$ the measurement noise. 
 * 
 * The linearization of the observation model about a state \f$ x_0 \f$ is
 * 
 * \f[
 * y = h(x_0,0) + H\delta x + Vv
 * \f]
 * with \f$ H \f$ being the partial derivative of \f$ h \f$ w.r.t. the state evaluated at \f$ x_0 \f$, \f$ V \f$ being the partial derivative of \f$ h \f$ w.r.t.
 * the noise evaluated at \f$ x_0 \f$, and \f$ \delta x \f$ being the geodesic distance between \f$ x \f$ and \f$ x_0 \f$.
 * 
 * A source must implement the observation fuction, calculate the Jacobians \f$ H \f$ and \f$ V \f$. Also, it must be able to compute the distance between two measurements. 
 * 
 * Since there can be many different sources, we allow users to create their own source type using the design methology curiously recurring template pattern. The 
 * base source acts as the API for all of the sources, and each child source must implement certain methods. 
 * 
 * When a measurement is received, the measurement indicates which source it came from using the source ID tag @see Meas.
 * 
 */ 




template<typename _SourceDerivedTraits, template<typename , MeasurementTypes , template <typename > typename > typename _Derived>
class SourceBase
{

public:

    typedef _SourceDerivedTraits DerivedTraits;                                             /**< The traits of the derived class. */
    typedef typename DerivedTraits::State State;                                            /**< The state of the target. @see State. */
    typedef typename DerivedTraits::DataType DataType;                                      /**< The scalar object for the data. Ex. float, double, etc. */
    typedef typename DerivedTraits::MatH MatH;                                              /**< The object type of the Jacobian H. */
    typedef typename DerivedTraits::MatV MatV;                                              /**< The object type of the Jacobians V. */
    typedef typename DerivedTraits::template Transformation<State> Transformation;          /**< The transformation used to transform the measurements and tracks. */
    typedef typename DerivedTraits::MatMeasCov MatMeasCov;                                  /**< The data type of the measurement covariance. */
    typedef typename DerivedTraits::VecMeas VecMeas;                                        /**< The data type of the measurement covariance. */
    static constexpr unsigned int meas_pose_rows_  = DerivedTraits::meas_pose_rows_;        /**< The number of rows in the pose measurement. */
    static constexpr unsigned int meas_pose_cols_  = DerivedTraits::meas_pose_cols_;        /**< The number of columns in the pose measurement. */
    static constexpr unsigned int meas_twist_rows_ = DerivedTraits::meas_twist_rows_;       /**< The number of rows in the twist measurement. */
    static constexpr unsigned int meas_twist_cols_ = DerivedTraits::meas_twist_cols_;       /**< The number of columns in the twist measurement. */
    static constexpr unsigned int meas_pose_dim_   = DerivedTraits::meas_pose_dim_;         /**< The measurement pose dimension. */
    static constexpr unsigned int meas_twist_dim_  = DerivedTraits::meas_twist_dim_;        /**< The measurement twist dimension. */
    static constexpr unsigned int total_meas_dim_  = DerivedTraits::total_meas_dim_;        /**< The total measurement dimension. */
    static constexpr bool has_vel_ = DerivedTraits::has_vel_;                               /**< Indicates if the measurement contains velocity.  */
    static constexpr MeasurementTypes measurement_type_ = DerivedTraits::measurement_type_; /**< The measurement type of the source. */
    typedef typename DerivedTraits::ModelCompatibility ModelCompatibility;                  /**< Indicates which model the source is compatible with. */
    typedef typename Transformation::TransformDataType TransformDataType;                   /**< The error type of the difference between two measurements. */
    typedef typename Transformation::Measurement Measurement;                               /**< The measurement data type. */
    typedef _Derived<State,measurement_type_, DerivedTraits::template Transformation> DerivedSource; /**< The derived source. */
                



    std::function<bool(const State&)> state_in_surveillance_region_callback_;   /**< A pointer to the function which determines if a target's state is inside the source's surveillance region. */
                                                     /**< The source parameters @see SourceParameters. */
    

    template <typename T>
    using SourceTemplate = _Derived<typename State::template StateTemplate<T>,measurement_type_, DerivedTraits::template Transformation>; /**< Used to create a source of a different data type. */



    static MatH H_;
    static MatV V_;                                                                          
 
    /**
     * Copy Constructor.
     */ 
    SourceBase(const SourceBase& other) : SourceBase() {
        params_ = other.params_;
        state_in_surveillance_region_callback_ = other.state_in_surveillance_region_callback_;
    }

    /**
     * Assignment operator
     */ 
    SourceBase& operator = (const SourceBase& other) {
        params_ = other.params_;
        state_in_surveillance_region_callback_ = other.state_in_surveillance_region_callback_;
        return *this;
    }

    /** 
     * Initializes the measurement source by setting the parameters using SetParameters, calculating the non user specified parameters,
     * setting up the callback function, and initializing the Jacobians.
     * @param[in] params The source parameters.
     * @param[in] state_in_surveillance_region_callback The callback function used to determine if a track's state is inside the source's surveillance region. 
    */
    void Init(const SourceParameters& params, std::function<bool(const State&)> state_in_surveillance_region_callback); 

    /** 
     * Initializes the measurement source by setting the parameters using SetParameters, calculating the non user specified parameters,
     * and initializing the Jacobians. By calling this constructor, it is assumed that a track is always inside the surveillance region.
     * This is a good assumption when all of the sources have the same surveillance region.
     * @param[in] params The source parameters.
     */ 
    void Init(const SourceParameters& params);     


    /**
     * Verify that the parameters are valid. If they are, the parameters are set. 
     * @param[in] params Source parameters.
     * \return returns true if the parameters were set; otherwise, false.
     */
    bool SetParameters(const SourceParameters& params); 

    /** 
     * Returns the jacobian of the observation function w.r.t. the states. 
     * If transformation data is provided, the state will be transformed before calculating the Jacobian.
     * Let H denote the Jacobian of the observation function without the state being transformed, and T
     * denote the Jacobian of the transformation. If the state is transformed, the new Jacobian becomes H*T,
     * which this function implements.
     * @param[in] state A state of the target.
     * @param[in] transform_state A flag used to indicate if the state needs to be transformed 
     * @param[in] transform_data The data needed to transform the state
     */
    static MatH GetLinObsMatState(const State& state, const bool transform_state, const TransformDataType& transform_data);
                 

    /** 
     * Returns the jacobian of the observation function w.r.t. the sensor noise.
     * If transformation data is provided, the state will be transformed before calculating
     * the Jacobian.
     * @param[in] state A state of the target.
     * @param[in] transform_state A flag used to indicate if the state needs to be transformed 
     * @param[in] transform_data The data needed to transform the state
     */
    static MatV GetLinObsMatSensorNoise(const State& state, const bool transform_state, const TransformDataType& transform_data);    

    /**
     *  Implements the observation function and returns an estimated measurement based on the state. 
     * If transformation data is provided, the state will be transformed first before calculating the
     * estimated measurement.
     * Currently, the measurement is only given a pose, twist, and measurement type. 
     * @param[in] state A state of the target.
     * @param[in] transform_state A flag used to indicate if the state needs to be transformed 
     * @param[in] transform_data The data needed to transform the state
     */
    static Measurement GetEstMeas(const State& state, const bool transform_state, const TransformDataType& transform_data);

    /**
     * Performs the OMinus operation between two measurement (m1 ominus m2) of the same type. In other words, this
     * method computes the geodesic distance between two measurements of the same type.
     * @param[in] m1 a measurement
     * @param[in] m2 a measurement
     */
    static VecMeas OMinus(const Measurement& m1, const Measurement& m2) {
#ifdef DEBUG_BUILD
    if(m1.type != measurement_type_ || m2.type !=measurement_type_) {
        throw std::runtime_error("SourceBase:: The measurements are not the right type.");
    }
#endif
        return DerivedSource::DerivedOMinus(m1, m2);
    } 


   /**
     * Generates a random measurement from a Gaussian distribution with mean defined by the state and standard deviation defined by meas_std. This
     * method is used primarily in simulations and tests. If transformation data is provided, the state will be transformed by the data before 
     * generating the estimated measurement.
     * @param[in] meas_std The measurement standard deviation.
     * @param[in] state    The state that serves as the mean of the Gaussian distribution.
     * @param[in] transform_state A flag used to indicate if the state needs to be transformed 
     * @param[in] transform_data The data needed to transform the state
     */ 
    Measurement GenerateRandomMeasurement(const MatMeasCov& meas_std, const State& state, const bool transform_state, const TransformDataType& transform_data) const ;

    /**
     * Returns true if the state is inside the source's surveillance region. The state can be transformed into another
     * frame if the transformation data is provided.  
     * @param[in] state A state of the target.
     * @param[in] transform_state A flag used to indicate if the state needs to be transformed 
     * @param[in] transform_data The data needed to transform the state
     */
    bool StateInsideSurveillanceRegion(const State& state, const bool transform_state, const TransformDataType& transform_data) const;

    /**
     * Calculates the temporal distance between two measurements.
     * @param[in] meas1 A measurement.
     * @param[in] meas2 A measurement.
     * @param[in] params The system parameters.
     * \return Returns temporal distance between two measurements
     */
   
    DataType GetTemporalDistance(const Measurement& meas1, const Measurement& meas2, const Parameters& params) const { return fabs(meas1.time_stamp - meas2.time_stamp); }

    /**
     * Calculates the geodesic distance between the pose of two measurements that have the same measurement space.
     * @param[in] meas1 A measurement.
     * @param[in] meas2 A measurement.
     * @param[in] params The system parameters.
     * \return Returns geodesic distance between pose of two measurements
     */
    
    DataType GetSpatialDistance(const Measurement& meas1, const Measurement& meas2, const Parameters& params) const {return gsd_ptr_[meas1.type][meas2.type](meas1,meas2,params);}

    /**
     * Finds the geodesic distance between the pose of two measurements of different time stamps normalized by the temporal distance. The measurements must have the same measurement space.
     * @param[in] meas1 A measurement.
     * @param[in] meas2 A measurement.
     * @param[in] params The system parameters.
     * \return Returns geodesic distance between two measurements
     */
    DataType GetVelocityDistance(const Measurement& meas1, const Measurement& meas2, const Parameters& params) const {
        if (meas1.time_stamp == meas2.time_stamp) {
            throw std::runtime_error("SourceBase::GetVelocityDistance Measurements have the same time stamp");
        } else {
            return GetSpatialDistance(meas1,meas2,params)/GetTemporalDistance(meas1,meas2,params);
        }
    }

    /**
     * Returns a const reference to the source parameters. 
     */ 
    const SourceParameters& GetParams() const {return params_;}

    /**
     * Verifies that the data in the measurement meets certain specifications.
     * @param measurement The measurement to be verified. 
     */ 
    bool IsAcceptableMeasurement(const Measurement& measurement);


// private:
     SourceBase()=default;
    ~SourceBase();
    friend DerivedSource; /**< The object that inherits this base class. */

protected:

     SourceParameters params_;  

private:

   

    /**
     * Ensure that the source parameters meet the specified criteria. If a parameter doesn't, a runtime error will be thrown.
     * @param params The source parameters needed to initialize the source. 
     * \return Returns true if the parameters were successfully set. 
     */ 
    bool VerifySourceParameters(const SourceParameters& params);

    /**
     * The Default callback function used with StateInsideSurveillanceRegion. It always returns true.  
     * @param[in] state A state of the target.
     */
    static bool StateInsideSurveillanceRegionDefaultCallback(const State& state) {
        return true;
    }

    bool source_init_ = false; /**< Indicates if the source has been initialized. */

    /**
     * This array of function pointers, holds pointers to the different methods of calculating the spatial distance
     * between two measurements. The index is determined by the measurement type. If two measurements are of different measurement
     * spaces, then their spatial distance cannot be calculated and a runtime error is thrown.
     */ 
    typedef DataType (*GSDFuncPTR)(const Measurement&, const Measurement&, const Parameters&);

    GSDFuncPTR **gsd_ptr_;

    static DataType GSD_RN_RN_POS(const Measurement& meas1, const Measurement& meas2, const Parameters& params) {return (meas1.pose - meas2.pose).norm();}
    static DataType GSD_SEN_SEN_POSE(const Measurement& meas1, const Measurement& meas2, const Parameters& params){return (State::Group::OMinus(meas1.pose,meas2.pose)).norm(); }
    static DataType GSD_SEN_SEN_POS(const Measurement& meas1, const Measurement& meas2, const Parameters& params){return (meas1.pose - meas2.pose).norm(); }
    static DataType GSD_NotImplemented(const Measurement& meas1, const Measurement& meas2, const Parameters& params){throw std::runtime_error("SourceBase::SpatialDistance Distance not implemented.");}
    static DataType GSD_SE3_CamDepth_SE3_CamDepth(const Measurement& meas1, const Measurement& meas2, const Parameters& params);
    static DataType GSD_R2_R3_Radar(const Measurement& meas1, const Measurement& meas2, const Parameters& params);
    

};

template<typename _SourceDerivedTraits, template<typename , MeasurementTypes , template <typename > typename > typename _Derived>
typename SourceBase<_SourceDerivedTraits,_Derived>::MatH SourceBase<_SourceDerivedTraits,_Derived>::H_;

template<typename _SourceDerivedTraits, template<typename , MeasurementTypes , template <typename > typename > typename _Derived>
typename SourceBase<_SourceDerivedTraits,_Derived>::MatV SourceBase<_SourceDerivedTraits,_Derived>::V_;


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                            Definitions
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename _SourceDerivedTraits, template<typename , MeasurementTypes , template <typename > typename > typename _Derived>
void SourceBase<_SourceDerivedTraits,_Derived>::Init(const SourceParameters& params, std::function<bool(const State&)> state_in_surveillance_region_callback) {
    this->Init(params);
    this->state_in_surveillance_region_callback_ = state_in_surveillance_region_callback;
}

//-------------------------------------------------------------------------------

template<typename _SourceDerivedTraits, template<typename , MeasurementTypes , template <typename > typename > typename _Derived>
void SourceBase<_SourceDerivedTraits,_Derived>::Init(const SourceParameters& params) {

    bool success = SetParameters(params); // Verifies the parameters. If there is an invalid parameter, an error will be thrown. Otherwise, the parameters are set.
    this->state_in_surveillance_region_callback_ = StateInsideSurveillanceRegionDefaultCallback;
    static_cast<DerivedSource*>(this)->DerivedInit(params_);


    if (!source_init_) {

        // Generate two dimensional array of function pointers.
        gsd_ptr_ = new GSDFuncPTR *[MeasurementTypes::NUM_TYPES];
        for(int i = 0; i < MeasurementTypes::NUM_TYPES; ++i)
        {
            gsd_ptr_[i] = new GSDFuncPTR[MeasurementTypes::NUM_TYPES];
        }

        // Set each function pointer to null
        for (int i = 0; i < MeasurementTypes::NUM_TYPES; ++i) {
            for (int j = 0; j < MeasurementTypes::NUM_TYPES; ++j) {
                gsd_ptr_[i][j] = &GSD_NotImplemented;
            }
        }
    
        gsd_ptr_[MeasurementTypes::RN_POS][MeasurementTypes::RN_POS]                 = &GSD_RN_RN_POS;
        gsd_ptr_[MeasurementTypes::RN_POS][MeasurementTypes::RN_POS_VEL]             = &GSD_RN_RN_POS;
        gsd_ptr_[MeasurementTypes::RN_POS_VEL][MeasurementTypes::RN_POS]             = &GSD_RN_RN_POS;
        gsd_ptr_[MeasurementTypes::RN_POS_VEL][MeasurementTypes::RN_POS_VEL]         = &GSD_RN_RN_POS;
        gsd_ptr_[MeasurementTypes::SEN_POSE][MeasurementTypes::SEN_POSE]             = &GSD_SEN_SEN_POSE;
        gsd_ptr_[MeasurementTypes::SEN_POSE][MeasurementTypes::SEN_POSE_TWIST]       = &GSD_SEN_SEN_POSE;
        gsd_ptr_[MeasurementTypes::SEN_POSE_TWIST][MeasurementTypes::SEN_POSE_TWIST] = &GSD_SEN_SEN_POSE;
        gsd_ptr_[MeasurementTypes::SEN_POSE_TWIST][MeasurementTypes::SEN_POSE]       = &GSD_SEN_SEN_POSE;
        gsd_ptr_[MeasurementTypes::SEN_POS][MeasurementTypes::SEN_POS]               = &GSD_SEN_SEN_POS;
        gsd_ptr_[MeasurementTypes::SEN_POS][MeasurementTypes::SEN_POS_VEL]           = &GSD_SEN_SEN_POS;
        gsd_ptr_[MeasurementTypes::SEN_POS_VEL][MeasurementTypes::SEN_POS_VEL]       = &GSD_SEN_SEN_POS;
        gsd_ptr_[MeasurementTypes::SEN_POS_VEL][MeasurementTypes::SEN_POS]           = &GSD_SEN_SEN_POS;
        gsd_ptr_[MeasurementTypes::SE3_CAM_DEPTH][MeasurementTypes::SE3_CAM_DEPTH]   = &GSD_SE3_CamDepth_SE3_CamDepth;
        gsd_ptr_[MeasurementTypes::R2_R3_RADAR][MeasurementTypes::R2_R3_RADAR]                           = &GSD_R2_R3_Radar;
        gsd_ptr_[MeasurementTypes::R2_R3_RADAR][MeasurementTypes::R2_R3_RADAR_DEPTH_DERIV]               = &GSD_R2_R3_Radar;
        gsd_ptr_[MeasurementTypes::R2_R3_RADAR_DEPTH_DERIV][MeasurementTypes::R2_R3_RADAR]               = &GSD_R2_R3_Radar;
        gsd_ptr_[MeasurementTypes::R2_R3_RADAR_DEPTH_DERIV][MeasurementTypes::R2_R3_RADAR_DEPTH_DERIV]   = &GSD_R2_R3_Radar;



    }
    


    source_init_ = true;
}   



//---------------------------------------------------

template<typename _SourceDerivedTraits, template<typename , MeasurementTypes , template <typename > typename > typename _Derived>
SourceBase<_SourceDerivedTraits,_Derived>::~SourceBase() {

    if(source_init_) {
        for (int i = 0; i < MeasurementTypes::NUM_TYPES; i++) {
            delete [] gsd_ptr_[i];
        }
        delete [] gsd_ptr_;
    }

}

//---------------------------------------------------

template<typename _SourceDerivedTraits, template<typename , MeasurementTypes , template <typename > typename > typename _Derived>
bool SourceBase<_SourceDerivedTraits,_Derived>::VerifySourceParameters(const SourceParameters& params) {

    bool success = true;
    unsigned int mult = 1;

 

    if (params.meas_cov_.rows() !=total_meas_dim_  ) { // Make sure that it is not empty
        throw std::runtime_error("SourceBase::VerifySourceParameters: Measurement covariance is not the right dimension");
        success = false;
    } 

    if (  ((params.meas_cov_ + params.meas_cov_.transpose())/2.0 - params.meas_cov_).norm() > 1e-12 ) {
        throw std::runtime_error("SourceBase::VerifySourceParameters: Measurement covariance is not symmetic. ");
        success = false;

    }

    Eigen::VectorXcd eigen_values = params.meas_cov_.eigenvalues();
    for (int ii =0; ii < eigen_values.rows(); ++ii){                       // positive definite
        if(std::real(eigen_values(ii)) <0) {
            throw std::runtime_error("SourceBase::VerifySourceParameters: Measurement covariance is not positive definite. ");
          success = false;

        }
    }


    // Expected number of false measurements
    if (params.spacial_density_of_false_meas_ <= 0 ) {
        throw std::runtime_error("SourceBase::VerifySourceParameters: The spacial_density_of_false_meas_ of false measurements must be greater than 0 . ");
        success = false;

    }

    // Verify the number of measurement types
    if(params.type_ != measurement_type_) {
        throw std::runtime_error("SourceBase::VerifySourceParameters: The measurement type doesn't match. ");
        success = false;
    }


    // Verify the probability of detection 
    if (params.probability_of_detection_ < 0 || params.probability_of_detection_ > 1) {
        throw std::runtime_error("SourceBase::VerifySourceParameters: probability_of_detection_ must be between 0 and 1. ");
        success = false;

    }

    // Verify the gate probability
    if (params.gate_probability_ < 0 || params.gate_probability_ > 1) {
        throw std::runtime_error("SourceBase::VerifySourceParameters: gate_probability_ must be between 0 and 1. ");
        success = false;

    }



    return success;

}

//------------------------------------------------------------------------------------------------------------------------
template<typename _SourceDerivedTraits, template<typename , MeasurementTypes , template <typename > typename > typename _Derived>
bool SourceBase<_SourceDerivedTraits,_Derived>::SetParameters(const SourceParameters& params) {
    bool success = VerifySourceParameters(params);
    if (success) {
        this->params_ = params;
        boost::math::chi_squared dist(total_meas_dim_);
        if (this->params_.gate_threshold_ <= 0) { // not set by user so calculate it
            this->params_.gate_threshold_ = boost::math::quantile(dist, this->params_.gate_probability_);
        }
        this->params_.vol_unit_hypershpere_ = pow(M_PI, static_cast<double>(total_meas_dim_)/2.0)/boost::math::tgamma(static_cast<double>(total_meas_dim_)/2.0 +1.0);
        this->params_.gate_threshold_sqrt_ = sqrt(this->params_.gate_threshold_ ); 

        this->params_.meas_pose_rows_ = meas_pose_rows_;       
        this->params_.meas_pose_cols_ = meas_pose_cols_;      
        this->params_.meas_twist_rows_ = meas_twist_rows_;      
        this->params_.meas_twist_cols_ = meas_twist_cols_;       
        this->params_.meas_pose_dim_ = meas_pose_dim_;         
        this->params_.meas_twist_dim_ = meas_twist_dim_;       
        this->params_.total_meas_dim_ = total_meas_dim_;       
        this->params_.has_vel_ = has_vel_;                     

    }
    return success;
}

//------------------------------------------------------------------------------------------------------------------------

template<typename _SourceDerivedTraits, template<typename , MeasurementTypes , template <typename > typename > typename _Derived>
typename SourceBase<_SourceDerivedTraits,_Derived>::Measurement SourceBase<_SourceDerivedTraits,_Derived>::GetEstMeas(const State& state, const bool transform_state, const TransformDataType& transform_data)  {
        
#ifdef DEBUG_BUILD
    if(transform_state && Transformation::is_null_transform_) {
        throw std::runtime_error("SourceBase::GetEstMeas Trying to transform the state when the transform object is TransformNULL");
    } else if(transform_state && transform_data.rows() == 0) {
        throw std::runtime_error("SourceBase::GetEstMeas Trying to transform the state when transform data is empty");

    }
#endif
    
    if(transform_state && !Transformation::is_null_transform_) {
        State state_transformed = Transformation::TransformState(state,transform_data);
        return DerivedSource::DerivedGetEstMeas(state_transformed);

    } else {
        return DerivedSource::DerivedGetEstMeas(state);
    }
} 


//------------------------------------------------------------------------------------------------------------------------
template<typename _SourceDerivedTraits, template<typename , MeasurementTypes , template <typename > typename > typename _Derived>
typename SourceBase<_SourceDerivedTraits,_Derived>::MatH SourceBase<_SourceDerivedTraits,_Derived>::GetLinObsMatState(const State& state, const bool transform_state, const TransformDataType& transform_data) {
    if(transform_state && !Transformation::is_null_transform_) {
        State state_transformed = Transformation::TransformState(state,transform_data);
        typename Transformation::MatCov transform_jacobian = Transformation::GetTransformationJacobian(state,transform_data);
        return DerivedSource::DerivedGetLinObsMatState(state_transformed)*transform_jacobian;

    } else {
        return DerivedSource::DerivedGetLinObsMatState(state);
    }        
}    

//------------------------------------------------------------------------------------------------------------------------

template<typename _SourceDerivedTraits, template<typename , MeasurementTypes , template <typename > typename > typename _Derived>
typename SourceBase<_SourceDerivedTraits,_Derived>::MatV SourceBase<_SourceDerivedTraits,_Derived>::GetLinObsMatSensorNoise(const State& state, const bool transform_state, const TransformDataType& transform_data) {
    if(transform_state && !Transformation::is_null_transform_) {
        State state_transformed = Transformation::TransformState(state,transform_data);
        return DerivedSource::DerivedGetLinObsMatSensorNoise(state_transformed);

    } else {
        return DerivedSource::DerivedGetLinObsMatSensorNoise(state);
    }
}  

//------------------------------------------------------------------------------------------------------------------------

template<typename _SourceDerivedTraits, template<typename , MeasurementTypes , template <typename > typename > typename _Derived>
typename SourceBase<_SourceDerivedTraits,_Derived>::Measurement SourceBase<_SourceDerivedTraits,_Derived>::GenerateRandomMeasurement(const MatMeasCov& meas_std, const State& state, const bool transform_state, const TransformDataType& transform_data) const {
    if(transform_state && !Transformation::is_null_transform_) {
        State state_transformed = Transformation::TransformState(state,transform_data);
        return static_cast<const DerivedSource*>(this)->DerivedGenerateRandomMeasurement(meas_std,state_transformed);

    } else {
        return static_cast<const DerivedSource*>(this)->DerivedGenerateRandomMeasurement(meas_std,state);
    }
}

//------------------------------------------------------------------------------------------------------------------------

template<typename _SourceDerivedTraits, template<typename , MeasurementTypes , template <typename > typename > typename _Derived>
bool SourceBase<_SourceDerivedTraits,_Derived>::StateInsideSurveillanceRegion(const State& state, const bool transform_state, const TransformDataType& transform_data) const {
    
    if(transform_state && !Transformation::is_null_transform_) {
        State state_transformed = Transformation::TransformState(state,transform_data);
        return state_in_surveillance_region_callback_(state_transformed);

    } else {
        return state_in_surveillance_region_callback_(state);
    }
}

//------------------------------------------------------------------------------------------------------------------------

template<typename _SourceDerivedTraits, template<typename , MeasurementTypes , template <typename > typename > typename _Derived>
bool SourceBase<_SourceDerivedTraits,_Derived>::IsAcceptableMeasurement(const Measurement& measurement) {

    bool success = true;

        if (measurement.type != measurement_type_) {
            throw std::runtime_error("SourceBase::IsAcceptableMeasurement Measurement type does not match the source's measurement type. Make sure the source index is correct.");
            success = false;            
        } else if (measurement.pose.rows() != meas_pose_rows_ || measurement.pose.cols() != meas_pose_cols_) {
            throw std::runtime_error("SourceBase::IsAcceptableMeasurement The pose of the measurement is not the correct dimension.");
            success = false;
        } else if (has_vel_ && (measurement.twist.rows() != meas_twist_rows_ || measurement.twist.cols() != meas_twist_cols_)) {
            throw std::runtime_error("SourceBase::IsAcceptableMeasurement The twist of the measurement is not the correct dimension.");
            success = false;
        } else if(measurement.transform_state) {
            if(!Transformation::IsAcceptableTransformData(measurement.transform_data_t_m)) {
                throw std::runtime_error("SourceBase::IsAcceptableMeasurement The transformation data is not correct.");
                success = false;
            }
        } else {
            success = true;
        }
    return success;
}

//------------------------------------------------------------------------------------------------------------------------

template<typename _SourceDerivedTraits, template<typename , MeasurementTypes , template <typename > typename > typename _Derived>
typename SourceBase<_SourceDerivedTraits,_Derived>::DataType SourceBase<_SourceDerivedTraits,_Derived>::GSD_SE3_CamDepth_SE3_CamDepth(const Measurement& meas1, const Measurement& meas2, const Parameters& params){
   
    typedef Eigen::Matrix<DataType,meas_pose_dim_,1> VecPose;
    DataType d = 0;
    VecPose pose1, pose2;

    if(meas1.transform_meas && meas2.transform_meas) {
        // Transform measurements into common frame
        pose1 = Transformation::TransformMeasurement(meas1,meas1.transform_data_m_t).pose;
        pose2 = Transformation::TransformMeasurement(meas2,meas2.transform_data_m_t).pose;
    } else if(meas1.transform_meas) {
        pose1 = Transformation::TransformMeasurement(meas1,meas1.transform_data_m_t).pose;
        pose2 = meas2.pose;
    } else if(meas2.transform_meas) {
        pose1 = meas1.pose;
        pose2 = Transformation::TransformMeasurement(meas2,meas2.transform_data_m_t).pose;
    } else {
        pose1 = meas1.pose;
        pose2 = meas2.pose;
    }

    d = (pose1.block(1,0,3,1)*pose1(0) - pose2.block(1,0,3,1)*pose2(0)).norm();

    return d;
}

//------------------------------------------------------------------------------------------------------------------------

template<typename _SourceDerivedTraits, template<typename , MeasurementTypes , template <typename > typename > typename _Derived>
typename SourceBase<_SourceDerivedTraits,_Derived>::DataType SourceBase<_SourceDerivedTraits,_Derived>::GSD_R2_R3_Radar(const Measurement& meas1, const Measurement& meas2, const Parameters& params){ 


    typedef Eigen::Matrix<DataType,meas_pose_dim_,1> VecPose;

    DataType d = 0;
    VecPose err, pose1, pose2;

   


    // The radar sensors are assumed to be stationary, so if they are from the same sensor,
    // we can compare them directly; otherwise, we must transform them into a common frame.
    if(meas1.transform_meas && meas2.transform_meas) {
        // Transform measurements into common frame
        pose1 = Transformation::TransformMeasurement(meas1,meas1.transform_data_m_t).pose;
        pose2 = Transformation::TransformMeasurement(meas2,meas2.transform_data_m_t).pose;
        
    } else if(meas1.transform_meas) {
        pose1 = Transformation::TransformMeasurement(meas1,meas1.transform_data_m_t).pose;
        pose2 = meas2.pose;
    } else if(meas2.transform_meas) {
        pose1 = meas1.pose;
        pose2 = Transformation::TransformMeasurement(meas2,meas2.transform_data_m_t).pose;
    } else {
        pose1 = meas1.pose;
        pose2 = meas2.pose;
    }

    // Transform pose to cartesian coordinates
    const DataType& r1 = pose1(0);
    const DataType& r2 = pose2(0);
    const DataType& azimuth1 = pose1(1);
    const DataType& azimuth2 = pose2(1);
    if (meas_pose_dim_ == 2) {
       d = sqrt(r1*r1 + r2*r2 - 2.0*r1*r2*cos(azimuth1-azimuth2));

    } else if(meas_pose_dim_ == 3) {
        const DataType& zenith1 = pose1(2);
        const DataType& zenith2 = pose2(2);

        d = sqrt( r1*r1 + r2*r2 - 2.0*r1*r2*( cos(zenith1)*cos(zenith2) + sin(zenith1)*sin(zenith2)*cos(azimuth1-azimuth2)));
    } else {
        throw std::runtime_error("SourceBase::GSD_R2_R3_Radar, measurement dimension isn't correct.");
    }



    return d;
}

} // namespace rransac



#endif // RRANSAC_COMMON_SOURCES_SOURCE_BASE_H_