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
    bool has_twist;                      /**< Indicates if the measurement has twist data in addition to pose. This is calculated from SourceParameters::type_. */
    double gate_threshold_sqrt_;         /**< The square root of the gate threshold. */    
    double vol_unit_hypershpere_;        /**< The volume of the unit hypershpere in the measurement space. This is determined by the dimension of the measurement source. */
    unsigned int meas_space_dim_;        /**< The dimension of the measurement space. */
    unsigned int meas_space_dim_mult_;   /**< Multiplier used on the measurement space incase the measurement space includes a derivative. */

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




template<typename tState, MeasurementTypes tMeasurementType, template <typename > typename tTransformation, int tMeasSpaceDim, template <typename , MeasurementTypes , template <typename > typename > typename tDerived>
class SourceBase
{

public:

    typedef tState State;                                                       /**< The state of the target. @see State. */
    typedef tTransformation<tState> Transformation;                             /**< The transformation used to transform the measurements and tracks. */
    typedef tDerived<tState,tMeasurementType,tTransformation> DerivedSource;    /**< The derived source. */
    typedef typename tState::DataType DataType;                                 /**< The scalar object for the data. Ex. float, double, etc. */
    typedef Eigen::Matrix<DataType,Eigen::Dynamic,Eigen::Dynamic> MatXd;        /**< The object type of the Jacobians. */
    static constexpr unsigned int meas_space_dim_ = tMeasSpaceDim;              /**< The Dimension of the measurement space. */
    static constexpr MeasurementTypes meas_type_ = tMeasurementType;
    static constexpr int meas_space_dim_mult_ = MeasDimMultiplier<tMeasurementType>::value; /**< a constant used when the measurement contains velocity. */
    static constexpr int has_vel_ = MeasHasVelocity<tMeasurementType>::value;               /**< Indicates if the measurement contains velocity.  */

    std::function<bool(const State&)> state_in_surveillance_region_callback_;   /**< A pointer to the function which determines if a target's state is inside the source's surveillance region. */
                                                     /**< The source parameters @see SourceParameters. */
    

    template <typename tDataType>
    using SourceTemplate = tDerived<typename tState::template StateTemplate<tDataType>,tMeasurementType,tTransformation>; /**< Used to create a source of a different data type. */



    static MatXd H_;
    static MatXd V_;                                                                          
 
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
    static MatXd GetLinObsMatState(const State& state, const bool transform_state, const MatXd& transform_data);
                 

    /** 
     * Returns the jacobian of the observation function w.r.t. the sensor noise.
     * If transformation data is provided, the state will be transformed before calculating
     * the Jacobian.
     * @param[in] state A state of the target.
     * @param[in] transform_state A flag used to indicate if the state needs to be transformed 
     * @param[in] transform_data The data needed to transform the state
     */
    static MatXd GetLinObsMatSensorNoise(const State& state, const bool transform_state, const MatXd& transform_data);    

    /**
     *  Implements the observation function and returns an estimated measurement based on the state. 
     * If transformation data is provided, the state will be transformed first before calculating the
     * estimated measurement.
     * Currently, the measurement is only given a pose, twist, and measurement type. 
     * @param[in] state A state of the target.
     * @param[in] transform_state A flag used to indicate if the state needs to be transformed 
     * @param[in] transform_data The data needed to transform the state
     */
    static Meas<DataType> GetEstMeas(const State& state, const bool transform_state, const MatXd& transform_data);

    /**
     * Performs the OMinus operation between two measurement (m1 ominus m2) of the same type. In other words, this
     * method computes the geodesic distance between two measurements of the same type.
     * @param[in] m1 a measurement
     * @param[in] m2 a measurement
     */
    static MatXd OMinus(const Meas<DataType>& m1, const Meas<DataType>& m2) {
#ifdef DEBUG_BUILD
    if(m1.type != tMeasurementType || m2.type !=tMeasurementType) {
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
    Meas<DataType> GenerateRandomMeasurement(const MatXd& meas_std, const State& state, const bool transform_state, const MatXd& transform_data) const ;

    /**
     * Returns true if the state is inside the source's surveillance region. The state can be transformed into another
     * frame if the transformation data is provided.  
     * @param[in] state A state of the target.
     * @param[in] transform_state A flag used to indicate if the state needs to be transformed 
     * @param[in] transform_data The data needed to transform the state
     */
    bool StateInsideSurveillanceRegion(const State& state, const bool transform_state, const MatXd& transform_data) const;

    /**
     * Calculates the temporal distance between two measurements.
     * @param[in] meas1 A measurement.
     * @param[in] meas2 A measurement.
     * @param[in] params The system parameters.
     * \return Returns temporal distance between two measurements
     */
   
    DataType GetTemporalDistance(const Meas<DataType>& meas1, const Meas<DataType>& meas2, const Parameters& params) const { return fabs(meas1.time_stamp - meas2.time_stamp); }

    /**
     * Calculates the geodesic distance between the pose of two measurements that have the same measurement space.
     * @param[in] meas1 A measurement.
     * @param[in] meas2 A measurement.
     * @param[in] params The system parameters.
     * \return Returns geodesic distance between pose of two measurements
     */
    
    DataType GetSpatialDistance(const Meas<DataType>& meas1, const Meas<DataType>& meas2, const Parameters& params) const {return gsd_ptr_[meas1.type][meas2.type](meas1,meas2,params);}

    /**
     * Finds the geodesic distance between the pose of two measurements of different time stamps normalized by the temporal distance. The measurements must have the same measurement space.
     * @param[in] meas1 A measurement.
     * @param[in] meas2 A measurement.
     * @param[in] params The system parameters.
     * \return Returns geodesic distance between two measurements
     */
    DataType GetVelocityDistance(const Meas<DataType>& meas1, const Meas<DataType>& meas2, const Parameters& params) const {
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
    typedef DataType (*GSDFuncPTR)(const Meas<DataType>&, const Meas<DataType>&, const Parameters&);

    GSDFuncPTR **gsd_ptr_;

    static DataType GSD_RN_RN_POS(const Meas<DataType>& meas1, const Meas<DataType>& meas2, const Parameters& params) {return (meas1.pose - meas2.pose).norm();}
    static DataType GSD_SEN_SEN_POSE(const Meas<DataType>& meas1, const Meas<DataType>& meas2, const Parameters& params){return (State::Group::OMinus(meas1.pose,meas2.pose)).norm(); }
    static DataType GSD_SEN_SEN_POS(const Meas<DataType>& meas1, const Meas<DataType>& meas2, const Parameters& params){return (meas1.pose - meas2.pose).norm(); }
    static DataType GSD_NotImplemented(const Meas<DataType>& meas1, const Meas<DataType>& meas2, const Parameters& params){throw std::runtime_error("SourceBase::SpatialDistance Distance not implemented.");}
    static DataType GSD_SE3_CamDepth_SE3_CamDepth(const Meas<DataType>& meas1, const Meas<DataType>& meas2, const Parameters& params);
    

};


template<typename tState, MeasurementTypes tMeasurementType, template <typename > typename tTransformation, int tMeasSpaceDim, template <typename , MeasurementTypes , template <typename > typename > typename tDerived>
Eigen::Matrix<typename tState::DataType,Eigen::Dynamic,Eigen::Dynamic> SourceBase<tState,tMeasurementType,tTransformation,tMeasSpaceDim,tDerived>::H_;

template<typename tState, MeasurementTypes tMeasurementType, template <typename > typename tTransformation, int tMeasSpaceDim, template <typename , MeasurementTypes , template <typename > typename > typename tDerived>
Eigen::Matrix<typename tState::DataType,Eigen::Dynamic,Eigen::Dynamic> SourceBase<tState,tMeasurementType,tTransformation,tMeasSpaceDim,tDerived>::V_;


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                            Definitions
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename tState, MeasurementTypes tMeasurementType, template <typename > typename tTransformation, int tMeasSpaceDim, template <typename , MeasurementTypes , template <typename > typename > typename tDerived>
void SourceBase<tState,tMeasurementType,tTransformation,tMeasSpaceDim,tDerived>::Init(const SourceParameters& params, std::function<bool(const State&)> state_in_surveillance_region_callback) {
    this->Init(params);
    this->state_in_surveillance_region_callback_ = state_in_surveillance_region_callback;
}

//-------------------------------------------------------------------------------

template<typename tState, MeasurementTypes tMeasurementType, template <typename > typename tTransformation, int tMeasSpaceDim, template <typename , MeasurementTypes , template <typename > typename > typename tDerived>
void SourceBase<tState,tMeasurementType,tTransformation,tMeasSpaceDim,tDerived>::Init(const SourceParameters& params) {

    bool success = SetParameters(params); // Verifies the parameters. If there is an invalid parameter, an error will be thrown. Otherwise, the parameters are set.
    this->state_in_surveillance_region_callback_ = StateInsideSurveillanceRegionDefaultCallback;
    static_cast<DerivedSource*>(this)->DerivedInit(params_);


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
    


    bool source_init_ = true;
}   



//---------------------------------------------------

template<typename tState, MeasurementTypes tMeasurementType, template <typename > typename tTransformation, int tMeasSpaceDim, template <typename , MeasurementTypes , template <typename > typename > typename tDerived>
SourceBase<tState,tMeasurementType,tTransformation,tMeasSpaceDim,tDerived>::~SourceBase() {

    if(source_init_) {
        for (int i = 0; i < MeasurementTypes::NUM_TYPES; i++) {
            delete [] gsd_ptr_[i];
        }
        delete [] gsd_ptr_;
    }

}

//---------------------------------------------------

template<typename tState, MeasurementTypes tMeasurementType, template <typename > typename tTransformation, int tMeasSpaceDim, template <typename , MeasurementTypes , template <typename > typename > typename tDerived>
bool SourceBase<tState,tMeasurementType,tTransformation,tMeasSpaceDim,tDerived>::VerifySourceParameters(const SourceParameters& params) {

    bool success = true;
    unsigned int mult = 1;

 

    if (params.meas_cov_.rows() != DerivedSource::meas_space_dim_*meas_space_dim_mult_ ) { // Make sure that it is not empty
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
    if(params.type_ != meas_type_) {
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
template<typename tState, MeasurementTypes tMeasurementType, template <typename > typename tTransformation, int tMeasSpaceDim, template <typename , MeasurementTypes , template <typename > typename > typename tDerived>
bool SourceBase<tState,tMeasurementType,tTransformation,tMeasSpaceDim,tDerived>::SetParameters(const SourceParameters& params) {
    bool success = VerifySourceParameters(params);
    if (success) {
        this->params_ = params;
        if(this->params_.type_ == MeasurementTypes::RN_POS_VEL || this->params_.type_ == MeasurementTypes::SEN_POS_VEL || this->params_.type_ == MeasurementTypes::SEN_POSE_TWIST) {
            this->params_.has_twist = true;
            boost::math::chi_squared dist(meas_space_dim_*2);
            if (this->params_.gate_threshold_ <= 0) { // not set by user so calculate it
                this->params_.gate_threshold_ = boost::math::quantile(dist, this->params_.gate_probability_);
            }
            this->params_.vol_unit_hypershpere_ = pow(M_PI, static_cast<double>(meas_space_dim_))/boost::math::tgamma(static_cast<double>(meas_space_dim_) +1.0);
        } else {
            this->params_.has_twist = false;
            boost::math::chi_squared dist(meas_space_dim_);
            if (this->params_.gate_threshold_ <= 0) { // not set by user so calculate it
                this->params_.gate_threshold_ = boost::math::quantile(dist, this->params_.gate_probability_);
            }
            this->params_.vol_unit_hypershpere_ = pow(M_PI, static_cast<double>(meas_space_dim_)/2.0)/boost::math::tgamma(static_cast<double>(meas_space_dim_)/2.0 +1.0);
        }
    
        this->params_.gate_threshold_sqrt_ = sqrt(this->params_.gate_threshold_ ); 
        this->params_.meas_space_dim_ = meas_space_dim_;
        this->params_.meas_space_dim_mult_ = meas_space_dim_mult_;
    }
    return success;
}

//------------------------------------------------------------------------------------------------------------------------

template<typename tState, MeasurementTypes tMeasurementType, template <typename > typename tTransformation, int tMeasSpaceDim, template <typename , MeasurementTypes , template <typename > typename > typename tDerived>
Meas<typename tState::DataType> SourceBase<tState,tMeasurementType,tTransformation,tMeasSpaceDim,tDerived>::GetEstMeas(const State& state, const bool transform_state, const MatXd& transform_data)  {
        
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
template<typename tState, MeasurementTypes tMeasurementType, template <typename > typename tTransformation, int tMeasSpaceDim, template <typename , MeasurementTypes , template <typename > typename > typename tDerived>
Eigen::Matrix<typename tState::DataType, Eigen::Dynamic, Eigen::Dynamic> SourceBase<tState,tMeasurementType,tTransformation,tMeasSpaceDim,tDerived>::GetLinObsMatState(const State& state, const bool transform_state, const MatXd& transform_data) {
    if(transform_state && !Transformation::is_null_transform_) {
        State state_transformed = Transformation::TransformState(state,transform_data);
        MatXd transform_jacobian = Transformation::GetTransformationJacobian(state,transform_data);
        return DerivedSource::DerivedGetLinObsMatState(state_transformed)*transform_jacobian;

    } else {
        return DerivedSource::DerivedGetLinObsMatState(state);
    }        
}    

//------------------------------------------------------------------------------------------------------------------------

template<typename tState, MeasurementTypes tMeasurementType, template <typename > typename tTransformation, int tMeasSpaceDim, template <typename , MeasurementTypes , template <typename > typename > typename tDerived>
Eigen::Matrix<typename tState::DataType, Eigen::Dynamic, Eigen::Dynamic> SourceBase<tState,tMeasurementType,tTransformation,tMeasSpaceDim,tDerived>::GetLinObsMatSensorNoise(const State& state, const bool transform_state, const MatXd& transform_data) {
    if(transform_state && !Transformation::is_null_transform_) {
        State state_transformed = Transformation::TransformState(state,transform_data);
        return DerivedSource::DerivedGetLinObsMatSensorNoise(state_transformed);

    } else {
        return DerivedSource::DerivedGetLinObsMatSensorNoise(state);
    }
}  

//------------------------------------------------------------------------------------------------------------------------

template<typename tState, MeasurementTypes tMeasurementType, template <typename > typename tTransformation, int tMeasSpaceDim, template <typename , MeasurementTypes , template <typename > typename > typename tDerived>
Meas<typename tState::DataType> SourceBase<tState,tMeasurementType,tTransformation,tMeasSpaceDim,tDerived>::GenerateRandomMeasurement(const MatXd& meas_std, const State& state, const bool transform_state, const MatXd& transform_data) const {
    if(transform_state && !Transformation::is_null_transform_) {
        State state_transformed = Transformation::TransformState(state,transform_data);
        return static_cast<const DerivedSource*>(this)->DerivedGenerateRandomMeasurement(meas_std,state_transformed);

    } else {
        return static_cast<const DerivedSource*>(this)->DerivedGenerateRandomMeasurement(meas_std,state);
    }
}

//------------------------------------------------------------------------------------------------------------------------

template<typename tState, MeasurementTypes tMeasurementType, template <typename > typename tTransformation, int tMeasSpaceDim, template <typename , MeasurementTypes , template <typename > typename > typename tDerived>
bool SourceBase<tState,tMeasurementType,tTransformation,tMeasSpaceDim,tDerived>::StateInsideSurveillanceRegion(const State& state, const bool transform_state, const MatXd& transform_data) const {
    
    if(transform_state && !Transformation::is_null_transform_) {
        State state_transformed = Transformation::TransformState(state,transform_data);
        return state_in_surveillance_region_callback_(state_transformed);

    } else {
        return state_in_surveillance_region_callback_(state);
    }
}

//------------------------------------------------------------------------------------------------------------------------

template<typename tState, MeasurementTypes tMeasurementType, template <typename > typename tTransformation, int tMeasSpaceDim, template <typename , MeasurementTypes , template <typename > typename > typename tDerived>
typename tState::DataType SourceBase<tState,tMeasurementType,tTransformation,tMeasSpaceDim,tDerived>::GSD_SE3_CamDepth_SE3_CamDepth(const Meas<DataType>& meas1, const Meas<DataType>& meas2, const Parameters& params){
    
    DataType d = 0;

    if(meas1.transform_state && meas2.transform_state) {

        d = (meas1.transform_data_m_t.block(0,0,3,3)*meas1.pose(0,0)*meas1.pose.block(1,0,3,1) + meas1.transform_data_m_t.block(0,3,3,1) - meas2.transform_data_m_t.block(0,0,3,3)*meas2.pose(0,0)*meas2.pose.block(1,0,3,1) - meas2.transform_data_m_t.block(0,3,3,1)).norm();

    } else if(meas1.transform_state) {

        d = (meas1.transform_data_m_t.block(0,0,3,3)*meas1.pose(0,0)*meas1.pose.block(1,0,3,1) + meas1.transform_data_m_t.block(0,3,3,1) - meas2.pose(0,0)*meas2.pose.block(1,0,3,1)).norm();

    } else if(meas2.transform_state) {
        
        d = (meas1.pose(0,0)*meas1.pose.block(1,0,3,1) - meas2.transform_data_m_t.block(0,0,3,3)*meas2.pose(0,0)*meas2.pose.block(1,0,3,1) - meas2.transform_data_m_t.block(0,3,3,1)).norm();

    } else {
        d = (meas1.pose - meas2.pose).norm();
    }

    return d;
    }


} // namespace rransac



#endif // RRANSAC_COMMON_SOURCES_SOURCE_BASE_H_