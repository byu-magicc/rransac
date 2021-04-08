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
    // TODO:: I should have the user specify if it has twist or not. Add it as a template parameter to optimize things. 
    bool has_twist;                      /**< Indicates if the measurement has twist data in addition to pose. This is calculated from SourceParameters::type_. */
    double gate_threshold_sqrt_;         /**< The square root of the gate threshold. */    
    double vol_unit_hypershpere_;        /**< The volume of the unit hypershpere in the measurement space. This is determined by the dimension of the measurement source. */

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

template<typename tState, typename tDerived>
class SourceBase
{

public:

    typedef tState State;                                                       /**< The state of the target. @see State. */
    typedef typename tState::DataType DataType;                                 /**< The scalar object for the data. Ex. float, double, etc. */
    typedef Eigen::Matrix<DataType,Eigen::Dynamic,Eigen::Dynamic> MatXd;        /**< The object type of the Jacobians. */
    static constexpr unsigned int meas_space_dim_ = tDerived::meas_space_dim_;  /**< The Dimension of the measurement space. */
   
    std::function<bool(const State&)> state_in_surveillance_region_callback_;   /**< A pointer to the function which determines if a target's state is inside the source's surveillance region. */
    SourceParameters params_;                                                   /**< The source parameters @see SourceParameters. */
    

                                                                                    
    MatXd H_; /**< The Jacobian of the observation function w.r.t. the state. */
    MatXd V_; /**< The Jacobian of the observation function w.r.t. the measurement noise. */

    /**
     * Copy Constructor.
     */ 
    SourceBase(const SourceBase& other) : SourceBase() {
        params_ = other.params_;
        H_ = other.H_;
        V_ = other.V_;
        state_in_surveillance_region_callback_ = other.state_in_surveillance_region_callback_;
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
     * This is an optimized function that returns the jacobian of the observation function w.r.t. the states. 
     * @param[in] state A state of the target.
    */
    MatXd GetLinObsMatState(const State& state) const {
        return static_cast<const tDerived*>(this)->DerivedGetLinObsMatState(state);
    }    

    /** 
     * Returns the jacobian of the observation function w.r.t. the states. 
     * @param[in] state A state of the target.
    */
    static MatXd GetLinObsMatState(const State& state, const MeasurementTypes type) {
        return tDerived::DerivedGetLinObsMatState(state, type);
    }                            

    /** 
     * This is an optimized function that returns the jacobian of the observation function w.r.t. the sensor noise.
     * @param[in] state A state of the target.
     */
    MatXd GetLinObsMatSensorNoise(const State& state) const {
        return static_cast<const tDerived*>(this)->DerivedGetLinObsMatSensorNoise(state);
    }       

    /** 
     * Returns the jacobian of the observation function w.r.t. the sensor noise.
     * @param[in] state A state of the target.
     */
    static MatXd GetLinObsMatSensorNoise(const State& state, const MeasurementTypes type)  {
        return tDerived::DerivedGetLinObsMatSensorNoise(state, type);
    }                      

    /** 
     * This is an optimized function that implements the observation function and returns an estimated measurement based on the state.
     * @param[in] state A state of the target.
     */
    Meas<DataType> GetEstMeas(const State& state) const {
        return static_cast<const tDerived*>(this)->DerivedGetEstMeas(state);
    } 

    /**
     *  Implements the observation function and returns an estimated measurement based on the state. 
     * @param[in] state A state of the target.
     */
    static Meas<DataType> GetEstMeas(const State& state, const MeasurementTypes type)  {
        return tDerived::DerivedGetEstMeas(state,type);
    } 

    /**
     * Performs the OMinus operation between two measurement (m1 ominus m2) of the same type. In other words, this
     * method computes the geodesic distance between two measurements of the same type.
     * @param[in] m1 a measurement
     * @param[in] m2 a measurement
     */
    static MatXd OMinus(const Meas<DataType>& m1, const Meas<DataType>& m2) {
        return tDerived::DerivedOMinus(m1, m2);
    } 


   /**
     * Generates a random measurement from a Gaussian distribution with mean defined by the state and standard deviation defined by meas_std. This
     * method is used primarily in simulations and tests.
     * @param[in] state    The state that serves as the mean of the Gaussian distribution.
     * @param[in] meas_std The measurement standard deviation.
     */ 
    Meas<DataType> GenerateRandomMeasurement(const State& state, const MatXd& meas_std){
        return static_cast<tDerived*>(this)->DerivedGenerateRandomMeasurement(state,meas_std);
    }

   /**
     * Generates a vector of random numbers from a Gaussian distribution of zero mean and unit standard deviation
     * @param[in] randn_nums The Gaussian random numbers to be generated
     */ 
    MatXd GaussianRandomGenerator(const int size);

    /**
     * Returns true if the state is inside the source's surveillance region. Note that the state is given in the global frame.  
     * @param[in] state A state of the target.
     */
    bool StateInsideSurveillanceRegion(const State& state) const {
        return state_in_surveillance_region_callback_(state);
    }

    /**
     * The Default callback function used with StateInsideSurveillanceRegion. It always returns true.  
     * @param[in] state A state of the target.
     */
    static bool StateInsideSurveillanceRegionDefaultCallback(const State& state) {
        return true;
    }

    /**
     * Calculates the temporal distance between two measurements.
     * @param[in] meas1 A measurement.
     * @param[in] meas2 A measurement.
     * @param[in] params The system parameters.
     * \return Returns temporal distance between two measurements
     */
   
    static DataType GetTemporalDistance(const Meas<DataType>& meas1, const Meas<DataType>& meas2, const Parameters& params) { return fabs(meas1.time_stamp - meas2.time_stamp); }

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
     * Verify that the parameters are valid. If they are, the parameters are set. 
     * @param[in] params Source parameters.
     * \return returns true if the parameters were set; otherwise, false.
     */
    bool SetParameters(const SourceParameters& params); 
    

// private:
     SourceBase();
    ~SourceBase();
    friend tDerived; /**< The object that inherits this base class. */

private:

    /**
     * Ensure that the source parameters meet the specified criteria. If a parameter doesn't, a runtime error will be thrown.
     * @param params The source parameters needed to initialize the source. 
     * \return Returns true if the parameters were successfully set. 
     */ 
    bool VerifySourceParameters(const SourceParameters& params);

    /**
     * This array of function pointers, holds pointers to the different methods of calculating the spatial distance
     * between two measurements. The index is determined by the measurement type. If two measurements are of different measurement
     * spaces, then their spatial distance cannot be calculated and a runtime error is thrown.
     */ 
    typedef double (*GSDFuncPTR)(const Meas<DataType>&, const Meas<DataType>&, const Parameters&);

    GSDFuncPTR **gsd_ptr_;

    static double GSD_RN_RN_POS(const Meas<DataType>& meas1, const Meas<DataType>& meas2, const Parameters& params) {return (meas1.pose - meas2.pose).norm();}
    static double GSD_SEN_SEN_POSE(const Meas<DataType>& meas1, const Meas<DataType>& meas2, const Parameters& params){return (State::Group::OMinus(meas1.pose,meas2.pose)).norm(); }
    static double GSD_SEN_SEN_POS(const Meas<DataType>& meas1, const Meas<DataType>& meas2, const Parameters& params){return (meas1.pose - meas2.pose).norm(); }
    static double GSD_NotImplemented(const Meas<DataType>& meas1, const Meas<DataType>& meas2, const Parameters& params){throw std::runtime_error("SourceBase::SpatialDistance Distance not implemented.");}

    

};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                            Definitions
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename tState, typename tDerived>
void SourceBase<tState,tDerived>::Init(const SourceParameters& params, std::function<bool(const State&)> state_in_surveillance_region_callback) {
    this->Init(params);
    this->state_in_surveillance_region_callback_ = state_in_surveillance_region_callback;
}

//-------------------------------------------------------------------------------

template<typename tState, typename tDerived>
void SourceBase<tState,tDerived>::Init(const SourceParameters& params) {

    bool success = SetParameters(params); // Verifies the parameters. If there is an invalid parameter, an error will be thrown. Otherwise, the parameters are set.
    this->state_in_surveillance_region_callback_ = StateInsideSurveillanceRegionDefaultCallback;
    static_cast<tDerived*>(this)->DerivedInit(params_);
}   

//-------------------------------------------------------------------------------

template<typename tState, typename tDerived>
Eigen::Matrix<typename tState::DataType,Eigen::Dynamic,Eigen::Dynamic>  SourceBase<tState,tDerived>::GaussianRandomGenerator(const int size){

    return utilities::GaussianRandomGenerator(size);
}

//-------------------------------------------------------------------------------

template<typename tState, typename tDerived>
SourceBase<tState,tDerived>::SourceBase() {

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



}

//---------------------------------------------------

template<typename tState, typename tDerived>
SourceBase<tState,tDerived>::~SourceBase() {

    for (int i = 0; i < MeasurementTypes::NUM_TYPES; i++) {
        delete [] gsd_ptr_[i];
    }
    delete [] gsd_ptr_;

}

//---------------------------------------------------

template<typename tState, typename tDerived>
bool SourceBase<tState,tDerived>::VerifySourceParameters(const SourceParameters& params) {

    bool success = true;
    unsigned int mult = 1;

    if(params.type_ == MeasurementTypes::RN_POS_VEL || params.type_ == MeasurementTypes::SEN_POS_VEL || params.type_ == MeasurementTypes::SEN_POSE_TWIST)
        mult = 2;


    if (params.meas_cov_.rows() != tDerived::meas_space_dim_*mult ) { // Make sure that it is not empty
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
    switch (params.type_)
    {
    case MeasurementTypes::RN_POS:
        break;    
    case MeasurementTypes::RN_POS_VEL:
        break;  
    case MeasurementTypes::SEN_POS:
        break;   
    case MeasurementTypes::SEN_POS_VEL:
        break;  
    case MeasurementTypes::SEN_POSE:
        break;  
    case MeasurementTypes::SEN_POSE_TWIST:
        break;  
    default:
        throw std::runtime_error("SourceBase::VerifySourceParameters: The measurement type is not known. ");
        success = false;

        break;
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

//---------------------------------------------------
template<typename tState, typename tDerived>
bool SourceBase<tState,tDerived>::SetParameters(const SourceParameters& params) {
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
    }
    return success;
}





} // namespace rransac



#endif // RRANSAC_COMMON_SOURCES_SOURCE_BASE_H_