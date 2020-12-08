#ifndef RRANSAC_COMMON_SOURCES_SOURCE_BASE_H_
#define RRANSAC_COMMON_SOURCES_SOURCE_BASE_H_

#include <Eigen/Core>
#include <cmath>
#include <typeinfo>
#include <random>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <functional>


#include "parameters.h"
#include "common/measurement/measurement_base.h"
#include "state.h"
#include "common/utilities.h"


namespace rransac
{


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

    MeasurementTypes type_;          /** < The measurement type @see MeasurementTypes */

    float probability_of_detection_; /**< The probability that the phenomenon of interest is detected by a source during
                                      a single scan. This value must be between 0 and 1.*/

    float gate_probability_;         /**< The probability that a true measurement will be inside the validation region of a track. This value must
                                          be between 0 and 1. */

    
    unsigned int source_index_;  /**< When a new source is added, it is added to the vector System::sources_. This is used to verify that the measurement corresponds to the proper source. */

    // These parameters are not defined by the user, but are calculated depending on the user specified parameters.
    bool has_twist;                      /**< Indicates if the measurement has velocity data in addition to position */
    double gate_threshold_;              /**< The gate threshold of the validation region */
    double gate_threshold_sqrt_;         /**< The square root of the gate threshold */    
    double vol_unit_hypershpere_;        /**< The Volume of the unit hypershpere */

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


template<typename tState, typename tDerived>
class SourceBase
{

public:

    typedef tState State;
    static constexpr unsigned int meas_dim_ = tDerived::meas_dim_;               /**< The Dimension of the measurement */
    std::function<bool(const State&)> state_in_surveillance_region_callback_;

    SourceParameters params_;  /** < The source parameters @see SourceParameters */

    SourceBase(const SourceBase& other) : SourceBase() {
        params_ = other.params_;
        H_ = other.H_;
        V_ = other.V_;
        state_in_surveillance_region_callback_ = other.state_in_surveillance_region_callback_;
    }

    /** Initializes the measurement source. This function must set the parameters.  */
    void Init(const SourceParameters& params, std::function<bool(const State&)> state_in_surveillance_region_callback); 

    /** Initializes the measurement source. This function must set the parameters.  */
    void Init(const SourceParameters& params);     


    /** Returns the jacobian of the observation function w.r.t. the states */
    Eigen::MatrixXd GetLinObsMatState(const State& state){
       
        // std::cout << typeid(M).name() << std::endl;
        return static_cast<tDerived*>(this)->DerivedGetLinObsMatState(state);
    }                              

    /** Returns the jacobian of the observation function w.r.t. the sensor noise */
    Eigen::MatrixXd GetLinObsMatSensorNoise(const State& state){
        return static_cast<tDerived*>(this)->DerivedGetLinObsMatSensorNoise(state);
    }                         

    /** Computes the estimated measurement given a state */
    Meas GetEstMeas(const State& state){
        return static_cast<tDerived*>(this)->DerivedGetEstMeas(state);
    } 

    /**
     * Performs the OMinus operation between two measurement (m1 ominus m2)
     * @param m1 a measurement
     * @param m2 a measurement
     */
    Eigen::MatrixXd OMinus(const Meas& m1, const Meas& m2) {
        return static_cast<tDerived*>(this)->DerivedOMinus(m1, m2);
    } 

    /**
     * Maps the pose to Euclidean space. If the pose is already Euclidean space, then it returns the pose; otherwise
     * it returns the mapped pose.
     * @param Meas The measurement whose pose needs to be transformed
     */
    Eigen::MatrixXd ToEuclidean(const Meas& m)  {
        return static_cast<tDerived*>(this)->DerivedToEuclidean(m); 
    }

   /**
     * Generates a random measurement from a Gaussian distribution with mean defined by the state and covariance defined by meas_cov
     * @param state The state that serves as the mean in the Gaussian distribution
     * @param meas_std The measurement covariance
     */ 
    Meas GenerateRandomMeasurement(const State& state, const Eigen::MatrixXd& meas_std){
        return static_cast<tDerived*>(this)->DerivedGenerateRandomMeasurement(state,meas_std);
    }

   /**
     * Generates a vector of random numbers from a Gaussian distribution of zero mean and 1 standard deviation
     * @param randn_nums The Gaussian random numbers to be generated
     */ 
    Eigen::MatrixXd GaussianRandomGenerator(const int size);

    /**
     * Returns true if the state is inside the source's surveillance region. Note that the state is given in the global frame.  
     */
    bool StateInsideSurveillanceRegion(const State& state) {
        return state_in_surveillance_region_callback_(state);
    }

    /**
     * The Default callback function used with StateInsideSurveillanceRegion. It always returns true.  
     */
    static bool StateInsideSurveillanceRegionDefaultCallback(const State& state) {
        return true;
    }

    /**
     * Calculates the temporal distance between two measurements.
     * @param[in] meas1 A measurement.
     * @param[in] meas2 A different measurement.
     * @param[in] params Contains all of the user defined parameters. A user can define a weight when calculating the distances.
     * \return Returns temporal distance between two measurements
     */
   
    static double GetTemporalDistance(const Meas& meas1, const Meas& meas2, const Parameters& params) { return fabs(meas1.time_stamp - meas2.time_stamp); }

    /**
     * Calculates the geodesic distance between two measurements depending on the type of measurement.
     * @param[in] meas1 A measurement.
     * @param[in] meas2 A different measurement.
     * @param[in] params Contains all of the user defined parameters. A user can define a weight when calculating the distances.
     * \return Returns spatial distance between two measurements
     */
    
    double GetSpatialDistance(const Meas& meas1, const Meas& meas2, const Parameters& params) const {return gsd_ptr_[meas1.type][meas2.type](meas1,meas2,params);}

    /**
     * Finds the geodesic distance between two measurements of different time stamps and normalizes by the temproal distance.
     * @param[in] meas1 A measurement.
     * @param[in] meas2 A different measurement.
     * @param[in] params Contains all of the user defined parameters. A user can define a weight when calculating the distances.
     * \return Returns spatial distance between two measurements
     */
    static double GetVelocityDistance(const Meas& meas1, const Meas& meas2, const Parameters& params) {
        if (meas1.time_stamp == meas2.time_stamp) {
            throw std::runtime_error("SourceBase::GetVelocityDistance Measurements have the same time stamp");
        } else {
            return GetSpatialDistance(meas1,meas2,params)/GetTemporalDistance(meas1,meas2,params);
        }
    }

// protected:
    Eigen::MatrixXd H_;
    Eigen::MatrixXd V_;
    

// private:
     SourceBase();
    ~SourceBase();
    friend tDerived;

private:

    /**
     * Ensure that the source parameters meet the specified criteria. If a parameter doesn't, and error will be thrown.
     * @param params The source parameters needed to initialize the source. 
     */ 
    void VerifySourceParameters(const SourceParameters& params);

    typedef double (*GSDFuncPTR)(const Meas&, const Meas&, const Parameters&);

    GSDFuncPTR **gsd_ptr_;

    static double GSD_RN_RN_POS(const Meas& meas1, const Meas& meas2, const Parameters& params) {return (meas1.pose - meas2.pose).norm();}
    static double GSD_SEN_SEN_POSE(const Meas& meas1, const Meas& meas2, const Parameters& params){return (State::g_type_::OMinus(meas1.pose,meas2.pose)).norm(); }
    static double GSD_SEN_SEN_POS(const Meas& meas1, const Meas& meas2, const Parameters& params){return (meas1.pose - meas2.pose).norm(); }
    static double GSD_NotImplemented(const Meas& meas1, const Meas& meas2, const Parameters& params){throw std::runtime_error("SourceBase::SpatialDistance Distance not implemented.");}

    

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

    VerifySourceParameters(params); // Verifies the parameters. If there is an invalid parameter, an error will be thrown.
    this->state_in_surveillance_region_callback_ = StateInsideSurveillanceRegionDefaultCallback;
    
    this->params_ = params;
    if(this->params_.type_ == MeasurementTypes::RN_POS_VEL || this->params_.type_ == MeasurementTypes::SEN_POS_VEL || this->params_.type_ == MeasurementTypes::SEN_POSE_TWIST) {
        this->params_.has_twist = true;
        boost::math::chi_squared dist(meas_dim_*2);
        this->params_.gate_threshold_ = boost::math::quantile(dist, this->params_.gate_probability_);
        this->params_.vol_unit_hypershpere_ = pow(M_PI, static_cast<double>(meas_dim_))/boost::math::tgamma(static_cast<double>(meas_dim_) +1.0);
    } else {
        this->params_.has_twist = false;
        boost::math::chi_squared dist(meas_dim_);
        this->params_.gate_threshold_ = boost::math::quantile(dist, this->params_.gate_probability_);
        this->params_.vol_unit_hypershpere_ = pow(M_PI, static_cast<double>(meas_dim_)/2.0)/boost::math::tgamma(static_cast<double>(meas_dim_)/2.0 +1.0);
    }
    
    this->params_.gate_threshold_sqrt_ = sqrt(this->params_.gate_threshold_ ); 
    

    static_cast<tDerived*>(this)->DerivedInit(params_);
}   

//-------------------------------------------------------------------------------

template<typename tState, typename tDerived>
Eigen::MatrixXd  SourceBase<tState, tDerived>::GaussianRandomGenerator(const int size){

    return utilities::GaussianRandomGenerator(size);
}

//-------------------------------------------------------------------------------

template<typename tState, typename tDerived>
SourceBase<tState, tDerived>::SourceBase() {

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
SourceBase<tState, tDerived>::~SourceBase() {

    // std::cerr << "Destructing" << std::endl;

    for (int i = 0; i < MeasurementTypes::NUM_TYPES; i++) {
        delete [] gsd_ptr_[i];
    }
    delete [] gsd_ptr_;

}

//---------------------------------------------------

template<typename tState, typename tDerived>
void SourceBase<tState, tDerived>::VerifySourceParameters(const SourceParameters& params) {

    // Measurement covariance
    if (params.meas_cov_fixed_) { // Make sure that the measurement covariance is set properly.

        if (params.meas_cov_.rows() ==0 ) { // Make sure that it is not empty
            throw std::runtime_error("SourceBase::VerifySourceParameters: Measurement covariance cannot be empty if the measurement covariance is fixed.");
        } 

        if (  ((params.meas_cov_ + params.meas_cov_.transpose())/2.0 - params.meas_cov_).norm() > 1e-12 ) {
            throw std::runtime_error("SourceBase::VerifySourceParameters: Measurement covariance is not symmetic. ");
        }

        Eigen::VectorXcd eigen_values = params.meas_cov_.eigenvalues();
        for (int ii =0; ii < eigen_values.rows(); ++ii){                       // positive definite
            if(std::real(eigen_values(ii)) <0) {
                throw std::runtime_error("SourceBase::VerifySourceParameters: Measurement covariance is not positive definite. ");
            }
        }

    }

    // Expected number of false measurements
    if (params.expected_num_false_meas_ < 0 || params.expected_num_false_meas_ > 1) {
        throw std::runtime_error("SourceBase::VerifySourceParameters: The expected number of false measurements must be between 0 and 1. ");
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
        break;
    }

    // Verify the probability of detection 
    if (params.probability_of_detection_ < 0 || params.probability_of_detection_ > 1) {
        throw std::runtime_error("SourceBase::VerifySourceParameters: probability_of_detection_ must be between 0 and 1. ");
    }

    // Verify the gate probability
    if (params.gate_probability_ < 0 || params.gate_probability_ > 1) {
        throw std::runtime_error("SourceBase::VerifySourceParameters: gate_probability_ must be between 0 and 1. ");
    }

}





} // namespace rransac



#endif // RRANSAC_COMMON_SOURCES_SOURCE_BASE_H_