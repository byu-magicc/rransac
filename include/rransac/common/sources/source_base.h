#ifndef RRANSAC_COMMON_SOURCES_SOURCE_BASE_H_
#define RRANSAC_COMMON_SOURCES_SOURCE_BASE_H_

#include <Eigen/Core>
#include "state.h"
#include "common/measurement/measurement_base.h"
#include "parameters.h"

namespace rransac
{



// /** \enum SourceTypes
//  * Indicates the source types that are implemented. The nomenclature of the source type 
//  * indicates the manifold the target moves on, and the type of measurement.
//  * @see MeasurementTypes
//  */ 
// enum class SourceTypes {
//     R2_R2_POSE,            // The target moves on R2 and the measurement type is MeasurementTypes::R2_POSE
//     R2_R2_POSE_TWIST,      // The target moves on R2 and the measurement type is MeasurementTypes::R2_POSE_TWIST
//     SE2_R2_POS,
// }



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

    
    unsigned int source_index_;  /**< When a new source is added, it is added to the vector System::sources_. This is used to verify that the measurement corresponds to the proper source. */


    // SourceParameters(MeasurementTypes type) : type_(type) {}
    SourceParameters()=default;

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


template<class S>
class SourceBase
{

public:

    

     SourceBase();
    ~SourceBase();

    SourceParameters params_;  /** < The source parameters @see SourceParameters */

    /** Initializes the measurement source. This function must set the parameters.  */
    virtual void Init(const SourceParameters& params) {params_ = params;}       

    /** Returns the jacobian of the observation function w.r.t. the states */
    virtual Eigen::MatrixXd GetLinObsMatState(const S& state){return H_;}                              

    /** Returns the jacobian of the observation function w.r.t. the sensor noise */
    virtual Eigen::MatrixXd GetLinObsMatSensorNoise(const S& state){return V_;}                         

    /** Computes the estimated measurement given a state */
    virtual Eigen::MatrixXd GetEstMeas(const S& state){return state.g_.data_;} /** Returns an estimated measurement according to the state. */

    /**
     * Calculates the temporal distance between two measurements.
     * @param[in] meas1 A measurement.
     * @param[in] meas2 A different measurement.
     * @param[in] params Contains all of the user defined parameters. A user can define a weight when calculating the distances.
     * \return Returns temporal distance between two measurements
     */
   
    static double GetTemporalDistance(const Meas& meas1, const Meas& meas2, const Parameters& params) { return fabs(meas1.time_stamp - meas2.time_stamp); }

    /**
     * Calculates the spatial distance between two measurements depending on the type of measurement.
     * @param[in] meas1 A measurement.
     * @param[in] meas2 A different measurement.
     * @param[in] params Contains all of the user defined parameters. A user can define a weight when calculating the distances.
     * \return Returns spatial distance between two measurements
     */
    
    double GetSpatialDistance(const Meas& meas1, const Meas& meas2, const Parameters& params) {return gsd_ptr_[meas1.type][meas2.type](meas1,meas2,params);}


// protected:
public:
    Eigen::MatrixXd H_;
    Eigen::MatrixXd V_;

    typedef double (*GSDFuncPTR)(const Meas&, const Meas&, const Parameters&);

    GSDFuncPTR **gsd_ptr_;

    static double GSD_RN_RN_POS(const Meas& meas1, const Meas& meas2, const Parameters& params) {return (meas1.pose - meas2.pose).norm();}
    static double GSD_SEN_SEN_POSE(const Meas& meas1, const Meas& meas2, const Parameters& params){return (S::g_type_::OMinus(meas1.pose,meas2.pose)).norm(); }
    static double GSD_SEN_SEN_POS(const Meas& meas1, const Meas& meas2, const Parameters& params){return (meas1.pose - meas2.pose).norm(); }
    static double GSD_NotImplemented(const Meas& meas1, const Meas& meas2, const Parameters& params){throw std::runtime_error("SourceBase::SpatialDistance Distance not implemented.");}

};



template< class S>
SourceBase<S>::SourceBase() {

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

template< class S>
SourceBase<S>::~SourceBase() {

    for (int i = 0; i < MeasurementTypes::NUM_TYPES; i++) {
        delete [] gsd_ptr_[i];
    }
    delete [] gsd_ptr_;

}


} // namespace rransac



#endif // RRANSAC_COMMON_SOURCES_SOURCE_BASE_H_