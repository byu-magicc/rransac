#ifndef RRANSAC_COMMON_SOURCES_SOURCE_CONTAINER_H_
#define RRANSAC_COMMON_SOURCES_SOURCE_CONTAINER_H_
#pragma once


#include <vector>
#include "rransac/common/sources/source_base.h"
#include "rransac/common/sources/source_null.h"
#include "rransac/common/measurement/measurement_base.h"
#include "rransac/common/utilities.h"

namespace rransac {


template<typename... Ts> struct CountSources;

template<typename _Source, typename... Ts>
struct CountSources<_Source, Ts...> { static constexpr int value = std::is_same<_Source,SourceNull<>>::value ? 0 :1 + CountSources<Ts...>::value; };

template <>
struct CountSources<> { static constexpr int value = 0;};




template <typename S0 = SourceNull<>, typename S1 = SourceNull<>, typename S2 = SourceNull<>, typename S3 = SourceNull<>, typename S4=SourceNull<> >
class SourceContainer {

public:

typedef S0 Source0;
typedef typename std::conditional<IsSourceNull<S1>::value,SourceNull<typename S0::State, S0::measurement_type_, S0::Base::DerivedTraits::template Transformation>,S1>::type Source1;
typedef typename std::conditional<IsSourceNull<S2>::value,SourceNull<typename S0::State, S0::measurement_type_, S0::Base::DerivedTraits::template Transformation>,S2>::type Source2;
typedef typename std::conditional<IsSourceNull<S3>::value,SourceNull<typename S0::State, S0::measurement_type_, S0::Base::DerivedTraits::template Transformation>,S3>::type Source3;
typedef typename std::conditional<IsSourceNull<S4>::value,SourceNull<typename S0::State, S0::measurement_type_, S0::Base::DerivedTraits::template Transformation>,S4>::type Source4;
// typedef S2 Source2;
// typedef S3 Source3;
// typedef S4 Source4;

typedef typename S0::State State;                                           /**< The state of the target. @see State. */
typedef typename S0::DataType DataType;                                     /**< The scalar object for the data. Ex. float, double, etc. */
typedef typename S0::Transformation Transformation;                         /**< The transformation data type. */
typedef typename S0::MatH MatH;                                             /**< The object type of the Jacobian of the observation function w.r.t. the states. */
typedef typename S0::MatV MatV;                                             /**< The object type of the Jacobian of the observation function w.r.t. the measurement noise. */
static constexpr int num_sources_ = CountSources<S0,S1,S2,S3,S4>::value;    /**< The number of user defined sources. */
typedef typename S0::ModelCompatibility ModelCompatibility;                 /**< Indicates which model this source is compatible with. */
typedef typename S0::TransformDataType TransformDataType;                   /**< The transform data type. */
typedef typename S0::Measurement Measurement;                               /**< The object type of the measurement. */
typedef Eigen::Matrix<DataType,Eigen::Dynamic,Eigen::Dynamic> MatXd;  
int sources_initialized_ = 0;       


// Ensure that all sources have the same model compatibility
static_assert( std::is_same<typename S0::ModelCompatibility, typename S1::ModelCompatibility>::value || IsSourceNull<S1>::value, "The sources are not compatible with the same model" );
static_assert( std::is_same<typename S0::ModelCompatibility, typename S2::ModelCompatibility>::value || IsSourceNull<S2>::value, "The sources are not compatible with the same model" );
static_assert( std::is_same<typename S0::ModelCompatibility, typename S3::ModelCompatibility>::value || IsSourceNull<S3>::value, "The sources are not compatible with the same model" );
static_assert( std::is_same<typename S0::ModelCompatibility, typename S4::ModelCompatibility>::value || IsSourceNull<S4>::value, "The sources are not compatible with the same model" );



// template <typename tScalar, template<typename> typename tStateTemplate>
// using ModelTemplate = ModelRN<tSourceContainer>; /**< Used to create a model of the state, source and transformation, but with a different DataType. This is needed to solve the 
//                                                                                      nonlinear log maximum likelihood estimation problem by Ceres. */

template<typename _DataType>
using SourceContainerTemplate = SourceContainer<typename Source0::template SourceTemplate<_DataType>, typename Source1::template SourceTemplate<_DataType>, typename Source2::template SourceTemplate<_DataType>, typename Source3::template SourceTemplate<_DataType>, typename Source4::template SourceTemplate<_DataType>>;



SourceContainer() : source_initialized_(num_sources_,false) {};

/**
 * This method creates a new source with the default surveillance region. Before adding a source,
 * it verifies that the source has a unique ID and that the member variables of the source
 * are set properly.
 * @param[in] source_params The required source parameters
 * @return Returns true if the source was added.
 */
bool AddSource(const SourceParameters& source_params);

/**
 * This method creates a new source with a custom default surveillance region. Before adding a source,
 * it verifies that the source has a unique ID and that the member variables of the source
 * are set properly.
 * @param[in] source_params The required source parameters
 * @param[in] state_in_surveillance_region_callback The callback function used to verify that a track's state in the in surveillance region of a source.
 * @return Returns true if the source was added.
 */
bool AddSource(const SourceParameters& source_params, std::function<bool(const State&)> state_in_surveillance_region_callback);


/**
 * Retruns a reference to the source parameters.
 * @param source_index The index to the source being inquired
 */ 
const SourceParameters& GetParams(const unsigned int source_index) const;

/**
 * Changes some of the parameters of the source. 
 * \param[in] new_source_params The new source parameters.
 * \param[in] source_index The index to the source whose parameters are to be changed.
 * \return If the parameters were successfully set (i.e. met all of the prescibed requirements), it returns true. Otherwise false.
 */
bool ChangeSourceParameters(const SourceParameters &new_source_params);

/** 
 * Returns the jacobian of the observation function w.r.t. the states.
 * \param[in] source_index The index to the source whose parameters are to be changed.
 * @param[in] state A state of the target.
 * @param[in] transform_state A flag used to indicate if the state needs to be transformed 
 * @param[in] transform_data The data needed to transform the state
*/
static MatH GetLinObsMatState(const unsigned int source_index, const State& state, const bool transform_state, const TransformDataType& transform_data);  

/** 
 * Returns the jacobian of the observation function w.r.t. the sensor noise.
 * \param[in] source_index The index to the source whose parameters are to be changed.
 * @param[in] state A state of the target.
 * @param[in] transform_state A flag used to indicate if the state needs to be transformed 
 * @param[in] transform_data The data needed to transform the state
 */
static MatV GetLinObsMatSensorNoise(const unsigned int source_index, const State& state, const bool transform_state, const TransformDataType& transform_data); 

/**
 *  Implements the observation function and returns an estimated measurement based on the state. 
 * \param[in] source_index The index to the source whose parameters are to be changed.
 * @param[in] state A state of the target.
 * @param[in] transform_state A flag used to indicate if the state needs to be transformed 
 * @param[in] transform_data The data needed to transform the state
 */
static Measurement GetEstMeas(const unsigned int source_index, const State& state, const bool transform_state, const TransformDataType& transform_data);

/**
 * Performs the OMinus operation between two measurement (m1 ominus m2) of the same type. In other words, this
 * method computes the geodesic distance between two measurements of the same type.
 * @param[in] source_index The index to the source whose parameters are to be changed.
 * @param[in] m1 a measurement
 * @param[in] m2 a measurement
 */
static MatXd OMinus(const unsigned int source_index, const Measurement& m1, const Measurement& m2); 

/**
 * Generates a random measurement from a Gaussian distribution with mean defined by the state and standard deviation defined by meas_std. This
 * method is used primarily in simulations and tests.
 * @param[in] meas_std The measurement standard deviation.
 * @param[in] state    The state that serves as the mean of the Gaussian distribution.
 * @param[in] transform_state A flag used to indicate if the state needs to be transformed 
 * @param[in] transform_data The data needed to transform the state
 */ 
Measurement GenerateRandomMeasurement(const unsigned int source_index, const MatXd& meas_std, const State& state, const bool transform_state, const TransformDataType& transform_data) const;

/**
 * Returns true if the state is inside the source's surveillance region. Note that the state is given in the global frame.  
 * @param[in] source_index The index to the source whose parameters are to be changed.
 * @param[in] state A state of the target.
 * @param[in] transform_state A flag used to indicate if the state needs to be transformed 
 * @param[in] transform_data The data needed to transform the state
 */
bool StateInsideSurveillanceRegion(const unsigned int source_index, const State& state, const bool transform_state, const TransformDataType& transform_data) const ;

/**
 * Calculates the temporal distance between two measurements.
 * @param[in] meas1 A measurement.
 * @param[in] meas2 A measurement.
 * @param[in] params The system parameters.
 * \return Returns temporal distance between two measurements
 */
DataType GetTemporalDistance(const Measurement& meas1, const Measurement& meas2, const Parameters& params) const { return std::get<0>(sources_).GetTemporalDistance(meas1,meas2,params); }

/**
 * Calculates the geodesic distance between the pose of two measurements that have the same measurement space.
 * @param[in] meas1 A measurement.
 * @param[in] meas2 A measurement.
 * @param[in] params The system parameters.
 * \return Returns geodesic distance between pose of two measurements
 */
DataType GetSpatialDistance(const Measurement& meas1, const Measurement& meas2, const Parameters& params) const {return std::get<0>(sources_).GetSpatialDistance(meas1,meas2,params);}

/**
 * Finds the geodesic distance between the pose of two measurements of different time stamps normalized by the temporal distance. The measurements must have the same measurement space.
 * @param[in] meas1 A measurement.
 * @param[in] meas2 A measurement.
 * @param[in] params The system parameters.
 * \return Returns geodesic distance between two measurements
 */
DataType GetVelocityDistance(const Measurement& meas1, const Measurement& meas2, const Parameters& params) const { return std::get<0>(sources_).GetVelocityDistance(meas1,meas2,params);}

/**
 * Verifies that the data in the measurement meets certain specifications.
 * @param measurement The measurement to be verified. 
 */ 
bool IsAcceptableMeasurement(const Measurement& measurement);


private:


std::tuple<Source0,Source1,Source2,Source3,Source4> sources_;

std::vector<bool> source_initialized_;

/**
 * Verify that the parameters are valid. If they are, the parameters are set. 
 * @param[in] params Source parameters.
 * @param source_index The index to the source being inquired
 * \return returns true if the parameters were set; otherwise, false.
 */
bool SetParameters(const SourceParameters& params); 

};

//-------------------------------------------------------------------------------------------------------------------------------------------------
//                                                           Definitions
//-------------------------------------------------------------------------------------------------------------------------------------------------
template <typename S0, typename S1 , typename S2, typename S3 , typename S4 >
bool SourceContainer<S0,S1,S2,S3,S4>::AddSource(const SourceParameters& source_params) {

    bool added = false;

    if (source_params.source_index_ >= num_sources_ || source_params.source_index_ < 0) {
        throw std::runtime_error("SourceContainer::AddSource The source index must be less than the number of sources, and greater than or equal to zero; i.e. in the range [0,num_sources).");
    } else if (source_initialized_[source_params.source_index_]) {
        throw std::runtime_error("SourceContainer::AddSource The source has already been initialized. ");
    } else {

        source_initialized_[source_params.source_index_] = true;
        sources_initialized_++;

        switch (source_params.source_index_)
        {
        case 0:
            std::get<0>(sources_).Init(source_params);
            added = true;
            break;
        case 1:
            std::get<1>(sources_).Init(source_params);
            added = true;
            break;
        case 2:
            std::get<2>(sources_).Init(source_params);
            added = true;
            break;
        case 3:
            std::get<3>(sources_).Init(source_params);
            added = true;
            break;
        case 4:
            std::get<4>(sources_).Init(source_params);
            added = true;
            break;      
        default:
            break;
        }
    }

    return added;

}

//-------------------------------------------------------------------------------


template <typename S0, typename S1 , typename S2, typename S3 , typename S4 >
bool SourceContainer<S0,S1,S2,S3,S4>::AddSource(const SourceParameters& source_params, std::function<bool(const State&)> state_in_surveillance_region_callback) {

    bool added = false;

    if (source_params.source_index_ >= num_sources_ || source_params.source_index_ < 0) {
        throw std::runtime_error("SourceContainer::AddSource The source index must be less than the number of sources, and greater than or equal to zero; i.e. in the range [0,num_sources).");
    } else if (source_initialized_[source_params.source_index_]) {
        throw std::runtime_error("SourceContainer::AddSource The source has already been initialized. ");
    } else {

        source_initialized_[source_params.source_index_] = true;
        sources_initialized_++;

        switch (source_params.source_index_)
        {
        case 0:
            std::get<0>(sources_).Init(source_params,state_in_surveillance_region_callback);
            added = true;
            break;
        case 1:
            std::get<1>(sources_).Init(source_params,state_in_surveillance_region_callback);
            added = true;
            break;
        case 2:
            std::get<2>(sources_).Init(source_params,state_in_surveillance_region_callback);
            added = true;
            break;
        case 3:
            std::get<3>(sources_).Init(source_params,state_in_surveillance_region_callback);
            added = true;
            break;
        case 4:
            std::get<4>(sources_).Init(source_params,state_in_surveillance_region_callback);
            added = true;
            break;      
        default:
            break;
        }
    }

    return added;

}

//-------------------------------------------------------------------------------
template <typename S0, typename S1 , typename S2, typename S3 , typename S4 >
const SourceParameters& SourceContainer<S0,S1,S2,S3,S4>::GetParams(const unsigned int source_index) const {


    switch (source_index)
    {
    case 0:
        return std::get<0>(sources_).GetParams();
        break;
    case 1:
        return std::get<1>(sources_).GetParams();
        break;
    case 2:
        return std::get<2>(sources_).GetParams();
        break;
    case 3:
        return std::get<3>(sources_).GetParams();
        break;
    case 4:
        return std::get<4>(sources_).GetParams();
        break;      
    default:
        throw std::runtime_error("SourceContainer::GetParams The source index must be greater than 0 and less than " + std::to_string(num_sources_)); 
        break;
    }



}

//-------------------------------------------------------------------------------
template <typename S0, typename S1 , typename S2, typename S3 , typename S4 >
bool SourceContainer<S0,S1,S2,S3,S4>::SetParameters(const SourceParameters& params) {

    switch (params.source_index_)
    {
    case 0:
        return std::get<0>(sources_).SetParameters(params);
        break;
    case 1:
        return std::get<1>(sources_).SetParameters(params);
        break;
    case 2:
        return std::get<2>(sources_).SetParameters(params);
        break;
    case 3:
        return std::get<3>(sources_).SetParameters(params);
        break;
    case 4:
        return std::get<4>(sources_).SetParameters(params);
        break;      
    default:
     throw std::runtime_error("SourceContainer::SetParameters The source index must be greater than 0 and less than " + std::to_string(num_sources_)); 
        break;
    }

} 

//-------------------------------------------------------------------------------

template <typename S0, typename S1 , typename S2, typename S3 , typename S4 >
bool SourceContainer<S0,S1,S2,S3,S4>::ChangeSourceParameters(const SourceParameters &new_params) {
    if (new_params.source_index_ < 0 || new_params.source_index_ >= num_sources_) {
        throw std::runtime_error("SourceContainer::ChangeSourceParameters The source index must be greater than 0 and less than " + std::to_string(num_sources_));
        return false;
    }
    else if (GetParams(new_params.source_index_).source_index_ != new_params.source_index_) {
        throw std::runtime_error("SourceContainer::ChangeSourceParameters Cannot change the source index.");
        return false;
    } else {
        return SetParameters(new_params);
    }
    
}

//-------------------------------------------------------------------------------

template <typename S0, typename S1 , typename S2, typename S3 , typename S4 >
typename SourceContainer<S0,S1,S2,S3,S4>::MatH SourceContainer<S0,S1,S2,S3,S4>::GetLinObsMatState(const unsigned int source_index, const State& state, const bool transform_state, const TransformDataType& transform_data) {

    switch (source_index)
    {
    case 0:
        return Source0::GetLinObsMatState(state, transform_state, transform_data);
        break;
    case 1:
        return Source1::GetLinObsMatState(state, transform_state, transform_data);
        break;
    case 2:
        return Source2::GetLinObsMatState(state, transform_state, transform_data);
        break;
    case 3:
        return Source3::GetLinObsMatState(state, transform_state, transform_data);
        break;
    case 4:
        return Source4::GetLinObsMatState(state, transform_state, transform_data);
        break;      
    default:
        throw std::runtime_error("SourceContainer::GetLinObsMatState The source index must be greater than 0 and less than " + std::to_string(num_sources_));
        break;
    }

}

//-------------------------------------------------------------------------------

template <typename S0, typename S1 , typename S2, typename S3 , typename S4 >
typename SourceContainer<S0,S1,S2,S3,S4>::MatV SourceContainer<S0,S1,S2,S3,S4>::GetLinObsMatSensorNoise(const unsigned int source_index, const State& state, const bool transform_state, const TransformDataType& transform_data) {
    switch (source_index)
    {
    case 0:
        return Source0::GetLinObsMatSensorNoise(state, transform_state, transform_data);
        break;
    case 1:
        return Source1::GetLinObsMatSensorNoise(state, transform_state, transform_data);
        break;
    case 2:
        return Source2::GetLinObsMatSensorNoise(state, transform_state, transform_data);
        break;
    case 3:
        return Source3::GetLinObsMatSensorNoise(state, transform_state, transform_data);
        break;
    case 4:
        return Source4::GetLinObsMatSensorNoise(state, transform_state, transform_data);
        break;      
    default:
        throw std::runtime_error("SourceContainer::GetLinObsMatSensorNoise The source index must be greater than 0 and less than " + std::to_string(num_sources_));
        break;
    }
}

//-------------------------------------------------------------------------------

template <typename S0, typename S1 , typename S2, typename S3 , typename S4 >
typename SourceContainer<S0,S1,S2,S3,S4>::Measurement SourceContainer<S0,S1,S2,S3,S4>::GetEstMeas(const unsigned int source_index, const State& state, const bool transform_state, const TransformDataType& transform_data) {
    switch (source_index)
    {
    case 0:
        return Source0::GetEstMeas(state, transform_state, transform_data);
        break;
    case 1:
        return Source1::GetEstMeas(state, transform_state, transform_data);
        break;
    case 2:
        return Source2::GetEstMeas(state, transform_state, transform_data);
        break;
    case 3:
        return Source3::GetEstMeas(state, transform_state, transform_data);
        break;
    case 4:
        return Source4::GetEstMeas(state, transform_state, transform_data);
        break;      
    default:
        throw std::runtime_error("SourceContainer::GetEstMeas The source index must be greater than 0 and less than " + std::to_string(num_sources_));
        break;
    }
}

//-------------------------------------------------------------------------------

template <typename S0, typename S1 , typename S2, typename S3 , typename S4 >
typename SourceContainer<S0,S1,S2,S3,S4>::MatXd SourceContainer<S0,S1,S2,S3,S4>::OMinus(const unsigned int source_index, const Measurement& m1, const Measurement& m2) {
    switch (source_index)
    {
    case 0:
        return Source0::OMinus(m1,m2);
        break;
    case 1:
        return Source1::OMinus(m1,m2);
        break;
    case 2:
        return Source2::OMinus(m1,m2);
        break;
    case 3:
        return Source3::OMinus(m1,m2);
        break;
    case 4:
        return Source4::OMinus(m1,m2);
        break;      
    default:
        throw std::runtime_error("SourceContainer::OMinus The source index must be greater than 0 and less than " + std::to_string(num_sources_));
        break;
    }
}

//-------------------------------------------------------------------------------
template <typename S0, typename S1 , typename S2, typename S3 , typename S4 >
typename SourceContainer<S0,S1,S2,S3,S4>::Measurement SourceContainer<S0,S1,S2,S3,S4>::GenerateRandomMeasurement(const unsigned int source_index, const MatXd& meas_std, const State& state, const bool transform_state, const TransformDataType& transform_data) const {
    switch (source_index)
    {
    case 0:
        return std::get<0>(sources_).GenerateRandomMeasurement(meas_std, state, transform_state, transform_data);
        break;
    case 1:
        return std::get<1>(sources_).GenerateRandomMeasurement(meas_std, state, transform_state, transform_data);
        break;
    case 2:
        return std::get<2>(sources_).GenerateRandomMeasurement(meas_std, state, transform_state, transform_data);
        break;
    case 3:
        return std::get<3>(sources_).GenerateRandomMeasurement(meas_std, state, transform_state, transform_data);
        break;
    case 4:
        return std::get<4>(sources_).GenerateRandomMeasurement(meas_std, state, transform_state, transform_data);
        break;      
    default:
        throw std::runtime_error("SourceContainer::GenerateRandomMeasurement The source index must be greater than 0 and less than " + std::to_string(num_sources_));
        break;
    }
}

//-------------------------------------------------------------------------------

template <typename S0, typename S1 , typename S2, typename S3 , typename S4 >
bool SourceContainer<S0,S1,S2,S3,S4>::StateInsideSurveillanceRegion(const unsigned int source_index, const State& state, const bool transform_state, const TransformDataType& transform_data) const {
    switch (source_index)
    {
    case 0:
        return std::get<0>(sources_).StateInsideSurveillanceRegion(state, transform_state, transform_data);
        break;
    case 1:
        return std::get<1>(sources_).StateInsideSurveillanceRegion(state, transform_state, transform_data);
        break;
    case 2:
        return std::get<2>(sources_).StateInsideSurveillanceRegion(state, transform_state, transform_data);
        break;
    case 3:
        return std::get<3>(sources_).StateInsideSurveillanceRegion(state, transform_state, transform_data);
        break;
    case 4:
        return std::get<4>(sources_).StateInsideSurveillanceRegion(state, transform_state, transform_data);
        break;      
    default:
        throw std::runtime_error("SourceContainer::GenerateRandomMeasurement The source index must be greater than 0 and less than " + std::to_string(num_sources_));
        break;
    }
}

//-------------------------------------------------------------------------------

template <typename S0, typename S1 , typename S2, typename S3 , typename S4 >
bool SourceContainer<S0,S1,S2,S3,S4>::IsAcceptableMeasurement(const Measurement& measurement) {
    switch (measurement.source_index)
    {
    case 0:
        return std::get<0>(sources_).IsAcceptableMeasurement(measurement);
        break;
    case 1:
        return std::get<1>(sources_).IsAcceptableMeasurement(measurement);
        break;
    case 2:
        return std::get<2>(sources_).IsAcceptableMeasurement(measurement);
        break;
    case 3:
        return std::get<3>(sources_).IsAcceptableMeasurement(measurement);
        break;
    case 4:
        return std::get<4>(sources_).IsAcceptableMeasurement(measurement);
        break;      
    default:
        throw std::runtime_error("SourceContainer::IsAcceptableMeasurement The source index must be greater than 0 and less than " + std::to_string(num_sources_));
        break;
    }
}


} // rransac

#endif // RRANSAC_COMMON_SOURCES_SOURCE_CONTAINER_H_