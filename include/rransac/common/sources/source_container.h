#ifndef RRANSAC_COMMON_SOURCES_SOURCE_CONTAINER_H_
#define RRANSAC_COMMON_SOURCES_SOURCE_CONTAINER_H_
#pragma once


#include <vector>
#include "rransac/common/sources/source_base.h"
#include "rransac/common/sources/source_null.h"
#include "rransac/common/measurement/measurement_base.h"

namespace rransac {


template<typename... Ts> struct CountSources;

template<typename tSource, typename... Ts>
struct CountSources<tSource, Ts...> { static constexpr int value = std::is_same<tSource,SourceNull<>>::value ? 0 :1 + CountSources<Ts...>::value; };

template <>
struct CountSources<> { static constexpr int value = 0;};




template <typename S0 = SourceNull<>, typename S1 = SourceNull<>, typename S2 = SourceNull<>, typename S3 = SourceNull<>, typename S4=SourceNull<> >
class SourceContainer {

public:

typedef S0 Source0;
typedef S1 Source1;
typedef S2 Source2;
typedef S3 Source3;
typedef S4 Source4;

typedef typename S0::State State;                                           /**< The state of the target. @see State. */
typedef typename S0::DataType DataType;                                     /**< The scalar object for the data. Ex. float, double, etc. */
typedef Eigen::Matrix<DataType,Eigen::Dynamic,Eigen::Dynamic> MatXd;        /**< The object type of the Jacobians. */



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
const SourceParameters& GetParams(const unsigned int source_index);

/**
 * Changes some of the parameters of the source. 
 * \param[in] new_source_params The new source parameters.
 * \param[in] source_index The index to the source whose parameters are to be changed.
 * \return If the parameters were successfully set (i.e. met all of the prescibed requirements), it returns true. Otherwise false.
 */
bool ChangeSourceParameters(const unsigned int source_index, const SourceParameters &new_source_params);

/** 
 * Returns the jacobian of the observation function w.r.t. the states.
 * \param[in] source_index The index to the source whose parameters are to be changed.
 * @param[in] state A state of the target.
 * @param[in] transform_state A flag used to indicate if the state needs to be transformed 
 * @param[in] transform_data The data needed to transform the state
*/
static MatXd GetLinObsMatState(const unsigned int source_index, const State& state, const bool transform_state, const MatXd& transform_data);  

/** 
 * Returns the jacobian of the observation function w.r.t. the sensor noise.
 * \param[in] source_index The index to the source whose parameters are to be changed.
 * @param[in] state A state of the target.
 * @param[in] transform_state A flag used to indicate if the state needs to be transformed 
 * @param[in] transform_data The data needed to transform the state
 */
static MatXd GetLinObsMatSensorNoise(const unsigned int source_index, const State& state, const bool transform_state, const MatXd& transform_data); 

/**
 *  Implements the observation function and returns an estimated measurement based on the state. 
 * \param[in] source_index The index to the source whose parameters are to be changed.
 * @param[in] state A state of the target.
 * @param[in] transform_state A flag used to indicate if the state needs to be transformed 
 * @param[in] transform_data The data needed to transform the state
 */
static Meas<DataType> GetEstMeas(const unsigned int source_index, const State& state, const bool transform_state, const MatXd& transform_data);

/**
 * Performs the OMinus operation between two measurement (m1 ominus m2) of the same type. In other words, this
 * method computes the geodesic distance between two measurements of the same type.
 * @param[in] source_index The index to the source whose parameters are to be changed.
 * @param[in] m1 a measurement
 * @param[in] m2 a measurement
 */
static MatXd OMinus(const unsigned int source_index, const Meas<DataType>& m1, const Meas<DataType>& m2); 

/**
 * Generates a random measurement from a Gaussian distribution with mean defined by the state and standard deviation defined by meas_std. This
 * method is used primarily in simulations and tests.
 * @param[in] meas_std The measurement standard deviation.
 * @param[in] state    The state that serves as the mean of the Gaussian distribution.
 * @param[in] transform_state A flag used to indicate if the state needs to be transformed 
 * @param[in] transform_data The data needed to transform the state
 */ 
Meas<DataType> GenerateRandomMeasurement(const unsigned int source_index, const MatXd& meas_std, const State& state, const bool transform_state, const MatXd& transform_data) const;

/**
 * Returns true if the state is inside the source's surveillance region. Note that the state is given in the global frame.  
 * @param[in] source_index The index to the source whose parameters are to be changed.
 * @param[in] state A state of the target.
 * @param[in] transform_state A flag used to indicate if the state needs to be transformed 
 * @param[in] transform_data The data needed to transform the state
 */
bool StateInsideSurveillanceRegion(const unsigned int source_index, const State& state, const bool transform_state, const MatXd& transform_data) const ;

/**
 * Calculates the temporal distance between two measurements.
 * @param[in] meas1 A measurement.
 * @param[in] meas2 A measurement.
 * @param[in] params The system parameters.
 * \return Returns temporal distance between two measurements
 */
DataType GetTemporalDistance(const Meas<DataType>& meas1, const Meas<DataType>& meas2, const Parameters& params) const { return std::get<0>(sources_).GetTemporalDistance(meas1,meas2,params); }

/**
 * Calculates the geodesic distance between the pose of two measurements that have the same measurement space.
 * @param[in] meas1 A measurement.
 * @param[in] meas2 A measurement.
 * @param[in] params The system parameters.
 * \return Returns geodesic distance between pose of two measurements
 */
DataType GetSpatialDistance(const Meas<DataType>& meas1, const Meas<DataType>& meas2, const Parameters& params) const {return std::get<0>(sources_).GetSpatialDistance(meas1,meas2,params);}

/**
 * Finds the geodesic distance between the pose of two measurements of different time stamps normalized by the temporal distance. The measurements must have the same measurement space.
 * @param[in] meas1 A measurement.
 * @param[in] meas2 A measurement.
 * @param[in] params The system parameters.
 * \return Returns geodesic distance between two measurements
 */
DataType GetVelocityDistance(const Meas<DataType>& meas1, const Meas<DataType>& meas2, const Parameters& params) const { return std::get<0>(sources_).GetVelocityDistance(meas1,meas2,params);}


private:

static const std::size_t num_sources_ = CountSources<S0,S1,S2,S3,S4>::value;

std::tuple<S0,S1,S2,S3,S4> sources_;

std::vector<bool> source_initialized_;

/**
 * Verify that the parameters are valid. If they are, the parameters are set. 
 * @param[in] params Source parameters.
 * @param source_index The index to the source being inquired
 * \return returns true if the parameters were set; otherwise, false.
 */
bool SetParameters(const unsigned int source_index, const SourceParameters& params); 

};

//-------------------------------------------------------------------------------------------------------------------------------------------------
//                                                           Definitions
//-------------------------------------------------------------------------------------------------------------------------------------------------
template <typename S0, typename S1 , typename S2, typename S3 , typename S4 >
bool SourceContainer<S0,S1,S2,S3,S4>::AddSource(const SourceParameters& source_params) {

    bool added = false;

    if (source_params.source_index_ >= num_sources_) {
        throw std::runtime_error("SourceContainer::AddSource The source index must be less than the number of sources, i.e. in the range [0,num_sources).");
    } else if (source_initialized_[source_params.source_index_]) {
        throw std::runtime_error("SourceContainer::AddSource The source has already been initialized. ");
    } else {

        source_initialized_[source_params.source_index_] = true;

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

    if (source_params.source_index_ >= num_sources_) {
        throw std::runtime_error("SourceContainer::AddSource The source index must be less than the number of sources, i.e. in the range [0,num_sources).");
    } else if (source_initialized_[source_params.source_index_]) {
        throw std::runtime_error("SourceContainer::AddSource The source has already been initialized. ");
    } else {

        source_initialized_[source_params.source_index_] = true;

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
const SourceParameters& SourceContainer<S0,S1,S2,S3,S4>::GetParams(const unsigned int source_index) {


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
bool SourceContainer<S0,S1,S2,S3,S4>::SetParameters(const unsigned int source_index, const SourceParameters& params) {

    switch (source_index)
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
bool SourceContainer<S0,S1,S2,S3,S4>::ChangeSourceParameters(const unsigned int source_index, const SourceParameters &new_params) {
    if (new_params.source_index_ < 0 || new_params.source_index_ >= num_sources_) {
        throw std::runtime_error("SourceContainer::ChangeSourceParameters The source index must be greater than 0 and less than " + std::to_string(num_sources_));
        return false;
    }
    else if (GetParams(new_params.source_index_).source_index_ != new_params.source_index_) {
        throw std::runtime_error("SourceContainer::ChangeSourceParameters Cannot change the source index.");
        return false;
    } else {
        return SetParameters(new_params.source_index_, new_params);
    }
    
}

//-------------------------------------------------------------------------------

template <typename S0, typename S1 , typename S2, typename S3 , typename S4 >
Eigen::Matrix<typename S0::DataType,Eigen::Dynamic,Eigen::Dynamic> SourceContainer<S0,S1,S2,S3,S4>::GetLinObsMatState(const unsigned int source_index, const State& state, const bool transform_state, const MatXd& transform_data) {

    switch (source_index)
    {
    case 0:
        return S0::GetLinObsMatState(state, transform_state, transform_data);
        break;
    case 1:
        return S1::GetLinObsMatState(state, transform_state, transform_data);
        break;
    case 2:
        return S2::GetLinObsMatState(state, transform_state, transform_data);
        break;
    case 3:
        return S3::GetLinObsMatState(state, transform_state, transform_data);
        break;
    case 4:
        return S4::GetLinObsMatState(state, transform_state, transform_data);
        break;      
    default:
        throw std::runtime_error("SourceContainer::GetLinObsMatState The source index must be greater than 0 and less than " + std::to_string(num_sources_));
        break;
    }

}

//-------------------------------------------------------------------------------

template <typename S0, typename S1 , typename S2, typename S3 , typename S4 >
Eigen::Matrix<typename S0::DataType,Eigen::Dynamic,Eigen::Dynamic> SourceContainer<S0,S1,S2,S3,S4>::GetLinObsMatSensorNoise(const unsigned int source_index, const State& state, const bool transform_state, const MatXd& transform_data) {
    switch (source_index)
    {
    case 0:
        return S0::GetLinObsMatSensorNoise(state, transform_state, transform_data);
        break;
    case 1:
        return S1::GetLinObsMatSensorNoise(state, transform_state, transform_data);
        break;
    case 2:
        return S2::GetLinObsMatSensorNoise(state, transform_state, transform_data);
        break;
    case 3:
        return S3::GetLinObsMatSensorNoise(state, transform_state, transform_data);
        break;
    case 4:
        return S4::GetLinObsMatSensorNoise(state, transform_state, transform_data);
        break;      
    default:
        throw std::runtime_error("SourceContainer::GetLinObsMatSensorNoise The source index must be greater than 0 and less than " + std::to_string(num_sources_));
        break;
    }
}

//-------------------------------------------------------------------------------

template <typename S0, typename S1 , typename S2, typename S3 , typename S4 >
Meas<typename S0::DataType> SourceContainer<S0,S1,S2,S3,S4>::GetEstMeas(const unsigned int source_index, const State& state, const bool transform_state, const MatXd& transform_data) {
    switch (source_index)
    {
    case 0:
        return S0::GetEstMeas(state, transform_state, transform_data);
        break;
    case 1:
        return S1::GetEstMeas(state, transform_state, transform_data);
        break;
    case 2:
        return S2::GetEstMeas(state, transform_state, transform_data);
        break;
    case 3:
        return S3::GetEstMeas(state, transform_state, transform_data);
        break;
    case 4:
        return S4::GetEstMeas(state, transform_state, transform_data);
        break;      
    default:
        throw std::runtime_error("SourceContainer::GetEstMeas The source index must be greater than 0 and less than " + std::to_string(num_sources_));
        break;
    }
}

//-------------------------------------------------------------------------------

template <typename S0, typename S1 , typename S2, typename S3 , typename S4 >
Eigen::Matrix<typename S0::DataType,Eigen::Dynamic,Eigen::Dynamic> SourceContainer<S0,S1,S2,S3,S4>::OMinus(const unsigned int source_index, const Meas<DataType>& m1, const Meas<DataType>& m2) {
    switch (source_index)
    {
    case 0:
        return S0::OMinus(m1,m2);
        break;
    case 1:
        return S1::OMinus(m1,m2);
        break;
    case 2:
        return S2::OMinus(m1,m2);
        break;
    case 3:
        return S3::OMinus(m1,m2);
        break;
    case 4:
        return S4::OMinus(m1,m2);
        break;      
    default:
        throw std::runtime_error("SourceContainer::OMinus The source index must be greater than 0 and less than " + std::to_string(num_sources_));
        break;
    }
}

//-------------------------------------------------------------------------------
template <typename S0, typename S1 , typename S2, typename S3 , typename S4 >
Meas<typename S0::DataType> SourceContainer<S0,S1,S2,S3,S4>::GenerateRandomMeasurement(const unsigned int source_index, const MatXd& meas_std, const State& state, const bool transform_state, const MatXd& transform_data) const {
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
bool SourceContainer<S0,S1,S2,S3,S4>::StateInsideSurveillanceRegion(const unsigned int source_index, const State& state, const bool transform_state, const MatXd& transform_data) const {
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

// #ifdef DEBUG_BUILD
//         std::cerr << "RRANSAC::AddSource Cannot add sources once measuremens have been added."
// #endif

} // rransac

#endif // RRANSAC_COMMON_SOURCES_SOURCE_CONTAINER_H_