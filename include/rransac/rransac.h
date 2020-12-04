#ifndef RRANSAC_RRANSAC_H_
#define RRANSAC_RRANSAC_H_

#include "system.h"
#include <functional>


namespace rransac {


/**
 * \class RRANSAC
 * This class is the interface for the user. 
 * 
 */ 

template<typename tModel>
class RRANSAC {

public:

    /**
     * This method creates a new source with the default surveillance region. Before adding a source,
     * it verifies that the source has a unique ID and that the member variables of the source
     * are set properly.
     * @param params The required source parameters
     * @return Returns true if the source was added.
     */
    bool AddSource(const SourceParameters& params);

    /**
     * This method creates a new source with a custom default surveillance region. Before adding a source,
     * it verifies that the source has a unique ID and that the member variables of the source
     * are set properly.
     * @param params The required source parameters
     * @param state_in_surveillance_region_callback The callback function used to verify that a track's state in the in surveillance region of a source.
     * @return Returns true if the source was added.
     */
    bool AddSource(const SourceParameters& params, std::function<bool(const State&)> state_in_surveillance_region_callback);

    /**
     * \detail Sets all of the parameters
     * \param[in] new_params The new parameters.
     * \return If the parameters were successfully set (i.e. met all of the prescibed requirements), it returns true. Otherwise false.
     */
    bool SetParameters(const Parameters &new_params) { return sys_.params_.SetParameters(new_params);  }


    /**
     * \detail Adds new measurements to RRANSAC. RRANSAC associates the measurements to tracks, clusters, or the data tree.
     * Tracks are updated with the new measurements.
     * @param new_measurements A vector that contains a vector of new measurements according to the source. That is, each vector of measurements corresponds to 
     * a unique source. 
     */
    void AddMeasurements(const std::vector<std::vector<Meas>>& new_measurements);

    /**
     * \detail Adds new measurements to RRANSAC and a new transformation. RRANSAC associates the measurements to tracks, clusters, or the data tree.
     * Tracks are propagated first, then transformed using the transformations. Measurements are associated with the track after the track has been propagated and then
     * transformed.
     * @param new_measurements A vector that contains a vector of new measurements according to the source. That is, each vector of measurements corresponds to 
     * @param transformation_data The data required to transform the tracks and measurements.
     * a unique source. 
     */
    void AddMeasurements(const std::vector<std::vector<Meas>>& new_measurements, const typename tModel::Transformation::Data& transformation_data);


    void run();

    const System<tModel>* GetSystemInformation() {return &sys_};

private:

    System<tModel> sys_;

};

}

#endif // RRANSAC_RRANSAC_H_