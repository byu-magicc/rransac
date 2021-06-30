#ifndef RRANSAC_RRANSAC_H_
#define RRANSAC_RRANSAC_H_
#pragma once

#include <functional>
#include <string>
#include <iostream>

#include "rransac/system.h"
#include "rransac/track_initialization/ransac.h"
#include "rransac/common/data_association/data_association_host.h"
#include "rransac/common/models/model_manager.h"
#include "rransac/common/sources/source_null.h"




namespace rransac {

/**
 * \class RRANSACTemplateParameters
 * Contains all of the template parameters for RRANSAC that used to initialize RRANSAC. See RRANSAC for more information
 */ 
template<typename tSourceContainer, template<typename> typename tModel, template <typename > typename tSeed, template<typename , template <typename > typename > typename tLMLEPolicy, template<class> typename tValidationRegionPolicy, template<class> typename tUpdateTrackLikelihoodPolicy, template<class> typename tMeasurementWeightPolicy>
struct RRANSACTemplateParameters {

    
    typedef tSourceContainer TSourceContainer;                                                       /**< The container of all of the sources. */
    typedef typename tSourceContainer::State TState;                                                 /**< The state of the target. @see State. */
    typedef typename TState::DataType TDataType;                                                     /**< The scalar object for the data. Ex. float, double, etc. */ 
    typedef typename tSourceContainer::Transformation TTransformation;                               /**< The object type of the measurement and track transformation. @see TransformBase*/
    typedef typename TTransformation::Data TTransformationData;                                      /**< The object type of the transformation data. @see TransformBase*/
    typedef tModel<TSourceContainer> TModel;                                                         /**< The object type of the model. @see ModelBase. */                          
    typedef ModelManager<TModel> TModelManager;                                                      /**< The object type of the model manager. @see ModelManager. */
    typedef Meas<TDataType> TMeas;                                                                   /**< The object type of the measurement. @see Meas */
    typedef Ransac<TModel,tSeed,tLMLEPolicy,tValidationRegionPolicy,tUpdateTrackLikelihoodPolicy,tMeasurementWeightPolicy> TRansac;               /**< The object type of Ransac or the track initializer. @see Ransac */
    typedef DataAssociationHost<TModel,tValidationRegionPolicy,tUpdateTrackLikelihoodPolicy,tMeasurementWeightPolicy> TDataAssociation;           /**< The object type of the data association method. @see DataAssociationHost */


};

/**
 * \class RRANSAC
 * This class serves as the user interface. The program that implements RRANSAC is designed to be very modular: to work with different models, measurement
 * sources, data association policies, etc. So the first step is to properly set up RRANSAC. The template parameters for RRANSAC are specified by the class RRANSACTemplateParameters.
 * It requires the type of state, source, transformaiton, model, model data association policy, cluster and data tree association policy, seed for the log maximum likelihood estimation 
 * optimization, and the log maximum likelihood estimation. For those familiar with RRANSAC, we use a policy method where each template parameter is a policy.
 * 
 * The type of state is the configuration manifold of the target which is a Lie group. The type of source is the measurement source that models observation function of the system model.
 * The source must be compatible with the state type. The transformation type specifies how to transform measurements and tracks in the case that the tracking frame changes. The type of model
 * describes the full system model: propagate, update, etc. The model must be compatible with the type of state and source. The type of data association policy determines how new measurements
 * will associated to existing tracks. There is only one current implementation which is the probabilistic data association filter. This policy is compatible with all models. The cluster
 * and data tree association policy determins how new measurements that were not associated to an existing track are to be associated with clusters and the data tree. There is only currently 
 * one implementation of this. The seed type determines how the nonlinear optimizer for the LMLE optimization will be seeded. Lastly, the log maximum likelihood estimation type determines 
 * what type of LMLE solver to use. Currently there are two types: one for linear systems and one for non linear systems. There are several compatibility checks but they are not comprehensive
 * so you must be carefull with putting the pieces together. 
 * 
 * All of the different pieces can be overwhelming at first but currently the few choices make it fairly simple. Since there is currently only one update track likelihood policy and one measurement
 * weight policy use TLI_IPDAFPolicy and MW_IPDAFPolicy. There are several choices for the validation region. See the description of each one to help you choose the one you want. 
 * There are only two LMLE solver policies: one for linear systems and one
 * for nonlinear systems. If your system is linear, use LinearLMLEPolicy; otherwise, use NonLinearLMLEPolicy. There are a few more choices with the seed policy; however, if your system is linear or
 * you want to initialize the optimizer with zeros, use the NULLSeedPolicy; otherwise, look at the other options.  * 
 * The transformation policy, transforms the measurements and tracks from the previous tracking frame to the current tracking frame. If the tracking frame doesn't move, then use TransformNULL; otherwise
 * look at the different options or make your own. 
 * 
 * There are several different states available. See State for your options. There are several measurement sources available. If the measurement space is Euclidean in cartesian coordinates and
 * the target's configuration manifold is RN, use SourceRN. If the target's configuration manifold is SE2 and only the position and its derivative are observable, use SourceSENPosVel. If the measurement space
 * is SEN and the target's configuration manifold is SEN, then use SourceSENPoseTwist. * 
 * There are several options for the model. If your target's configuration manifold is RN, use ModelRN. If it is SEN and the corresponding measurement source is SourceSENPosVel, use ModelSENPosVel. Lastly, if
 * your target's configuration manifold is SEN and the measurement space is SEN, use ModelSENPoseTwist. 
 * 
 * 
 * For example, suppose that the target's configuration manifold is R2 and the measurement space is R2 in cartesian coordinates. Also suppose that the tracking frame is stationary then RRANSAC template parameters are
 * RRANSACTemplateParameters<lie_groups::R2_r2,SourceRN,TransformNULL,ModelRN,NULLSeedPolicy,LinearLMLEPolicy,ModelPDFPolicy,DataTreeClusterAssociationPolicy>
 * 
 * Another example. Suppose that the target's configuration manifold is SE2, the measurement source is a camera with only position and velocity being observe, the tracking frame is the camera frame and thus changes,
 * but it's transformation can be computed using the homography, then the rRANSAC template parameters are
 * RRANSACTemplateParameters<lie_groups::SE2_se2,SourceSENPosVel,TransformHomography,ModelSENPosVel,SE2PosSeedPolicy,NonLinearLMLEPolicy,ModelPDFPolicy,DataTreeClusterAssociationPolicy>
 * 
 * If there are questions, please ask. 
 * 
 * 
 * Once the policies (template parameters are selected) we need to setup RRANSAC by adding system parameters and measurement sources. Measurement sources are added using the member functions AddSource(const SourceParameters& params) 
 * and AddSource(const SourceParameters& params, std::function<bool(const State_&)> state_in_surveillance_region_callback). See SourceParameters for information regarding the parameters. Measurement sources must be added
 * in order according to their source index starting from 0 and incrementing by 1. Once a source is added, it cannot be removed; however, their parameters can be changed by ChangeSourceParameters. The system parameters
 * can be set by SetSystemParameters, and can be changed at any time calling the same method. 
 * 
 * New measurements are adding using the member methods AddMeasurements(const std::list<Meas_>& new_measurements) or AddMeasurements(const std::list<Meas_>& new_measurements, const TransformationData_& transformation_data) if
 * the tracking frame changes and the measurements and models need to be transformed. All of the measurements added must have the same time stamp and measurements must be given in chronological order. When measurements are given,
 * the tracks are propagated to the current time, measurements and tracks are transformed if transformation data is provided, the new measurements are then associated to tracks, clusters, and the data tree. Finally, the
 * measurements associated to the tracks used to update the tracks. See Meas for information regarding the measurements and the member vairables the user must specify.
 * 
 * Track initialization is the process of generating new tracks from measurements not associated with any current tracks. To run track initialization, call RunTrackInitialization. Track management merges similar tracks, 
 * prunes poor tracks, and ranks the tracks. To run track management, use RunTrackManagement. 
 * 
 * As a user, you will want to know at least which tracks are good. Good tracks are stored in System::good_models_. To get access to all of the system information, use GetSystemInformation.
 * 
 * Please let me know if  you have any questions, or if you want to contribute to making this R-RANSAC version the best it can be. 
 * 
 * 
 * 
 * 
 * @see State
 * @see SourceBase
 * @see ModelBase
 * @see DataAssociationHost
 * @see Ransac
 * @see TransformBase
 * @see MeasurementBase
 * @see Meas
 * 
 */ 


template<typename tRRANSACTemplateParameters>
class RRANSAC {

public:


    typedef typename tRRANSACTemplateParameters::TSourceContainer TSourceContainer;           /**< The container of all of the sources. */
    typedef typename tRRANSACTemplateParameters::TState TState;                               /**< The state of the target. @see State. */
    typedef typename tRRANSACTemplateParameters::TDataType TDataType;                         /**< The scalar object for the data. Ex. float, double, etc. */ 
    typedef typename tRRANSACTemplateParameters::TTransformation TTransformation;             /**< The object type of the measurement and track transformation. @see TransformBase*/
    typedef typename tRRANSACTemplateParameters::TTransformationData TTransformationData;     /**< The object type of the transformation data. @see TransformBase*/
    typedef typename tRRANSACTemplateParameters::TModel TModel;                               /**< The object type of the model. @see ModelBase. */                          
    typedef typename tRRANSACTemplateParameters::TModelManager TModelManager;                 /**< The object type of the model manager. @see ModelManager. */
    typedef typename tRRANSACTemplateParameters::TMeas TMeas;                                 /**< The object type of the measurement. @see Meas */
    typedef typename tRRANSACTemplateParameters::TRansac TRansac;                             /**< The object type of Ransac or the track initializer. @see Ransac */
    typedef typename tRRANSACTemplateParameters::TDataAssociation TDataAssociation;           /**< The object type of the data association method. @see DataAssociationHost */





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
    bool AddSource(const SourceParameters& source_params, std::function<bool(const TState&)> state_in_surveillance_region_callback);

    /**
     * Sets all of the system parameters that are not source parameters. If the process noise covariance changes, then the model manager is used to
     * update the tracks process noise covariance. 
     * \param[in] new_params The new parameters.
     * \return If the parameters were successfully set (i.e. met all of the prescibed requirements), it returns true. Otherwise false.
     */
    bool SetSystemParameters(const Parameters &new_params) { 
        bool update_track_params = false;
        if (system_parameters_set_ && new_params.process_noise_covariance_ != sys_.params_.process_noise_covariance_) {
            update_track_params = true;
        }
        system_parameters_set_ = sys_.params_.SetParameters(new_params);

        if (update_track_params)
            TModelManager::SetModelParameters(sys_);

        return  system_parameters_set_; 
        }

    /**
     * Changes some of the parameters of the source. 
     * \param[in] new_source_params The new source parameters.
     * \param[in] source_index The index to the source whose parameters are to be changed.
     * \return If the parameters were successfully set (i.e. met all of the prescibed requirements), it returns true. Otherwise false.
     */
    bool ChangeSourceParameters(const SourceParameters &new_source_params);


    /**
     * Performs data management by first propagating the tracks to the current time stamp and
     * then associating the measurements to tracks, clusters, or the data tree. Measurements associated to tracks are then used to 
     * update the tracks.
     * @param[in] new_measurements A list that contains new measurements with the same time stamp. 
     * @param[in] current_time There are times when you can call AddMeasurements without giving any measurements. In this case, the time stamp on the 
     *                         measurements cannot be used to get the current time; thus this parameter is needed. 
     * a unique source. 
     */
    void AddMeasurements(const std::list<TMeas>& new_measurements, const double current_time);

    /**
     * \detail Performes data management by first propagating the tracks to the current time stamp, transforming the tracks,
     * and previous measurements to the current tracking frame, and
     * then associating the measurements to tracks, clusters, or the data tree. Measurements associated to tracks are then used to 
     * update the tracks.
     * @param[in] new_measurements A list that contains new measurements with the same time stamp. 
     * @param[in] current_time There are times when you can call AddMeasurements without giving any measurements. In this case, the time stamp on the 
     *                         measurements cannot be used to get the current time; thus this parameter is needed. 
     * @param[in] transformation_data The data required to transform the tracks and measurements.
     * a unique source. 
     */
    void AddMeasurements(const std::list<TMeas>& new_measurements, const double current_time, const TTransformationData& transformation_data);


    /**
     * Runs RANSAC on every cluster in order to initialize new tracks. 
     */ 
    void RunTrackInitialization();

    /**
     * Prunes, merges, and ranks the tracks. 
     */ 
    void RunTrackManagement() {
        if (!system_parameters_set_)
            throw std::runtime_error("System parameters are not set. ");
        
        TModelManager::ManageModels(sys_,sys_.current_time_-sys_.params_.meas_time_window_);
        }
        

    /**
     * Returns a constant pointer to the system which contains all of the R-RANSAC data. 
     */ 
    const System<TModel>* GetSystemInformation() {return &sys_;};

private:

    System<TModel> sys_;                      /**< Contains all of the data for R-RANSAC. */
    TRansac ransac_;                          /**< The track initializer. */
    TDataAssociation data_association_host_;  /**< The data association host. Responsible for associating new measurements to tracks and the data tree. */
    bool transform_data_ = false;             /**< A flag used to indicate if the measurements and tracks should be transformed. */
    bool system_parameters_set_ = false;      /**< A flag used to indicate if the system parameters have been set. true indicates that they have been set. */

    /**
     * When the build type is Debug, checks will be done on the measurements to verify 
     * that they are set correctly. The checks include but not limited to: source index, measurement type, the pose and twist data.
     * If a single measurement doesn't pass a test, then a runtime error will be thrown and non of the measurements will be added.
     * \param[in] new_measurements The list of new measurements that will be verified.
     * \return Returns true if all of the measurements pass the test; otherwise, false. 
     */
    bool VerifyMeasurements(const std::list<TMeas>& new_measurements);

};


//-------------------------------------------------------------------------------------------------------------------------------------------------
//                                                           Definitions
//-------------------------------------------------------------------------------------------------------------------------------------------------

template<typename tRRANSACTemplateParameters>
bool RRANSAC<tRRANSACTemplateParameters>::AddSource(const SourceParameters& source_params) {

    return sys_.source_container_.AddSource(source_params);

}

//-------------------------------------------------------------------------------------------------------------------------------------------------

template<typename tRRANSACTemplateParameters>
bool RRANSAC<tRRANSACTemplateParameters>::AddSource(const SourceParameters& source_params, std::function<bool(const TState&)> state_in_surveillance_region_callback) {

    return sys_.source_container_.AddSource(source_params,state_in_surveillance_region_callback);

}

//-------------------------------------------------------------------------------------------------------------------------------------------------

template<typename tRRANSACTemplateParameters>
bool RRANSAC<tRRANSACTemplateParameters>::ChangeSourceParameters(const SourceParameters &new_params) {
    if (new_params.source_index_ < 0 || new_params.source_index_ >= sys_.source_container_.num_sources_) {
        throw std::runtime_error("RANSAC::ChangeSourceParameters The source index must be greater than 0 and less than " + std::to_string(sys_.source_container_.num_sources_));
        return false;
    }
    else if (sys_.source_container_.GetParams(new_params.source_index_).source_index_ != new_params.source_index_) {
        throw std::runtime_error("RANSAC::ChangeSourceParameters Cannot change the source index.");
        return false;
    } else {
        return sys_.source_container_.ChangeSourceParameters(new_params);
    }
    
}

//-------------------------------------------------------------------------------------------------------------------------------------------------

template<typename tRRANSACTemplateParameters>
bool RRANSAC<tRRANSACTemplateParameters>::VerifyMeasurements(const std::list<TMeas>& new_measurements) {

     double time_stamp = new_measurements.begin()->time_stamp;
     bool success = true;


    if (sys_.time_set_) {
        if (time_stamp < sys_.current_time_) {
            throw std::runtime_error("RRANSAC::VerifyMeasurements The measurement time stamp is less than the current system time stamp. Measurements must be provided in chronological order");
        success = false;
        }
    }

    for (auto meas_iter = new_measurements.begin(); meas_iter != new_measurements.end(); ++meas_iter) {

        if (meas_iter->time_stamp != time_stamp) {
            throw std::runtime_error("RANSAC::VerifyMeasurements All of the measurements must have the same time stamp.");
            success = false;
        } else if (!sys_.source_container_.IsMeasurementAcceptable(*meas_iter)) {
            success = false;
        }
    }

    return success;

}


//-------------------------------------------------------------------------------------------------------------------------------------------------

template<typename tRRANSACTemplateParameters>
void RRANSAC<tRRANSACTemplateParameters>::AddMeasurements(const std::list<TMeas>& new_measurements, const double current_time) {

    if (!system_parameters_set_)
        throw std::runtime_error("System parameters are not set. ");

    if (new_measurements.size() >= 0) {

        sys_.dt_ = current_time - sys_.current_time_;

#ifdef DEBUG_BUILD
        bool correct = VerifyMeasurements(new_measurements);
        if (sys_.dt_ < 0)
            throw std::runtime_error("Measurements must be provided in chronological order. The current time stamp is: " + std::to_string(sys_.current_time_)+ " The measurement time stamp is: " + std::to_string(new_measurements.begin()->time_stamp));

#endif

    
        sys_.current_time_ = current_time; 
        sys_.data_tree_.PruneDataTree(sys_,sys_.current_time_-sys_.params_.meas_time_window_);

        sys_.new_meas_ = new_measurements;

        if (sys_.dt_ > 0) {
            TModelManager::PropagateModels(sys_,sys_.dt_);
        }

// Calculate the innovation covariances used to compute the validation region. This is only for visualization purposes. 
#if RRANSAC_VIZ_HOOKS
        for (auto& track: sys_.models_) {
            track.S_validation_.clear();
            for (auto& source : sys_.sources_)
                track.S_validation_.push_back(track.GetInnovationCovariance(sys_.sources_, source.params_.source_index_));
        }
#endif


        if (transform_data_) {
            sys_.data_tree_.TransformMeasurements(sys_.transformaion_);
            TModelManager::TransformModels(sys_);
            transform_data_ = false;
        }

        data_association_host_.AssociateNewMeasurements(sys_);

        TModelManager::UpdateModels(sys_);

        sys_.data_tree_.ConstructClusters(sys_);


        if (!sys_.time_set_)
            sys_.time_set_ = true;


    } else {
        if (transform_data_) {
            sys_.data_tree_.TransformMeasurements(sys_.transformaion_);
            TModelManager::TransformModels(sys_);
            transform_data_ = false;
        }
        sys_.current_time_ = current_time;
    }

}

//-------------------------------------------------------------------------------------------------------------------------------------------------

template<typename tRRANSACTemplateParameters>
void RRANSAC<tRRANSACTemplateParameters>::AddMeasurements(const std::list<TMeas>& new_measurements, const double current_time, const TTransformationData& transformation_data) {

    sys_.transformaion_.SetData(transformation_data);
    transform_data_ = true;
    AddMeasurements(new_measurements,current_time);

}

//-------------------------------------------------------------------------------------------------------------------------------------------------

template<typename tRRANSACTemplateParameters>
void RRANSAC<tRRANSACTemplateParameters>::RunTrackInitialization() {

    if (!system_parameters_set_)
        throw std::runtime_error("System parameters are not set. ");

    ransac_.Run(sys_);


}

} // namespace rransac

#endif // RRANSAC_RRANSAC_H_