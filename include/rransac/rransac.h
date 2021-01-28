#ifndef RRANSAC_RRANSAC_H_
#define RRANSAC_RRANSAC_H_

#include "system.h"
#include "track_initialization/ransac.h"
#include "common/data_association/data_association_host.h"
#include "common/models/model_manager.h"
#include <functional>
#include <string>


namespace rransac {

/**
 * \class RRANSACTemplateParameters
 * Contains all of the template parameters for RRANSAC and used to initialize a RRANSAC object
 */ 
template<typename tModel, template<class> typename tModelDataAssociationPolicyClass, template<class> typename tClusterDataTreeAssociationPolicyClass, template <typename > typename tSeed, template<typename , template <typename > typename > typename tLMLEPolicy>
struct RRANSACTemplateParameters {

    typedef DataAssociationHost<tModel,tModelDataAssociationPolicyClass,tClusterDataTreeAssociationPolicyClass> DataAssociation_;
    typedef typename tModel::DataType DataType_;
    typedef typename tModel::Source Source_;
    typedef typename tModel::State State_;
    typedef typename tModel::Transformation Transformation_;
    typedef typename tModel::Transformation::Data TransformationData_;
    typedef ModelManager<tModel> ModelManager_;
    typedef Ransac<tModel,tSeed,tLMLEPolicy,tModelDataAssociationPolicyClass> Ransac_;
    typedef Meas<DataType_> Meas_;

};

/**
 * \class RRANSAC
 * This class is the interface for the user. 
 * 
 */ 

// Model
// -- State
// -- Source
// -- Transformation
// -- DataType
// tModelDataAssociationPolicyClass
// tClusterDataTreeAssociationPolicyClass
// ransac
// -- tSeed
// -- tLMLEPolicy

template<typename tRRANSACTemplateParameters>
class RRANSAC {

public:

    typedef typename tRRANSACTemplateParameters::DataAssociation_ DataAssociation_;
    typedef typename tRRANSACTemplateParameters::DataType_ DataType_;
    typedef typename tRRANSACTemplateParameters::Source_ Source_;
    typedef typename tRRANSACTemplateParameters::State_ State_;
    typedef typename tRRANSACTemplateParameters::Transformation_ Transformation_;
    typedef typename tRRANSACTemplateParameters::TransformationData_ TransformationData_;
    typedef typename tRRANSACTemplateParameters::ModelManager_ ModelManager_;
    typedef typename tRRANSACTemplateParameters::Ransac_ Ransac_;
    typedef typename tRRANSACTemplateParameters::Meas_ Meas_;

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
    bool AddSource(const SourceParameters& params, std::function<bool(const State_&)> state_in_surveillance_region_callback);

    /**
     * \detail Sets all of the system parameters
     * \param[in] new_params The new parameters.
     * \return If the parameters were successfully set (i.e. met all of the prescibed requirements), it returns true. Otherwise false.
     */
    bool SetSystemParameters(const Parameters &new_params) { return sys_.params_.SetParameters(new_params);  }

    /**
     * \detail Changes some of the parameters of the source. 
     * \param[in] new_params The new source parameters.
     * \param[in] source_index The index to the source whose parameters are to be changed.
     * \return If the parameters were successfully set (i.e. met all of the prescibed requirements), it returns true. Otherwise false.
     */
    bool ChangeSourceParameters(const SourceParameters &new_params);


    /**
     * \detail Performs data management by first propagating the tracks to the current time stamp and
     * then associating the measurements to tracks, clusters, or the data tree. Measurements associated to tracks are then used to 
     * update the tracks.
     * @param new_measurements A list that contains new measurements with the same time stamp. 
     * a unique source. 
     */
    void AddMeasurements(const std::list<Meas_>& new_measurements);

    /**
     * \detail Performes data management by first propagating the tracks to the current time stamp, transforming the tracks,
     * and previous measurements to the current tracking frame, and
     * then associating the measurements to tracks, clusters, or the data tree. Measurements associated to tracks are then used to 
     * update the tracks.
     * @param new_measurements A list that contains new measurements with the same time stamp. 
     * @param transformation_data The data required to transform the tracks and measurements.
     * a unique source. 
     */
    void AddMeasurements(const std::list<Meas_>& new_measurements, const TransformationData_& transformation_data);


    /**
     * \detail Runs RANSAC on every cluster in order to initialize new tracks. 
     */ 
    void RunTrackInitialization();

    /**
     * \detail Prunes, merges, and ranks the tracks. 
     */ 
    void RunTrackManagement() {ModelManager_::ManageModels(sys_,sys_.current_time_-sys_.params_.meas_time_window_);}

    /**
     * Returns a constant pointer to the system which contains all of the R-RANSAC data. 
     */ 
    const System<tModel>* GetSystemInformation() {return &sys_};

private:

    System<tModel> sys_;
    Ransac_ ransac_;

        /**
     * \detail When the build type is Debug, checks will be done on the measurements to verify 
     * that they are set correctly. The checks include but not limited to: source index, measurement type, the pose and twist data.
     * If a single measurement doesn't pass a test, then a runtime error will be thrown and non of the measurements will be added.
     * \param[in] new_measurements The list of new measurements that will be verified.
     * \return Returns true if all of the measurements pass the test; otherwise, false. 
     */
    bool VerifyMeasurements(const std::list<Meas_>& new_measurements);

};


//-------------------------------------------------------------------------------------------------------------------------------------------------
//                                                           Definitions
//-------------------------------------------------------------------------------------------------------------------------------------------------

template<typename tRRANSACTemplateParameters>
bool RRANSAC<tRRANSACTemplateParameters>::AddSource(const SourceParameters& params) {

    if (params.source_index_ != sys_.source_index_counter_) {
        throw std::runtime_error("RRANSAC::AddSource Sources must be added in consecutive sequential order of their source index starting from 0. 
        Your source index is " + std::to_string(params.source_index_) + " when it should be " + std::to_string(sys_.source_index_counter_));
        return false;
    } else {
        Source_ source;
        source.init(params);
        sys_.sources_.push_back(source);
        ++sys_.source_index_counter_;
        return true;
        
    }

}

//-------------------------------------------------------------------------------------------------------------------------------------------------

template<typename tRRANSACTemplateParameters>
bool RRANSAC<tRRANSACTemplateParameters>::AddSource(const SourceParameters& params, std::function<bool(const State_&)> state_in_surveillance_region_callback) {

    if (params.source_index_ != sys_.source_index_counter_) {
        throw std::runtime_error("RRANSAC::AddSource Sources must be added in consecutive sequential order of their source index starting from 0. 
        Your source index is " + std::to_string(params.source_index_) + " when it should be " + std::to_string(sys_.source_index_counter_));
        return false;
    } else {
        Source_ source;
        source.init(params,state_in_surveillance_region_callback);
        sys_.sources_.push_back(source);
        ++sys_.source_index_counter_;
        return true;
        
    }

}

//-------------------------------------------------------------------------------------------------------------------------------------------------

template<typename tRRANSACTemplateParameters>
bool RRANSAC<tRRANSACTemplateParameters>::ChangeSourceParameters(const SourceParameters &new_params) {
    if (new_params.source_index < 0 || source_index >= sys_.sources_.size())
        throw std::runtime_error("RANSAC::ChangeSourceParameters The source index must be greater than 0 and less than " + std::to_string(sys_.sources_.size()));
        return false;
    else if (sys_.sources_[source_index].params_.source_index_ != new_params.source_index) {
        throw std::runtime_error("RANSAC::ChangeSourceParameters Cannot change the source index.")
        return false;
    } else {
        return sys_.sources_[source_index].ChangeSourceParameters(new_params);
    }
    
}

//-------------------------------------------------------------------------------------------------------------------------------------------------

template<typename tRRANSACTemplateParameters>
bool RRANSAC<tRRANSACTemplateParameters>::VerifyMeasurements(const std::list<Meas_>& new_measurements) {

     double time_stamp = new_measurements.begin()->time_stamp;
     bool success = true;


    if (sys_.time_set_) {
        if (time_stamp < sys_.current_time_) {
            throw std::runtime_error("RRANSAC::VerifyMeasurements The measurement time stamp is less than the current system time stamp. Measurements must be provided in chronological order")
        success = false;
        }
    }

    for (auto meas_iter = new_measurements.begin(); meas_iter != new_measurements.end(); ++meas_iter) {

        if (meas_iter->time_stamp != time_stamp) {
            throw std::runtime_error("RANSAC::VerifyMeasurements All of the measurements must have the same time stamp.");
            success = false;
        }

        if (meas_iter->source_index < 0 || meas_iter->source_index > sys_.sources_.size()) {
            throw std::runtime_error("RANSAC::VerifyMeasurements Source index of measurement is out of scope.");
            success = false;
        }

        if (meas_iter->type != sys_.sources_[meas_iter->source_index].type_) {
            throw std::runtime_error("RANSAC::VerifyMeasurements Measurement type does not match the source's measurement type. Make sure the source index is correct.");
            success = false;            
        }

        if (meas_iter->pose.rows() != Source_::meas_dim_) {
            throw std::runtime_error("RANSAC::VerifyMeasurements The pose of the measurement is not the correct dimension.");
            success = false;
        }

        if (sys_.source_[meas_iter->source_index].params_.has_twist && meas_iter->twist.rows() != Source_::meas_dim_) {
            throw std::runtime_error("RANSAC::VerifyMeasurements The twist of the measurement is not the correct dimension.");
            success = false;
        }

    }


}


//-------------------------------------------------------------------------------------------------------------------------------------------------

template<typename tRRANSACTemplateParameters>
void RRANSAC<tRRANSACTemplateParameters>::AddMeasurements(const std::list<Meas_>& new_measurements) {

    if (new_measurements.size() > 0) {

#ifdef DEBUG_BUILD
    bool correct = VerifyMeasurements(new_measurements);
#endif

    double dt = new_measurements.begin()->time_stamp - sys_.current_time_;
    sys_.current_time_ = new_measurements.begin()->time_stamp; 
    sys_.new_meas_ = new_measurements;
    if (dt > 0)
        ModelManager_::PropagateModels(sys_,dt);
    sys_.data_tree_.PruneDataTree(sys_,sys_.current_time_-sys_.params_.meas_time_window_);
    DataAssociation_::AssociateNewMeasurements(sys_);
    ModelManager_::UpdateModels(sys_);
    sys_.data_tree_.ConstructClusters(sys_);

    if (!sys_.time_set_)
        sys_.time_set_ = true;


    }

}

//-------------------------------------------------------------------------------------------------------------------------------------------------

template<typename tRRANSACTemplateParameters>
void RRANSAC<tRRANSACTemplateParameters>::AddMeasurements(const std::list<Meas_>& new_measurements, const TransformationData_& transformation_data) {


    sys_.transformaion_.SetData(transformation_data);
    sys_.data_tree_.TransformMeasurements(sys_.transformaion_);
    ModelManager_::TransformModels(sys_);
    AddMeasurements(new_measurements);
}

//-------------------------------------------------------------------------------------------------------------------------------------------------

template<typename tRRANSACTemplateParameters>
void RRANSAC<tRRANSACTemplateParameters>::RunTrackInitialization() {

    ransac_.Run(sys_);

}

} // namespace rransac

#endif // RRANSAC_RRANSAC_H_