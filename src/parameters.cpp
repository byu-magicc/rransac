#include "rransac/parameters.h"

#include <iostream>
#include <stdexcept>

namespace rransac {

//------------------------------------------------------------------------------------
Parameters::Parameters()
{
    meas_time_window_ = 5;
    transform_consensus_set_ = false;
    
    cluster_time_threshold_ = 0.5;
    cluster_velocity_threshold_ = 2;
    cluster_position_threshold_ = 1;
    cluster_min_size_requirement_ = 10;

    RANSAC_max_iters_ = 50;
    RANSAC_score_stopping_criteria_ = 10;
    RANSAC_score_minimum_requirement_ = 5;
    RANSAC_minimum_subset_ = 3;
    
    good_model_threshold_ = 100;
    max_missed_detection_time_ = 5;
    similar_tracks_threshold_ = 1;
    max_num_models_ = 10;
    
    NonLinearInnovCovId_ = false;
    NonLinearLMLECeresThreads_= 1;
    NonLinearLMLECeresMaxNumIters_ = 50;
}

//------------------------------------------------------------------------------------
Parameters::Parameters(const Parameters &new_params) {
    SetParameters(new_params);
}

//------------------------------------------------------------------------------------

void Parameters::operator= (const Parameters &new_params) {
    SetParameters(new_params);
}

//------------------------------------------------------------------------------------
bool Parameters::SetParameters(const Parameters &new_params) {

    bool successfull = true;
    transform_consensus_set_ = new_params.transform_consensus_set_;
    NonLinearInnovCovId_ = new_params.NonLinearInnovCovId_;

    if (new_params.max_num_models_ <= 0 ) {
        throw std::runtime_error("Parameters::SetParameters The provided value of max_num_models_ has not been initialized");
        successfull = false;
    } else {
        max_num_models_ = new_params.max_num_models_;
    }


    if (new_params.similar_tracks_threshold_ <= 0 ) {
        throw std::runtime_error("Parameters::SetParameters The provided value of similar_tracks_threshold_ has not been initialized");
        successfull = false;
    } else {
        similar_tracks_threshold_ = new_params.similar_tracks_threshold_;
    }

    if (new_params.max_missed_detection_time_ <= 0 ) {
        throw std::runtime_error("Parameters::SetParameters The provided value of max_missed_detection_time_ has not been initialized");
        successfull = false;
    } else {
        max_missed_detection_time_ = new_params.max_missed_detection_time_;
    }

    if (new_params.good_model_threshold_ <= 0 ) {
        throw std::runtime_error("Parameters::SetParameters The provided value of good_model_threshold_ has not been initialized");
        successfull = false;
    } else {
        good_model_threshold_ = new_params.good_model_threshold_;
    }
    
    if (new_params.process_noise_covariance_.rows() <= 0 ) {
        throw std::runtime_error("Parameters::SetParameters The provided value of process_noise_covariance_ has not been initialized");
        successfull = false;
    } else {
        process_noise_covariance_ = new_params.process_noise_covariance_;
    }

    if (new_params.NonLinearLMLECeresMaxNumIters_ <= 0 ) {
        throw std::runtime_error("Parameters::SetParameters The provided value of NonLinearLMLECeresMaxNumIters_ is less than or equal to 0. It must be greater than 0.");
        successfull = false;
    } else {
        NonLinearLMLECeresMaxNumIters_ = new_params.NonLinearLMLECeresMaxNumIters_;
    }

    if (new_params.NonLinearLMLECeresThreads_ <= 0 ) {
        throw std::runtime_error("Parameters::SetParameters The provided value of NonLinearLMLECeresThreads_ is less than or equal to 0. It must be greater than 0.");
        successfull = false;
    } else {
        NonLinearLMLECeresThreads_ = new_params.NonLinearLMLECeresThreads_;
    }

    if (new_params.cluster_min_size_requirement_ <= 0 ) {
        throw std::runtime_error("Parameters::SetParameters The provided value of cluster_min_size_requirement_ is less than or equal to 0. It must be greater than 0.");
        successfull = false;
    } else {
        cluster_min_size_requirement_ = new_params.cluster_min_size_requirement_;
    }

    if (new_params.RANSAC_score_minimum_requirement_ <= 0 ) {
        throw std::runtime_error("Parameters::SetParameters The provided value of RANSAC_score_minimum_requirement_ is less than or equal to 0. It must be greater than 0.");
        successfull = false;
    } else {
        RANSAC_score_minimum_requirement_ = new_params.RANSAC_score_minimum_requirement_;
    }

    if (new_params.RANSAC_minimum_subset_ <= 0 ) {
        throw std::runtime_error("Parameters::SetParameters The provided value of RANSAC_minimum_subset_ is less than or equal to 0. It must be greater than 0.");
        successfull = false;
    } else {
        RANSAC_minimum_subset_ = new_params.RANSAC_minimum_subset_;
    }

    if (new_params.cluster_velocity_threshold_ <= 0 ) {
        throw std::runtime_error("Parameters::SetParameters The provided value of cluster_velocity_threshold_ is less than or equal to 0. It must be greater than 0.");
        successfull = false;
    } else {
        cluster_velocity_threshold_ = new_params.cluster_velocity_threshold_;
    }


    if (new_params.cluster_time_threshold_ <= 0 ) {
        throw std::runtime_error("Parameters::SetParameters The provided value of cluster_time_threshold_ is less than or equal to 0. It must be greater than 0.");
        successfull = false;
    } else {
        cluster_time_threshold_ = new_params.cluster_time_threshold_;
    }


    if (new_params.cluster_position_threshold_ <= 0 ) {
        throw std::runtime_error("Parameters::SetParameters The provided value of cluster_position_threshold_ is less than or equal to 0. It must be greater than 0.");
        successfull = false;
    } else {
        cluster_position_threshold_ = new_params.cluster_position_threshold_;
    }
    

    //Ensure that the new measurement time window is greater than 0.
    if (new_params.meas_time_window_ <= 0)
    {
        throw std::runtime_error("Parameters::SetParameters The provided value of meas_time_window_ is less than or equal to 0. It must be greater than 0.");
        successfull = false;
    }
    else
    {
        meas_time_window_ = new_params.meas_time_window_;
    }

    // Ensure that max_RANSAC_iters_ has a realistic value.
    if (new_params.RANSAC_max_iters_ <=0)
    {
        throw std::runtime_error("Parameters::SetParameters The provided value of max_RANSAC_iters_ is not greater than 0.");
        successfull = false;
    }
    else
    {
        RANSAC_max_iters_ = new_params.RANSAC_max_iters_;
    }

    // Ensure that RANSAC_stopping_criteria_ has a realistic value.
    if (new_params.RANSAC_score_stopping_criteria_ <=0)
    {
        throw std::runtime_error("Parameters::SetParameters The provided value of RANSAC_score_stopping_criteria_ must be greater than 1.");
        successfull = false;
    }
    else
    {
        RANSAC_score_stopping_criteria_ = new_params.RANSAC_score_stopping_criteria_;
    }


    return successfull;

}

//------------------------------------------------------------------------------------
Parameters::~Parameters() = default;





}
