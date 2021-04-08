#include "rransac/parameters.h"


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
    
    track_good_model_threshold_ = 0.8;
    track_max_missed_detection_time_ = 5;
    track_similar_tracks_threshold_ = 1;
    track_max_num_tracks_ = 10;
    
    nonlinear_innov_cov_id_ = false;
    nonlinear_LMLE_Ceres_threads_= 1;
    nonlinear_LMLE_Ceres_max_num_iters_ = 50;

    set_initial_error_covariance_to_id_ = true;
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
    nonlinear_innov_cov_id_ = new_params.nonlinear_innov_cov_id_;
    set_initial_error_covariance_to_id_ = new_params.set_initial_error_covariance_to_id_;

    if (new_params.track_max_num_tracks_ <= 0 ) {
        throw std::runtime_error("Parameters::SetParameters The provided value of track_max_num_tracks_ has not been initialized");
        successfull = false;
    } else {
        track_max_num_tracks_ = new_params.track_max_num_tracks_;
    }


    if (new_params.track_similar_tracks_threshold_ <= 0 ) {
        throw std::runtime_error("Parameters::SetParameters The provided value of track_similar_tracks_threshold_ has not been initialized");
        successfull = false;
    } else {
        track_similar_tracks_threshold_ = new_params.track_similar_tracks_threshold_;
    }

    if (new_params.track_max_missed_detection_time_ <= 0 ) {
        throw std::runtime_error("Parameters::SetParameters The provided value of track_max_missed_detection_time_ has not been initialized");
        successfull = false;
    } else {
        track_max_missed_detection_time_ = new_params.track_max_missed_detection_time_;
    }

    if (new_params.track_good_model_threshold_ <= 0 || new_params.track_good_model_threshold_ > 1) {
        throw std::runtime_error("Parameters::SetParameters The provided value of track_good_model_threshold_ must be greater than 0 and less than or equal to 1");
        successfull = false;
    } else {
        track_good_model_threshold_ = new_params.track_good_model_threshold_;
    }
    
    if (new_params.process_noise_covariance_.rows() <= 0 ) {
        throw std::runtime_error("Parameters::SetParameters The provided value of process_noise_covariance_ has not been initialized");
        successfull = false;
    } else {
        process_noise_covariance_ = new_params.process_noise_covariance_;
    }

    
    if (!new_params.set_initial_error_covariance_to_id_) {
        if (new_params.initial_error_covariance_.rows() <= 0 ) {
            throw std::runtime_error("Parameters::SetParameters The provided value of initial_error_covariance_ has not been initialized");
            successfull = false;
        } else if(new_params.initial_error_covariance_.rows() != new_params.initial_error_covariance_.cols()) {
            throw std::runtime_error("Parameters::SetParameters The provided value of initial_error_covariance_ is not a square matrix");
            successfull = false;
        } else if( (new_params.initial_error_covariance_ - new_params.initial_error_covariance_.transpose() ).norm() > 1e-9   ) {
            throw std::runtime_error("Parameters::SetParameters The provided value of initial_error_covariance_ is not symmetric");
            successfull = false;
        } else {
            
            Eigen::VectorXcd eigenvalues = new_params.initial_error_covariance_.eigenvalues();
            for (int ii=0; ii < eigenvalues.rows(); ++ii) {
                if ( std::real(eigenvalues(ii)) <= 0) {
                    throw std::runtime_error("Parameters::SetParameters The provided value of initial_error_covariance_ is not positive definite");
                    successfull = false;
                }
            }

            if (successfull)
                initial_error_covariance_ = new_params.initial_error_covariance_;
        }
    }

    if (new_params.nonlinear_LMLE_Ceres_max_num_iters_ <= 0 ) {
        throw std::runtime_error("Parameters::SetParameters The provided value of nonlinear_LMLE_Ceres_max_num_iters_ is less than or equal to 0. It must be greater than 0.");
        successfull = false;
    } else {
        nonlinear_LMLE_Ceres_max_num_iters_ = new_params.nonlinear_LMLE_Ceres_max_num_iters_;
    }

    if (new_params.nonlinear_LMLE_Ceres_threads_ <= 0 ) {
        throw std::runtime_error("Parameters::SetParameters The provided value of nonlinear_LMLE_Ceres_threads_ is less than or equal to 0. It must be greater than 0.");
        successfull = false;
    } else {
        nonlinear_LMLE_Ceres_threads_ = new_params.nonlinear_LMLE_Ceres_threads_;
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
