#include "rransac/parameters.h"

#include <iostream>
#include <stdexcept>

namespace rransac {

//------------------------------------------------------------------------------------
Parameters::Parameters()
{
    fixed_time_interval_ = false;
    time_interval_ = 1;
    transform_consensus_set_ = false;
    meas_time_window_ = 1;
    RANSAC_max_iters_ = 5000;
    RANSAC_score_stopping_criteria_ = 10;
    RANSAC_score_minimum_requirement_ = 5;
    cluster_min_size_requirement_ = 10;
    NonLinearLMLECeresThreads_= 1;
    NonLinearInnovCovId_ = false;
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

    cluster_position_threshold_ = new_params.cluster_position_threshold_;
    cluster_time_threshold_ = new_params.cluster_time_threshold_;
    cluster_velocity_threshold_ = new_params.cluster_velocity_threshold_;
    RANSAC_minimum_subset_ = new_params.RANSAC_minimum_subset_;
    process_noise_covariance_ = new_params.process_noise_covariance_;
    RANSAC_score_minimum_requirement_ = new_params.RANSAC_score_minimum_requirement_;
    cluster_min_size_requirement_ = new_params.cluster_min_size_requirement_;
    NonLinearLMLECeresThreads_ = new_params.NonLinearLMLECeresThreads_;
    NonLinearInnovCovId_ = new_params.NonLinearInnovCovId_;
    NonLinearLMLECeresMaxNumIters_ = new_params.NonLinearLMLECeresMaxNumIters_;

    bool successfull = true;

    fixed_time_interval_ = new_params.fixed_time_interval_;

    // Ensure that the new time interval is greater than 0.
    if (new_params.time_interval_ <= 0)
    {
        throw std::runtime_error("Parameters::SetParameters The provided value of time_interval_ is less than or equal to 0. It must be greater than 0.");
        
        successfull = false;
    }
    else
    {
        time_interval_ = new_params.time_interval_;
    }

    transform_consensus_set_ = new_params.transform_consensus_set_;

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
