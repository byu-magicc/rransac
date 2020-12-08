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
    RANSAC_stopping_criteria_ = 0.5;
}

//------------------------------------------------------------------------------------
Parameters::Parameters(const Parameters &new_params) {
    SetParameters(new_params);
}

//------------------------------------------------------------------------------------
bool Parameters::SetParameters(const Parameters &new_params) {

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
    if (new_params.RANSAC_stopping_criteria_ <=0 || new_params.RANSAC_stopping_criteria_ >=1)
    {
        throw std::runtime_error("Parameters::SetParameters The provided value of RANSAC_stopping_criteria_ is not between 0 and 1.");
        successfull = false;
    }
    else
    {
        RANSAC_stopping_criteria_ = new_params.RANSAC_stopping_criteria_;
    }


    return successfull;

}

//------------------------------------------------------------------------------------
Parameters::~Parameters() = default;





}
