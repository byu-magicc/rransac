#include <iostream>

#include "common/param.h"


namespace rransac {

//------------------------------------------------------------------------------------
Parameters::Parameters()
{
    fixed_time_interval_ = false;
    time_interval_ = 1;
    transform_consensus_set_ = false;
    meas_time_window_ = 1;
    probability_of_detection_ = 0.8;
    expected_num_false_meas_ = 10;
    max_RANSAC_iters_ = 5000;
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
        std::cerr << "Parameters: The provided value of time_interval_ is less than or equal to 0. It must be greater than 0. Keeping the current value of " << time_interval_ << std::endl;
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
        std::cerr << "Parameters: The provided value of meas_time_window_ is less than or equal to 0. It must be greater than 0. Keeping the current value of " << meas_time_window_ << std::endl;
        successfull = false;
    }
    else
    {
        meas_time_window_ = new_params.meas_time_window_;
    }

    // Ensure that the probability of detection is between 0 and 1
    if (new_params.probability_of_detection_ <=0 || new_params.probability_of_detection_ >=1)
    {
        std::cerr << "Parameters: The provided value of probability_of_detection_ is not between 0 and 1. Keeping the current value of " << probability_of_detection_ << std::endl;
        successfull = false;
    }
    else
    {
        probability_of_detection_ = new_params.probability_of_detection_;
    }

    // Ensure that expected_num_false_meas_ has a realistic value.
    if (new_params.expected_num_false_meas_ <0)
    {
        std::cerr << "Parameters: The provided value of expected_num_false_meas_ is not greater than or equal to 0. Keeping the current value of " << expected_num_false_meas_ << std::endl;
        successfull = false;
    }
    else
    {
        expected_num_false_meas_ = new_params.expected_num_false_meas_;
    }

    // Ensure that max_RANSAC_iters_ has a realistic value.
    if (new_params.max_RANSAC_iters_ <=0)
    {
        std::cerr << "Parameters: The provided value of max_RANSAC_iters_ is not greater than 0. Keeping the current value of " << max_RANSAC_iters_ << std::endl;
        successfull = false;
    }
    else
    {
        max_RANSAC_iters_ = new_params.max_RANSAC_iters_;
    }
    
    // Ensure that RANSAC_stopping_criteria_ has a realistic value.
    if (new_params.RANSAC_stopping_criteria_ <=0 || new_params.RANSAC_stopping_criteria_ >=1)
    {
        std::cerr << "Parameters: The provided value of RANSAC_stopping_criteria_ is not between 0 and 1. Keeping the current value of " << RANSAC_stopping_criteria_ << std::endl;
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