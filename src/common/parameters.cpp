#include "rransac/common/parameters.h"

#include <iostream>

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


//------------------------------------------------------------------------------------
bool Parameters::AddSource(const SourceParameters& params) {


// See if source already exists
for (SourceBase source : sources_ ) {
    if (params.source_id_ == source.params_.source_id_) {
        throw std::runtime_error("Parameters::AddSource Source ID already exists. Cannot add it");
        return false;
    }
}

// The new source id must be +1 than the previous id starting at 0
if (sources_.size() == 0) {
    if (params.source_id_ != 0) {
        throw std::runtime_error("The first source id should be 0. You set it to " + std::to_string(params.source_id_));
        return false;
    }
} else {
    if ( sources_.back().params_.source_id_ +1 != params.source_id_) {
        throw std::runtime_error("The previous id was " + std::to_string(sources_[sources_.size()].params_.source_id_) + ". The new id should be " + std::to_string(sources_[sources_.size()].params_.source_id_+1));
        return false;
    }
}

if (params.expected_num_false_meas_ <0) {
    throw std::runtime_error("Parameters::AddSource Expected number of false measurements cannot be negative.");
    return false;
}

// If the measurement covariance is fixed, you must supply a measurement
// covariance
if (params.meas_cov_fixed_) {
    if (params.meas_cov_.size() == 0) {
        throw std::runtime_error("Parameters::AddSource Measurement covariance must be specified when the measurement covariance is fixed.");
        return false;
    }
}

// If the source type exists, create the new source.
SourceBase new_source;
switch (params.type_)
{
case SourceTypes::R2_POS:
    new_source.Init<SourceTypes::R2_POS>(params);
    break;
case SourceTypes::R2_POS_VEL:
    new_source.Init<SourceTypes::R2_POS_VEL>(params);
    break;
default:
    throw std::runtime_error("Parameters::AddSource Source type does not exist");
    return false;
    break;
}

// Source does not exist so add it
sources_.emplace_back(new_source);
return true;
}



}
