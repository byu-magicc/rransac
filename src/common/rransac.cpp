#include "rransac.h"



namespace rransac
{

bool RRANSAC::AddSource(const SourceParameters& params) {

    // See if source already exists
    for (SourceBase source : sys_.sources_ ) {
        if (params.source_index_ == source.params_.source_index_) {
            throw std::runtime_error("RRANSAC::AddSource Source ID already exists. Cannot add it");
            return false;
        }
    }

    // The new source id must be +1 than the previous id starting at 0
    if (sys_.sources_.size() == 0) {
        if (params.source_index_ != 0) {
            throw std::runtime_error("RRANSAC::AddSource The first source id should be 0. You set it to " + std::to_string(params.source_index_));
            return false;
        }
    } else {
        if ( sys_.sources_.back().params_.source_index_ +1 != params.source_index_) {
            throw std::runtime_error("RRANSAC::AddSource The previous id was " + std::to_string(sys_.sources_[sys_.sources_.size()].params_.source_index_) + ". The new id should be " + std::to_string(sys_.sources_[sys_.sources_.size()].params_.source_index_+1));
            return false;
        }
    }

    if (params.expected_num_false_meas_ <0) {
        throw std::runtime_error("RRANSAC::AddSource Expected number of false measurements cannot be negative.");
        return false;
    }

    if(params.probability_of_detection_ < 0 || params.probability_of_detection_ >1) {
        throw std::runtime_error("RRANSAC::AddSource The probability of detection must be between 0 and 1.");
        return false;
    }

    // If the measurement covariance is fixed, you must supply a measurement
    // covariance
    if (params.meas_cov_fixed_) {
        if (params.meas_cov_.size() == 0) {
            throw std::runtime_error("RRANSAC::AddSource Measurement covariance must be specified when the measurement covariance is fixed.");
            return false;
        }
    }

    // If the measurement type exists, create the new source.
    SourceBase new_source();
    switch (params.type_)
    {
    case MeasurementTypes::R2_POSE:
        new_source.Init<MeasurementTypes::R2_POSE>(params);
        break;
    case MeasurementTypes::R2_POSE_TWIST:
        new_source.Init<MeasurementTypes::R2_POSE_TWIST>(params);
        break;
    default:
        throw std::runtime_error("RRANSAC::AddSource Source type does not exist");
        return false;
        break;
    }

    // Source does exist so add it
    sys_.sources_.emplace_back(new_source);
    return true;
    }

    
} // namespace rransac


