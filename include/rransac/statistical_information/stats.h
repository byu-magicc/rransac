#ifndef RRANSAC_STATISTICAL_INFORMATION_H_
#define RRANSAC_STATISTICAL_INFORMATION_H_

#include <Eigen/Core>
#include "rransac/system.h"
#include "rransac/common/measurement/measurement_base.h"
#include <numeric>
#include <algorithm>
#include <iostream>
#include <fstream>

namespace rransac {


/**
 * \class Stats
 * This class generates statistical data to serve as performance measures of R-RANSAC.
 * The performance measures are defined in the paper Gorji, A. A., Tharmarasa, R., & Kirubarajan, T. (2011). 
 * Performance measures for multiple target tracking problems. Fusion 2011 - 14th International Conference on Information Fusion, 1560–1567.
 * 
 * Some of the performance measures requires truth data to be provided. 
 * 
 */ 
template<typename tModel>
class Stats {

public:

// Track cardinality measures
std::vector<double> nvt_;             /**< The number of valid tracks. A track is valid if it is assigned to only one target and also, the assigned target is not associated with any other track.*/
std::vector<double> nft_;             /**< The number of false tracks. A track is detected as a false one if it isn't assigned to any target. */
std::vector<double> nmt_;             /**< The number of missed targets. A target is missed if it is not associated with any track. */
std::vector<double> nst_;             /**< The number of spurious tracks. A track is spurious is it is assigned to more than one target. */
std::vector<double> moc_;             /**< The measure of completeness is the number of valid tracks divided by the number of targets. */
std::vector<double> cm_time_;         /**< The time at which the cardinality measures were calculated. */

// General data
std::vector<double> num_good_tracks_;  /**< The number of good tracks currently existing. */
std::vector<double> num_targets_;      /**< The number of targets. */
std::vector<double> num_tracks_;       /**< The number of tracks currently existing. */
std::vector<double> total_num_tracks_; /**< The accumulative number of tracks that have existed since R-RANSAC began. */

// Track Accuracy
std::vector<double> average_euclidean_err_; /**< On each iteration, if a target is associated with a track, it's error will be computed. The sum
                                                 of the norm of the error divided by the number of tracks. */


// Accumulative measures.
std::vector<double> acc_nvt_;                    /**< The number of valid tracks. A track is valid if it is assigned to only one target and also, the assigned target is not associated with any other track.*/
std::vector<double> acc_nft_;                    /**< The number of false tracks. A track is detected as a false one if it isn't assigned to any target. */
std::vector<double> acc_nmt_;                    /**< The number of missed targets. A target is missed if it is not associated with any track. */
std::vector<double> acc_nst_;                    /**< The number of spurious tracks. A track is spurious is it is assigned to more than one target. */
std::vector<double> acc_moc_;                    /**< The measure of completeness is the number of valid tracks divided by the number of targets. */
std::vector<double> acc_cm_time_;                /**< The time at which the cardinality measures were calculated. */
std::vector<double> acc_num_good_tracks_;        /**< The number of good tracks currently existing. */
std::vector<double> acc_num_targets_;            /**< The number of targets. */
std::vector<double> acc_num_tracks_;             /**< The number of tracks currently existing. */
std::vector<double> acc_total_num_tracks_;       /**< The accumulative number of tracks that have existed since R-RANSAC began. */
std::vector<double> acc_average_euclidean_err_;  /**< On each iteration, if a target is associated with a track, it's error will be computed. The sum
                                                  of the norm of the error divided by the number of tracks. */



double tpd_;  /**< Track probability of detection for one monte carlo simulation. The percentage of times a target was detected. */
double taee_; /**< The total average euclidean error. */
double tant_; /**< The total average number of tracks. */

// Other variables
bool first_run_ = true;
double num_monte_carl_runs_ = 0;


/**
 * Clears all of the performance measures;
 */ 
void Reset();

/**
 * Gets ready for the next monte carlo simulation by adding the performance measures to 
 * the accumulative measures and clearing the performance measures. 
 */ 
void Next();

/**
 * Calculates all of the performance measures at each time step of the monte carlo simulation. 
 * @param sys          A pointer to the R-RANSAC system.
 * @param targets      The true tracks.
 * @param source_index The index to the source used to measure the error between the targets and tracks. 
 */ 
void CalculatePerformanceMeasures(const System<tModel>* sys, const std::vector<tModel>& targets, int source_index );

/**
 * Averages the accumulatve measures by the number of monte carlo simulations, and 
 * computes some other performance measures.
 * 
 */ 
void CalculateTotalPerformanceMeasures();

/**
 * Writes the accumulative performance measures and tpd_, taee_, and tant_ to a file that can be read in matlab.
 * The file name is absolute_file_path/output.csv
 * @param absolute_file_path The absolute file path where the data will be saved. 
 */ 
void WriteData(const std::string absolute_file_path);

private:

/**
 * The association matrix assigns a target to either a track or no track. If a target
 * is not in the validation region of a good track according to the source index provided, 
 * it is assigned to no track. If a target is in the validation region of multiple tracks,
 * it is assigned to the closest track according to the metric of the validation region
 * 
 * Let N denote the number of targets, and M the number of good tracks. The
 * association matrix is a Nx(M+1) matrix with binary entries {0,1}. A 1 indicates 
 * the association between a target and a track. The last column is used to indicate that
 * the target was not associated with any track. 
 * @param sys          A pointer to the R-RANSAC system.
 * @param targets      The true tracks.
 * @param source_index The index to the source used to measure the error between the targets and tracks. 
 */ 
Eigen::MatrixXd ComputeAssociationMatrix(const System<tModel>* sys, const std::vector<tModel>& targets, int source_index );

/**
 * Adds the elements of v2 to v1 element wise.
 */ 
template<typename T>
void AddVectorElementWise(std::vector<T>& v1, const std::vector<T>& v2);

}; 

//--------------------------------------------------------------------------------------------------------

template<typename tModel>
void Stats<tModel>::Reset() {
    nvt_.clear();
    nft_.clear();
    nmt_.clear();
    nst_.clear();
    cm_time_.clear();
    moc_.clear();
    num_good_tracks_.clear();
    num_targets_.clear();
    num_tracks_.clear();
    total_num_tracks_.clear();
    average_euclidean_err_.clear();

    acc_nvt_.clear();
    acc_nft_.clear();
    acc_nmt_.clear();
    acc_nst_.clear();
    acc_cm_time_.clear();
    acc_moc_.clear();
    acc_num_good_tracks_.clear();
    acc_num_targets_.clear();
    acc_num_tracks_.clear();
    acc_total_num_tracks_.clear();
    acc_average_euclidean_err_.clear();

    first_run_ = true;
    num_monte_carl_runs_ = 0;


}

//--------------------------------------------------------------------------------------------------------

template<typename tModel>
void Stats<tModel>::Next() {


    if(first_run_) {

        acc_nvt_= nvt_;
        acc_nft_=nft_;
        acc_nmt_=nmt_;
        acc_nst_=nst_;
        acc_cm_time_=cm_time_;
        acc_moc_=moc_;
        acc_num_good_tracks_=num_good_tracks_;
        acc_num_targets_=num_targets_;
        acc_num_tracks_=num_tracks_;
        acc_total_num_tracks_ = total_num_tracks_;
        acc_average_euclidean_err_=average_euclidean_err_;
        first_run_ = false;
        


    } else {

        AddVectorElementWise(acc_nvt_, nvt_);
        AddVectorElementWise(acc_nft_, nft_);
        AddVectorElementWise(acc_nmt_, nmt_);
        AddVectorElementWise(acc_nst_, nst_);
        // AddVectorElementWise(acc_cm_time_, cm_time_  );
        AddVectorElementWise(acc_moc_, moc_);
        AddVectorElementWise(acc_num_good_tracks_, num_good_tracks_ );
        AddVectorElementWise(acc_num_targets_, num_targets_ );
        AddVectorElementWise(acc_num_tracks_, num_tracks_ );
        AddVectorElementWise(acc_total_num_tracks_, total_num_tracks_ );
        AddVectorElementWise(acc_average_euclidean_err_, average_euclidean_err_ );

    }

    nvt_.clear();
    nft_.clear();
    nmt_.clear();
    nst_.clear();
    cm_time_.clear();
    moc_.clear();
    num_good_tracks_.clear();
    num_targets_.clear();
    num_tracks_.clear();
    average_euclidean_err_.clear();

    


}

//--------------------------------------------------------------------------------------------------------

template<typename tModel>
void Stats<tModel>::CalculateTotalPerformanceMeasures() {


    auto avg = [this](double& a){return a/num_monte_carl_runs_;};

    std::transform(acc_nvt_.begin(), acc_nvt_.end(), acc_nvt_.begin(), avg );
    std::transform(acc_nft_.begin(), acc_nft_.end(), acc_nft_.begin(), avg );
    std::transform(acc_nmt_.begin(), acc_nmt_.end(), acc_nmt_.begin(), avg );
    std::transform(acc_nst_.begin(), acc_nst_.end(), acc_nst_.begin(), avg );
    std::transform(acc_moc_.begin(), acc_moc_.end(), acc_moc_.begin(), avg );
    std::transform(acc_num_good_tracks_.begin(), acc_num_good_tracks_.end(), acc_num_good_tracks_.begin(), avg );
    std::transform(acc_num_targets_.begin(), acc_num_targets_.end(), acc_num_targets_.begin(), avg );
    std::transform(acc_num_tracks_.begin(), acc_num_tracks_.end(), acc_num_tracks_.begin(), avg );
    std::transform(acc_total_num_tracks_.begin(), acc_total_num_tracks_.end(), acc_total_num_tracks_.begin(), avg );
    std::transform(acc_average_euclidean_err_.begin(), acc_average_euclidean_err_.end(), acc_average_euclidean_err_.begin(), avg );

    double sum1 = std::accumulate(acc_num_targets_.begin(), acc_num_targets_.end(),0.0);
    double sum2 = std::accumulate(acc_nmt_.begin(), acc_nmt_.end(),0.0);
    double sum3 = std::accumulate(acc_average_euclidean_err_.begin(), acc_average_euclidean_err_.end(),0.0);

    tpd_ = (sum1-sum2)/sum1;
    taee_ = sum3 / acc_average_euclidean_err_.size();
    tant_ = acc_total_num_tracks_.back();

}

//--------------------------------------------------------------------------------------------------------

template<typename tModel>
void Stats<tModel>::CalculatePerformanceMeasures(const System<tModel>* sys, const std::vector<tModel>& targets, int source_index ) {
    
    // This is the first iteration in the monte carlo run. 
    if (nvt_.size() ==0) {
        num_monte_carl_runs_++;
    }

    cm_time_.push_back(sys->current_time_);
    Eigen::MatrixXd C = ComputeAssociationMatrix(sys,targets,source_index);
    Eigen::MatrixXd errs(targets.size(),1);
    errs.setZero();

    double nvt = 0;
    double nft = 0;
    double nmt = C.col(sys->good_models_.size()).sum(); // The last column represent that a target is not associated with a track
    double nst = 0;
    double moc = 0;

    for (int ii =0; ii < sys->good_models_.size(); ++ii) {
        if(C.col(ii).sum() == 1) {          // Valid track. Associated with only one target
            nvt++; 
        } else if(C.col(ii).sum() == 0) {   // False track. Not associated with any targets
            nft++;
        } else {                            // Spurious track. Associated with multiple targets. 
            nst++;
        }
    }

    Eigen::MatrixXd err;
    Meas<double> meas;
    meas.type = sys->sources_[source_index].params_.type_;
    meas.source_index = source_index;

    for (int ii = 0; ii < targets.size(); ++ii) {
        for (int jj=0; jj < sys->good_models_.size(); ++jj) {
            if(C(ii,jj)==1) {
                meas = sys->sources_[source_index].GetEstMeas(targets[ii].state_);
                err = sys->sources_[source_index].OMinus(meas, sys->sources_[source_index].GetEstMeas(sys->good_models_[jj]->state_));
                errs(ii) = err.norm();
            }       
        }
    }




    moc = nvt/targets.size();

    nvt_.push_back(nvt);
    nft_.push_back(nft);
    nmt_.push_back(nmt);
    nst_.push_back(nst);
    moc_.push_back(moc);
    num_good_tracks_.push_back(sys->good_models_.size());
    num_targets_.push_back(targets.size());
    num_tracks_.push_back(sys->models_.size());
    total_num_tracks_.push_back(sys->accumulative_number_of_tracks_);
    average_euclidean_err_.push_back(errs.norm()/errs.size());
    

    


}

//--------------------------------------------------------------------------------------------------------

template<typename tModel>
Eigen::MatrixXd Stats<tModel>::ComputeAssociationMatrix(const System<tModel>* sys, const std::vector<tModel>& targets, int source_index ) {

    if (sys->sources_.size() == 0 || source_index >= sys->sources_.size()) {
        throw std::runtime_error("Stats::ComputeAssociationMatrix: Source index out of bounds.");
    }

    Eigen::MatrixXd Cd(targets.size(), sys->good_models_.size()+1);
    Eigen::MatrixXd C(targets.size(), sys->good_models_.size()+1);
    Cd.setZero();
    C.setZero();
    Eigen::MatrixXd err;
    Eigen::MatrixXd S; // Innovation covariance
    double distance = 0; 
    double det_S = 0;
    Meas<double> meas;
    meas.type = sys->sources_[source_index].params_.type_;
    meas.source_index = source_index;

    // Associate targets to tracks using the validation region. If a track
    // is associated with a target, place the distance based on the validation
    // region in the association matrix.
    for (int ii = 0; ii < targets.size(); ++ii) {
        meas = sys->sources_[source_index].GetEstMeas(targets[ii].state_);
        for (int jj=0; jj < sys->good_models_.size(); ++jj) {
            S = sys->good_models_[jj]->GetInnovationCovariance(sys->sources_, source_index);
            err = sys->sources_[source_index].OMinus(meas, sys->sources_[source_index].GetEstMeas(sys->good_models_[jj]->state_));
            distance = (err.transpose()*S.inverse()*err)(0,0);
            if (distance <= sys->sources_[source_index].params_.gate_threshold_) {
                Cd(ii,jj) = distance;
                C(ii,jj) = 1;
            }         
        }
    }


    // We only want to assign a target to at most one track. If a 
    // target has been associated to multiple tracks, associate it with
    // the one whose distance is less than the others. 
    int min_distance_index = 0;
    double min_distance =0;
    bool min_distance_set = false;
    for (int ii = 0; ii < targets.size(); ++ii) {
        min_distance_index = sys->good_models_.size();
        min_distance_set = false;
        min_distance = 0;
        for (int jj=0; jj < sys->good_models_.size(); ++jj) {
            if (C(ii,jj) == 1) {
                if (!min_distance_set) {
                    min_distance_set = true;
                    min_distance = Cd(ii,jj);
                    min_distance_index = jj;
                } else if(Cd(ii,jj) < min_distance) {
                    min_distance = Cd(ii,jj);
                    min_distance_index = jj;
                }

                
            } 
        }
        C.row(ii).setZero();
        C(ii,min_distance_index) = 1;

    }

    return C;

}

//-------------------------------------------------------------------------------------------------------

template<typename tModel>
template<typename T>
void Stats<tModel>::AddVectorElementWise(std::vector<T>& v1, const std::vector<T>& v2) {

    for (long unsigned int ii =0; ii < v1.size(); ++ii) {
        v1[ii] += v2[ii];
    }

}

//-------------------------------------------------------------------------------------------------------
template<typename tModel>
void Stats<tModel>::WriteData(const std::string absolute_file_path) {
    std::string file_name = absolute_file_path + "/output.csv";
    std::ofstream file;
    file.open(file_name);

    auto write = [](const std::vector<double> data, const std::string name, std::ofstream& file){
        file << name << ",";
        for (long int ii = 0; ii < data.size() -1; ++ii) {
            file << data[ii] << ",";
        }
        file << data.back() << "\n";
    };
    
    if (file.is_open()) {

        write(acc_nvt_, "Number Valid Tracks", file);
        write(acc_nft_, "Number False Tracks", file);
        write(acc_nmt_, "Number Missed Targets", file);
        write(acc_nst_, "Number Spurious Tracks", file);
        write(acc_cm_time_, "Time", file);
        write(acc_moc_, "Measure of Completeness", file);
        write(acc_num_good_tracks_, "Number of Good Tracks", file);
        write(acc_num_targets_, "Number of Targets", file);
        write(acc_num_tracks_, "Number Tracks", file);
        write(acc_total_num_tracks_, "Total Number Tracks", file);
        write(acc_average_euclidean_err_, "Average Euclidean Error", file);
        file << "Track Probability of Detection " << "," << tpd_ << std::endl;
        file << "Total Average Euclidean Error " << "," << taee_ << std::endl;
        file << "Total Average number of Tracks " << "," << tant_ << std::endl;




        file.close();


    } else {
        std::cout << "Stats::WriteData: Could not open file" << std::endl;
    }

}

} // namespace

#endif //RRANSAC_STATISTICAL_INFORMATION_H_

    