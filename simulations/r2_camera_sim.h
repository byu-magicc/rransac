#pragma once


#include <stdlib.h>
#include <time.h>
#include <Eigen/Dense>
#include <string.h>
#include <chrono>
#include <vector>
#include <random>

#include "lie_groups/state.h"
#include "rransac/system.h"
#include "rransac/common/models/model_SEN_pos_vel.h"
#include "rransac/common/sources/source_SEN_pos_vel.h"
#include "rransac/common/models/model_RN.h"
#include "rransac/common/sources/source_RN.h"
#include "rransac/common/measurement/measurement_base.h"
#include "rransac/parameters.h"
#include "rransac/track_initialization/lmle_policies/linear_lmle_policy.h"
#include "rransac/track_initialization/seed_policies/null_seed_policy.h"
#include "rransac/common/transformations/trans_homography.h"
#include "rransac/common/transformations/transformation_null.h"
#include "rransac/track_initialization/ransac.h"
#include "rransac/common/data_association/validation_region_policies/validation_region_innov_policy.h"
#include "rransac/common/data_association/track_likelihood_info_policies/tli_ipdaf_policy.h"
#include "rransac/common/data_association/measurement_weight_policies/mw_ipdaf_policy.h"
#include "rransac/rransac.h"
#include "rransac/common/utilities.h"
#include "rransac/visualization/visualization_host.h"
#include "rransac/visualization/draw_meas_policies/draw_meas_R2_SE2_pos_policy.h"
#include "rransac/visualization/draw_track_policies/draw_track_policy_R2.h"
#include "rransac/statistical_information/stats.h"
#include "rransac/common/sources/source_container.h"



using namespace lie_groups;
using namespace rransac;


class CameraData {


public:

int width;                      // number of pixels along the x axis
int height;                     // number of pixels along the z axis
double fov;                     // degrees
double altitude;                // height above ground
Eigen::Matrix<double,4,4> R;              // measurement covariance matrix
// Eigen::Matrix<double,2,2> R;              // measurement covariance matrix
Eigen::Matrix3d K;              // camera matrix
double minx, maxx, miny, maxy;  // Surveillance region boundaries
double volume;
int num_false_meas;
double lambda;

// CameraData(int width, int height, double fov, double altitude, const Eigen::Matrix<double,2,2>& R, int num_false_meas) : width(width), height(height), fov(fov), altitude(altitude), R(R), num_false_meas(num_false_meas) {
CameraData(int width, int height, double fov, double altitude, const Eigen::Matrix<double,4,4>& R, int num_false_meas) : width(width), height(height), fov(fov), altitude(altitude), R(R), num_false_meas(num_false_meas) {

double f = width/(2.0*tan(fov*M_PI/180.0/2.0));
K << f, 0, width/2.0, 0, f ,height/2.0, 0,0,1;
Eigen::Matrix<double,3,1> tmp;
tmp << width, height,1;
tmp = K.inverse()*tmp;
maxx = tmp(0);
maxy = tmp(1);
tmp << 0,0,1;
tmp = K.inverse()*tmp;
minx = tmp(0);
miny = tmp(1);
volume = (maxx-minx)*(maxy-miny);
lambda = num_false_meas/volume;
if (lambda <= 0) {
    lambda = 0.1;
}

std::cout << "K: " << std::endl << K << std::endl;
std::cout << "minx: " << minx << std::endl;
std::cout << "maxx: " << maxx << std::endl;
std::cout << "miny: " << miny << std::endl;
std::cout << "maxy: " << maxy << std::endl;

}




};

//----------------------------------------------------------------------------------------------------------------------------------

class CameraSimR2 {

public:

typedef SourceSENPosVel<SE2_se2,MeasurementTypes::SEN_POS_VEL,TransformNULL> SourceSE2PosVelNull;
typedef SourceContainer<SourceSE2PosVelNull> SourceContainerSE2PosVelNull;

typedef SourceRN<R2_r2,MeasurementTypes::RN_POS_VEL,TransformNULL> SourceR2PosVelNull;
typedef SourceContainer<SourceR2PosVelNull> SourceContainerR2PosVelNull;

typedef ModelSENPosVel<SourceContainerSE2PosVelNull> TargetModel_;
typedef typename TargetModel_::State TargetState_;
typedef typename TargetState_::Algebra TargetAlgebra_;

typedef ModelRN<SourceContainerR2PosVelNull> TrackingModel_;
typedef typename TrackingModel_::Transformation Transformation_;
typedef typename TrackingModel_::Transformation::TransformDataType TransformMatData_;
typedef typename TrackingModel_::State TrackingState_;
typedef typename TrackingState_::Algebra TrackingAlgebra_;
typedef SourceSE2PosVelNull TargetSource_;
typedef RRANSACTemplateParameters<SourceContainerR2PosVelNull,ModelRN,NULLSeedPolicy,LinearLMLEPolicy,ValidationRegionInnovPolicy, TLI_IPDAFPolicy, MW_IPDAFPolicy> RRANSACParameters;
typedef RRANSAC<RRANSACParameters> RRANSAC_;
typedef typename RRANSACParameters::_Ransac RANSAC_;
typedef Eigen::Matrix<double,4,4> MatR_;
static constexpr MeasurementTypes MeasurementType= MeasurementTypes::RN_POS_VEL;
typedef Eigen::Matrix<double,4,4> ProcessNoiseCov_;
typedef Eigen::Matrix<double,2,1> VecU_;
typedef typename TargetModel_::Base::Measurement Measurement;




CameraData camera_data_;
std::vector<TargetModel_> tracks_;
double dt_;
double end_time_;
double start_time_;
std::vector<std::vector<double>> e_history_;
TransformMatData_ t_data_;
Transformation_ transformation_;
static constexpr bool transform_data_ = false;
RRANSAC_ rransac_;
const System<TrackingModel_>* sys_;
std::vector<TargetSource_> sources_;
Measurement m_;
Eigen::Matrix<double,TargetAlgebra_::dim_,TargetAlgebra_::dim_> noise_mat_;
double process_noise_;
double meas_noise_;
VisualizationHost<TrackingModel_, DrawMeasR2SE2PosPolicy, DrawTrackPolicyR2> viz_;
std::default_random_engine gen_;

CameraSimR2(const CameraData& camera_data, double dt, double end_time, int num_tracks);

TargetState_ GenerateRandomState();

bool InsurveillanceRegion(const TargetState_& state);

void Propagate(double start_time, double end_time, Stats<TrackingModel_>& stats);

};