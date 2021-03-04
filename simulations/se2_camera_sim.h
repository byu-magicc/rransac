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
#include "rransac/common/measurement/measurement_base.h"
#include "rransac/parameters.h"
#include "rransac/track_initialization/seed_policies/SE2_pos_seed_policy.h"
#include "rransac/track_initialization/lmle_policies/nonlinear_lmle_policy.h"
#include "rransac/common/transformations/trans_homography.h"
#include "rransac/track_initialization/ransac.h"
#include "rransac/common/data_association/model_policies/model_pdf_policy.h"
#include "rransac/common/data_association/cluster_data_tree_policies/data_tree_cluster_association_policy.h"
#include "rransac/rransac.h"
#include "rransac/common/utilities.h"
#include "rransac/visualization/visualization_host.h"
#include "rransac/visualization/draw_meas_policies/draw_meas_R2_SE2_pos_policy.h"
#include "rransac/visualization/draw_track_policies/draw_track_policy_SE2.h"



using namespace lie_groups;
using namespace rransac;


class CameraData {


public:

int width;                      // number of pixels along the x axis
int height;                     // number of pixels along the z axis
double fov;                     // degrees
double altitude;                // height above ground
// Eigen::Matrix<double,4,4> R;              // measurement covariance matrix
Eigen::Matrix<double,2,2> R;              // measurement covariance matrix
Eigen::Matrix3d K;              // camera matrix
double minx, maxx, miny, maxy;  // Surveillance region boundaries
double volume;
int num_false_meas;
double lambda;

CameraData(int width, int height, double fov, double altitude, const Eigen::Matrix<double,2,2>& R, int num_false_meas) : width(width), height(height), fov(fov), altitude(altitude), R(R), num_false_meas(num_false_meas) {
// CameraData(int width, int height, double fov, double altitude, const Eigen::Matrix<double,4,4>& R, int num_false_meas) : width(width), height(height), fov(fov), altitude(altitude), R(R), num_false_meas(num_false_meas) {

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

class CameraSimSE2 {

typedef ModelSENPosVel<SE2_se2, TransformHomography,SourceSENPosVel> Model_;
typedef typename Model_::Transformation Transformation_;
typedef typename Model_::Transformation::MatData TransformMatData_;
typedef typename Model_::State State_;
typedef typename State_::Algebra Algebra_;
typedef typename Model_::Source Source_;
typedef Ransac<Model_, SE2PosSeedPolicy, NonLinearLMLEPolicy, ModelPDFPolicy> RANSAC_;
typedef RRANSACTemplateParameters<SE2_se2,SourceSENPosVel,TransformHomography,ModelSENPosVel,SE2PosSeedPolicy,NonLinearLMLEPolicy,ModelPDFPolicy,DataTreeClusterAssociationPolicy> RRANSACParameters;
typedef RRANSAC<RRANSACParameters> RRANSAC_;
// typedef Eigen::Matrix<double,4,4> MatR_;
typedef Eigen::Matrix<double,2,2> MatR_;
static constexpr MeasurementTypes MeasurementType= MeasurementTypes::SEN_POS;
// static constexpr MeasurementTypes MeasurementType= MeasurementTypes::SEN_POS_VEL;
typedef Eigen::Matrix<double,5,5> ProcessNoiseCov_;
typedef Eigen::Matrix<double,3,1> VecU_;


public:

CameraData camera_data_;
std::vector<Model_> tracks_;
double dt_;
double end_time_;
double start_time_;
std::vector<std::vector<double>> e_history_;
TransformMatData_ t_data_;
Transformation_ transformation_;
static constexpr bool transform_data_ = false;
RRANSAC_ rransac_;
const System<Model_>* sys_;
std::vector<Source_> sources_;
Meas<double> m_;
Eigen::Matrix<double,Algebra_::dim_,Algebra_::dim_> noise_mat_;
double process_noise_;
double meas_noise_;
VisualizationHost<Model_, DrawMeasR2SE2PosPolicy, DrawTrackPolicySE2> viz_;
std::default_random_engine gen_;

CameraSimSE2(const CameraData& camera_data, double dt, double end_time, int num_tracks);

State_ GenerateRandomState();

bool InsurveillanceRegion(const State_& state);

void Propagate(double start_time, double end_time);

};