
#include <opencv2/core.hpp>
#include <opencv2/plot.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <vector>
#include <Eigen/Core>

#include "common/models/model_SEN_pos_vel.h"
#include "common/transformations/transformation_null.h"
#include "state.h"
#include "visualization/draw_track_policies/draw_track_SE2_policy.h"
#include "visualization/draw_meas_policies/draw_meas_R2_SE2_pos_policy.h"
#include "lie_algebras/so2.h"

#include "visualization/draw_info.h"
#include "system.h"
#include "common/measurement/measurement_base.h"
#include "common/sources/source_base.h"

cv::Scalar ScalarHSV2BGR(const cv::Scalar& scalar) {
    cv::Mat bgr;
    cv::Mat hsv(1,1, CV_8UC3, scalar);
    cvtColor(hsv, bgr, cv::COLOR_HSV2BGR);
    return cv::Scalar(bgr.data[0], bgr.data[1], bgr.data[2]);
}

cv::Scalar ScalarBGR2HSV(const cv::Scalar& scalar) {
    cv::Mat bgr(1,1,CV_8UC3,scalar);
    cv::Mat hsv;
    cv::cvtColor(bgr,hsv,cv::COLOR_BGR2HSV);
    return cv::Scalar(hsv.data[0], hsv.data[1], hsv.data[2]);
}

int main() {

typedef rransac::ModelSENPosVel<lie_groups::SE2_se2,rransac::TransformNULL> Model;
rransac::System<Model> sys;
typename Model::Source source;
rransac::SourceParameters source_params;
Eigen::Matrix<double,1,1> rot;
rot << 0.45;
Eigen::Matrix2d R = lie_groups::so2<double>::Exp(rot);
source_params.gate_probability_ = 0.95;
source_params.meas_cov_ = Eigen::Matrix2d::Identity();
source_params.meas_cov_ << 0.1, 0, 0, 0.1;
source_params.meas_cov_ = R*source_params.meas_cov_*R.transpose();
std::cerr << "source_params.meas_cov_: " << std::endl << source_params.meas_cov_ << std::endl;

source_params.type_ = rransac::MeasurementTypes::SEN_POS;
source_params.source_index_ = 0;
source.Init(source_params);
sys.sources_.push_back(source);

sys.params_.cluster_position_threshold_ = 5;
sys.params_.cluster_velocity_threshold_ = 4;
sys.params_.process_noise_covariance_ = Eigen::Matrix<double,5,5>::Identity()*0.01;
// sys.params_.process_noise_covariance_ = rot*source_params.meas_cov_*rot.transpose();

unsigned int w = 500;
cv::Mat img = cv::Mat( w, w, CV_8UC3 , cv::Scalar(255,255,255));

std::vector<rransac::Meas<double>> measurements;
rransac::Meas<double> m;
m.pose = Eigen::Matrix<double,2,1>::Zero();
m.twist = m.pose;

m.pose << 50,50;
m.twist << 1,5;

measurements.push_back(m);
m.pose << -50, 50;
m.twist << -1,5;
measurements.push_back(m);
m.pose << -50,-50;
m.twist << -1,-5;
measurements.push_back(m);
m.pose << 50, -50;
m.twist << 1,-5;
measurements.push_back(m);


rransac::DrawInfo draw_info;
draw_info.img_center = cv::Point(w/2.0, w/2.0);
draw_info.color_pos = cv::Scalar(0,0,255);
cv::Scalar covert = ScalarBGR2HSV(draw_info.color_pos);
covert[1]*=0.5;
draw_info.color_pos = ScalarHSV2BGR(covert);
// draw_info.color_pos[3]=0.1;
draw_info.color_vel = cv::Scalar(0,0,255,1);
draw_info.scale_draw_pos = 2;
draw_info.scale_draw_vel = 5;
draw_info.line_thickness = 1;
draw_info.scale_drawing = 3;


Model track;
track.Init(sys.params_);
track.err_cov_*=0.01;

rot << 3.14159;
track.state_ = Model::State::Random();
track.state_.g_.R_ = lie_groups::so2<double>::Exp(rot/4.0);
track.state_.g_.t_*=20;
track.state_.u_.p_*=10;
track.state_.u_.th_*=-1;

std::cout << "state g: " << std::endl << track.state_.g_.data_ << std::endl;
std::cout << "state u: " << std::endl << track.state_.u_.data_ << std::endl;


rransac::DrawSE2Policy<Model>::DrawTrackPolicy(img,track,sys,draw_info);
for (auto& meas : measurements) 
    rransac::DrawMeasR2SE2PosPolicy<Model>::DrawMeasPolicy(img,meas,sys,draw_info);


// cv::resize(img,img,cv::Size(1000,1000));

cv::imshow("plot",img);
cv::waitKey(0);

    return 0;
}