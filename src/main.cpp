
#include <opencv2/core.hpp>
#include <opencv2/plot.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <vector>

#include "common/models/model_SEN_pose_twist.h"
#include "common/transformations/transformation_null.h"
#include "state.h"
#include "visualization/draw_track_policies/draw_SE2_policy.h"
#include "lie_algebras/so2.h"

#include "visualization/draw_info.h"

int main() {

unsigned int w = 500;
cv::Mat img = cv::Mat( w, w, CV_8UC3 , cv::Scalar(255,255,255));

rransac::DrawInfo draw_info;
draw_info.img_center = cv::Point(w/2.0, w/2.0);
draw_info.color_pos = cv::Scalar(255,0,0);
draw_info.color_vel = cv::Scalar(0,0,255);
draw_info.scale_pos = 5;
draw_info.scale_vel = 5;
draw_info.line_thickness = 1;

typedef rransac::ModelSENPoseTwist<lie_groups::SE2_se2,rransac::TransformNULL> Model;
Model track;
Eigen::Matrix<double,1,1> rot;
rot << 3.14159;
track.state_ = Model::State::Random();
track.state_.g_.R_ = lie_groups::so2<double>::Exp(rot/4.0);
track.state_.g_.t_*=20;
track.state_.u_.p_*=10;
track.state_.u_.th_*=-1;

std::cout << "state g: " << std::endl << track.state_.g_.data_ << std::endl;
std::cout << "state u: " << std::endl << track.state_.u_.data_ << std::endl;


rransac::DrawSE2Policy<Model>::DrawTrackPolicy(img,track,draw_info);

// cv::resize(img,img,cv::Size(1000,1000));

cv::imshow("plot",img);
cv::waitKey(0);

    return 0;
}