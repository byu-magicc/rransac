#ifndef RRANSAC_VISUALIZATION_DRAW_TRACK_POLICIES_DRAW_SE2_POLICY_H_
#define RRANSAC_VISUALIZATION_DRAW_TRACK_POLICIES_DRAW_SE2_POLICY_H_


#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <vector>
#include <Eigen/Core>
#include "visualization/draw_info.h"
#include "lie_algebras/so2.h"
#include <cmath>

namespace rransac {

template <typename tModel>
class DrawSE2Policy {

public:

static void DrawTrackPolicy(cv::Mat& img, const tModel& model, const DrawInfo& draw_info);


private:

std::vector<cv::Point> points_;

};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                        Definitions
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename tModel>
void DrawSE2Policy<tModel>::DrawTrackPolicy(cv::Mat& img, const tModel& model, const DrawInfo& draw_info) {

    // The negations on some values are need in order to transform the tracking frame to the frame for drawing the image.
    
    int num_points = 4;
    cv::Size radius(1*draw_info.scale_pos,1*draw_info.scale_pos);
    cv::Point points[1][num_points];
    Eigen::Matrix<double,3,1> pt;
    double body_length = 5;

    // Draw body
    pt << body_length*2, 0, 0;
    pt*= draw_info.scale_pos;
    pt(2) = 1;
    pt = model.state_.g_.data_*pt;
    pt(1) = -pt(1);
    points[0][0] = cv::Point(pt(0),pt(1))+draw_info.img_center;
    
    pt << -body_length, -body_length, 0;
    pt*= draw_info.scale_pos;
    pt(2) = 1;
    pt = model.state_.g_.data_*pt;
    pt(1) = -pt(1);
    points[0][1] = cv::Point(pt(0),pt(1))+draw_info.img_center;

    pt << 0, 0, 0;
    pt*= draw_info.scale_pos;
    pt(2) = 1;
    pt = model.state_.g_.data_*pt;
    pt(1) = -pt(1);
    points[0][2] = cv::Point(pt(0),pt(1))+draw_info.img_center;

    pt << -body_length, body_length, 0;
    pt*= draw_info.scale_pos;
    pt(2) = 1;
    pt = model.state_.g_.data_*pt;
    pt(1) = -pt(1);
    points[0][3] = cv::Point(pt(0),pt(1))+draw_info.img_center;

    const cv::Point* ppt[1] = {points[0]};
    cv::fillPoly(img,ppt,&num_points,1, draw_info.color_pos,cv::LINE_AA);

    // Draw Velocity
    Eigen::Matrix<double,2,1> vel;
    vel = model.state_.u_.p_;
    vel = model.state_.g_.R_*vel*draw_info.scale_vel;
    cv::Point vel_point;
    vel_point.x = vel(0);
    vel_point.y = -vel(1);

    double start_angle = -lie_groups::so2<double>::Log(model.state_.g_.R_)(0)*180/M_PI;
    double end_angle = start_angle - model.state_.u_.th_(0)*180/M_PI*draw_info.scale_vel;

    if (fabs(start_angle-end_angle)>=300) {
        end_angle = start_angle - model.state_.u_.th_(0)*300/M_PI/fabs(model.state_.u_.th_(0));
    }

    cv::arrowedLine(img, points[0][2],vel_point+draw_info.img_center,draw_info.color_vel,draw_info.line_thickness, cv::LINE_AA);

    cv::ellipse(img,points[0][2],radius,0,start_angle, end_angle,draw_info.color_vel,draw_info.line_thickness, cv::LINE_AA);


}


} // namespace rransac

#endif //RRANSAC_VISUALIZATION_DRAW_TRACK_POLICIES_DRAW_SE2_POLICY_H_