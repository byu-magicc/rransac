#ifndef RRANSAC_VISUALIZATION_DRAW_TRACK_POLICIES_DRAW_POLICY_SE2_H_
#define RRANSAC_VISUALIZATION_DRAW_TRACK_POLICIES_DRAW_POLICY_SE2_H_


#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <cmath>

#include "lie_groups/lie_algebras/so2.h"
#include "rransac/visualization/draw_info.h"
#include "rransac/system.h"

namespace rransac {


template<typename tModel>
class DrawTrackPolicySE2 {

public:

static void DrawTrackPolicy(cv::Mat& img, const tModel& model, const System<tModel>* sys, const DrawInfo& draw_info);


};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                        Definitions
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename tModel>
void DrawTrackPolicySE2<tModel>::DrawTrackPolicy(cv::Mat& img, const tModel& model, const System<tModel>* sys, const DrawInfo& draw_info) {

    // The negations on some values are need in order to transform the tracking frame to the frame for drawing the image.
    
    int num_points = 4;
    cv::Size radius(2*draw_info.scale_draw_pos,2*draw_info.scale_draw_pos);
    cv::Point points[1][num_points];
    Eigen::Matrix<double,3,1> pt;
    double body_length = 5;

    // Draw body
    pt << body_length*2, 0, 0;
    pt*= draw_info.scale_draw_pos;
    pt(2) = draw_info.scale_drawing;
    pt = model.state_.g_.data_*pt;
    pt(1) = -pt(1);
    points[0][0] = cv::Point(pt(0),pt(1))+draw_info.img_center;
    
    pt << -body_length, -body_length, 0;
    pt*= draw_info.scale_draw_pos;
    pt(2) = draw_info.scale_drawing;
    pt = model.state_.g_.data_*pt;
    pt(1) = -pt(1);
    points[0][1] = cv::Point(pt(0),pt(1))+draw_info.img_center;

    pt << 0, 0, 0;
    pt*= draw_info.scale_draw_pos;
    pt(2) = draw_info.scale_drawing;
    pt = model.state_.g_.data_*pt;
    pt(1) = -pt(1);
    points[0][2] = cv::Point(pt(0),pt(1))+draw_info.img_center;

    pt << -body_length, body_length, 0;
    pt*= draw_info.scale_draw_pos;
    pt(2) = draw_info.scale_drawing;
    pt = model.state_.g_.data_*pt;
    pt(1) = -pt(1);
    points[0][3] = cv::Point(pt(0),pt(1))+draw_info.img_center;

    const cv::Point* ppt[1] = {points[0]};
    cv::fillPoly(img,ppt,&num_points,1, draw_info.color_pos,cv::LINE_AA);

    // Draw translational and angular velocity
    Eigen::Matrix<double,2,1> vel;
    vel = model.state_.u_.p_;
    vel = model.state_.g_.R_*vel*draw_info.scale_draw_vel*draw_info.scale_drawing;
    cv::Point vel_point;
    vel_point.x = vel(0);
    vel_point.y = -vel(1);

    double start_angle = -lie_groups::so2<double>::Log(model.state_.g_.R_)(0)*180/M_PI;
    double end_angle = start_angle - model.state_.u_.th_(0)*180/M_PI*draw_info.scale_draw_vel;

    if (fabs(start_angle-end_angle)>=300) {
        end_angle = start_angle - model.state_.u_.th_(0)*300/M_PI/fabs(model.state_.u_.th_(0));
    }

    cv::arrowedLine(img, points[0][2],points[0][2]+vel_point,draw_info.color_vel,draw_info.line_thickness, cv::LINE_AA);

    cv::ellipse(img,points[0][2],radius,0,start_angle, end_angle,draw_info.color_vel,draw_info.line_thickness, cv::LINE_AA);

    // Draw validation region for position only
    if (draw_info.draw_validation_region) {

        unsigned int source_index = -1;
        for (auto& source : sys->sources_) {
            if (source.params_.type_ == MeasurementTypes::SEN_POS) {
                source_index = source.params_.source_index_;
                break;
            }
        }
        if (source_index < 0);
            source_index =0;

        auto& source = sys->sources_[source_index];
        Eigen::MatrixXd S = model.GetInnovationCovariance(sys->sources_,source_index);
        Eigen::Matrix2d S_pos = S.block(0,0,2,2)*source.params_.gate_threshold_;
        Eigen::EigenSolver<Eigen::Matrix2d> eigen_solver;
        eigen_solver.compute(S_pos);
        Eigen::Vector2cd eigen_values = eigen_solver.eigenvalues();
        Eigen::Matrix2cd eigen_vectors = eigen_solver.eigenvectors();
        double th = 0;
        // Make sure that the x component is positive
        if (std::real(eigen_vectors(0,0)) < 0) {
            eigen_vectors.block(0,0,2,1)*=-1;
        }
        if (std::real(eigen_vectors(0,1)) < 0) {
            eigen_vectors.block(0,1,2,1)*=-1;
        }


        double scale = draw_info.scale_drawing;
        if (std::real(eigen_vectors(0,0)*eigen_vectors(1,0)) < 0) {
            th = std::real(eigen_vectors(0,0))*180/M_PI;
            cv::Size size(std::sqrt(std::norm(eigen_values(0)))*scale,std::sqrt(std::norm(eigen_values(1)))*scale);
            cv::ellipse(img,points[0][2],size,th, 0,360,draw_info.color_vel,draw_info.line_thickness, cv::LINE_AA);
        } else {
            th = std::real(eigen_vectors(0,1))*180/M_PI;
            cv::Size size(std::sqrt(std::norm(eigen_values(1)))*scale,std::sqrt(std::norm(eigen_values(0)))*scale);
            cv::ellipse(img,points[0][2],size,th, 0,360,draw_info.color_vel,draw_info.line_thickness, cv::LINE_AA);
        }
    }



}


} // namespace rransac

#endif //RRANSAC_VISUALIZATION_DRAW_TRACK_POLICIES_DRAW_POLICY_SE2_H_