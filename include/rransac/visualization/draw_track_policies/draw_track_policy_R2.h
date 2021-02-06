#ifndef RRANSAC_VISUALIZATION_DRAW_TRACK_POLICIES_DRAW_SE2_POLICY_H_
#define RRANSAC_VISUALIZATION_DRAW_TRACK_POLICIES_DRAW_SE2_POLICY_H_


#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <cmath>

#include "visualization/draw_info.h"
#include "system.h"

namespace rransac {


template<typename tModel>
class DrawTrackPolicyR2 {

public:

static void DrawTrackPolicy(cv::Mat& img, const tModel& model, const System<tModel>* sys, const DrawInfo& draw_info);


};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                        Definitions
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename tModel>
void DrawTrackPolicyR2<tModel>::DrawTrackPolicy(cv::Mat& img, const tModel& model, const System<tModel>* sys, const DrawInfo& draw_info) {

    // The negations on some values are need in order to transform the tracking frame to the frame for drawing the image.
    int radius = 5;
    Eigen::Matrix<double,2,1> pos = model.state_.g_.data_*draw_info.scale_drawing;
    Eigen::Matrix<double,2,1> vel = model.state_.u_.data_*draw_info.scale_drawing*draw_info.scale_draw_vel;

    cv::Point p(pos(0),-pos(1));
    p += draw_info.img_center;
    cv::Point v(vel(0),-vel(1));

    cv::circle(img,p,radius*draw_info.scale_draw_pos*draw_info.scale_drawing,draw_info.color_pos,-1,cv::LINE_AA);

    cv::arrowedLine(img,p,p+v,draw_info.color_vel,draw_info.line_thickness, cv::LINE_AA);

    

    // // Draw validation region for position only
    // if (draw_info.draw_validation_region) {
    //     for (auto& source : sys->sources_) {
    //         if (source.params_.type_ == MeasurementTypes::RN_POS) {
    //             unsigned int source_index = source.params_.source_index_;
    //             Eigen::MatrixXd S = model.GetInnovationCovariance(sys->sources_,source_index);
    //             Eigen::EigenSolver<Eigen::Matrix2d> eigen_solver;
    //             eigen_solver.compute(S);
    //             Eigen::Vector2cd eigen_values = eigen_solver.eigenvalues();
    //             Eigen::Matrix2cd eigen_vectors = eigen_solver.eigenvectors();
    //             double th = 0;
    //             // Make sure that the x component is positive
    //             if (std::real(eigen_vectors(0,0)) < 0) {
    //                 eigen_vectors.block(0,0,2,1)*=-1;
    //             }
    //             if (std::real(eigen_vectors(0,1)) < 0) {
    //                 eigen_vectors.block(0,1,2,1)*=-1;
    //             }


    //             double scale = std::sqrt(source.params_.gate_threshold_)*draw_info.scale_drawing;
    //             if (std::real(eigen_vectors(0,0)*eigen_vectors(1,0)) < 0) {
    //                 th = std::real(eigen_vectors(0,0))*180/M_PI;
    //                 cv::Size size(std::sqrt(std::norm(eigen_values(0)))*scale,std::sqrt(std::norm(eigen_values(1)))*scale);
    //                 cv::ellipse(img,p,size,th, 0,360,draw_info.color_vel,draw_info.line_thickness, cv::LINE_AA);
    //             } else {
    //                 th = std::real(eigen_vectors(0,1))*180/M_PI;
    //                 cv::Size size(std::sqrt(std::norm(eigen_values(1)))*scale,std::sqrt(std::norm(eigen_values(0)))*scale);
    //                 cv::ellipse(img,p,size,th, 0,360,draw_info.color_vel,draw_info.line_thickness, cv::LINE_AA);
    //             }
    //             break;
    //         }

    //     }
    // }



}


} // namespace rransac

#endif //RRANSAC_VISUALIZATION_DRAW_TRACK_POLICIES_DRAW_SE2_POLICY_H_