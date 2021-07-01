#ifndef RRANSAC_VISUALIZATION_DRAW_TRACK_POLICIES_DRAW_SE2_POLICY_H_
#define RRANSAC_VISUALIZATION_DRAW_TRACK_POLICIES_DRAW_SE2_POLICY_H_


#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <cmath>

#include "rransac/visualization/draw_info.h"
#include "rransac/system.h"

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

    bool transform_state = false;
    Eigen::MatrixXd EmptyMat;

    // The negations on some values are need in order to transform the tracking frame to the frame for drawing the image.
    int radius = 2;
    Eigen::Matrix<double,2,1> pos = model.state_.g_.data_*draw_info.scale_drawing;
    Eigen::Matrix<double,2,1> vel = model.state_.u_.data_*draw_info.scale_drawing*draw_info.scale_draw_vel;

    cv::Point p(pos(0),draw_info.flip_y*pos(1));
    p += draw_info.img_center;

    if (p.x <= img.cols && p.x >= 0 && p.y <= img.rows && p.y>=0){ 

        cv::Point v(vel(0),draw_info.flip_y*vel(1));

        cv::circle(img,p,radius*draw_info.scale_draw_pos,draw_info.color_pos,-1,cv::LINE_AA);

        cv::arrowedLine(img,p,p+v,draw_info.color_vel,draw_info.line_thickness, cv::LINE_AA);

    

        // Draw validation region for the first source
        if (draw_info.draw_validation_region) {

            unsigned int source_index = -1;
            for (unsigned int ii = 0; ii < sys->source_container_.num_sources_; ++ii) {
                if (sys->source_container_.GetParams(ii).type_ == MeasurementTypes::RN_POS) {
                    source_index = ii;
                    break;
                }
            }
            if (source_index < 0);
                source_index =0;

#if RRANSAC_VIZ_HOOKS
            Eigen::MatrixXd S;
            if (source_index > model.S_validation_.size()) {
                std::cerr << "DrawTrackPolicySE2::DrawTrackPolicy source index is larger than the than the model's member variable size S_validation_" << std::endl;
                std::cerr << "source index: " << source_index << std::endl;
                std::cerr << "Size of S_validation_: " << model.S_validation_.size() << std::endl;
                S = model.S_validation_[source_index];
            }
            else {
                S = model.GetInnovationCovariance(sys->source_container_,source_index);
            }
            
#else
            Eigen::MatrixXd S = model.GetInnovationCovariance(sys->source_container_,source_index,transform_state,EmptyMat);
#endif
            Eigen::Matrix2d S_pos = S.block(0,0,2,2)*sys->source_container_.GetParams(source_index).gate_threshold_;
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

            double scale = static_cast<double>(draw_info.scale_drawing);
            if (std::real(eigen_vectors(0,0)*eigen_vectors(1,0)) < 0) {
                th = std::real(eigen_vectors(0,0))*180/M_PI;
                cv::Size size(scale*std::sqrt(std::fabs(eigen_values(0).real())),scale*std::sqrt(std::fabs(eigen_values(1).real())));
                std::cerr << "size: " << size << std::endl;
                cv::ellipse(img,p,size,th, 0,360,draw_info.color_vel,draw_info.line_thickness, cv::LINE_AA);
            } else {
                th = std::real(eigen_vectors(0,1))*180/M_PI;
                cv::Size size(scale*std::sqrt(std::fabs(eigen_values(1).real())),scale*std::sqrt(std::fabs(eigen_values(0).real())));
                std::cerr << "size: " << size << std::endl;
                cv::ellipse(img,p,size,th, 0,360,draw_info.color_vel,draw_info.line_thickness, cv::LINE_AA);
            }
        }
    }


}


} // namespace rransac

#endif //RRANSAC_VISUALIZATION_DRAW_TRACK_POLICIES_DRAW_SE2_POLICY_H_