
#include "common/models/model_RN.h"
#include "common/sources/source_base.h"
#include "common/sources/source_RN.h"
#include "common/sources/source_SEN_pose_twist.h"
#include "state.h"
#include "common/transformations/transformation_null.h"
#include "common/measurement/measurement_base.h"
#include "system.h"
#include "common/transformations/transformation_null.h"

#include "visualization/draw_info.h"
#include "visualization/draw_track_policies/draw_track_policy_R2.h"





int main() {

    typedef rransac::ModelRN<lie_groups::R2_r2,rransac::TransformNULL,rransac::SourceRN> Model;
    rransac::ModelRN<lie_groups::R2_r2,rransac::TransformNULL,rransac::SourceRN> model;
    cv::Mat img( 500,500, CV_8UC3 , cv::Scalar(255,255,255));
    model.state_.g_.data_ << 50,50;
    model.state_.u_.data_ << -1,2;
    model.err_cov_.setIdentity();
    
    rransac::SourceParameters source_params;
    source_params.source_index_ = 0;
    source_params.type_ = rransac::MeasurementTypes::RN_POS_VEL;
    source_params.meas_cov_ = Eigen::Matrix4d::Identity()*2;
    
    rransac::SourceRN<lie_groups::R2_r2> source;
    source.Init(source_params);
    rransac::System<Model> sys;
    sys.sources_.push_back(source);

    rransac::DrawInfo draw_info;

    draw_info.color_pos = cv::Scalar(255,0,0);
    draw_info.color_vel = cv::Scalar(0,255,0);
    draw_info.draw_validation_region = true;
    draw_info.scale_draw_pos = 1;
    draw_info.scale_draw_vel = 2;
    draw_info.scale_drawing = 3;
    draw_info.line_thickness = 2;
    draw_info.img_center = cv::Point(250,250);

    rransac::DrawTrackPolicyR2<Model>::DrawTrackPolicy(img, model, &sys, draw_info);

    cv::imshow("blah", img);
    cv::waitKey(0);





    return 0;
}