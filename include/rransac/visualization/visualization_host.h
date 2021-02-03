#ifndef RRANSAC_VISUALIZATION_VISUALIZATION_HOST_INFO_H_
#define RRANSAC_VISUALIZATION_VISUALIZATION_HOST_INFO_H_

#include "system.h"
#include "draw_info.h"
#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <vector>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <random>

namespace rransac
{
    
template<typename tModel, template<typename > typename tDrawMeasurementPolicy, template<typename> typename tDrawTrackPolicy>
class VisualizationHost : tDrawMeasurementPolicy<tModel>, tDrawTrackPolicy<tModel> {

public:

VisualizationHost(const std::vector<int>& img_dimensions, const double drawing_scale);

void DrawSystem(const System<tModel>* sys);

void DrawSystemAndTrueTracks(const System<tModel>* sys, const std::vector<tModel>& tracks);

const cv::Mat* GetImg();

private:

std::vector<int> img_dimensions_;
double drawing_scale_(drawing_scale_);
cv::Mat img_;
cv::Scalar true_track_color_;
cv::Scalar estimated_poor_track_color_;
cv::Scalar estimated_good_track_color_;
cv::Scalar velocity_color_;
cv::Scalar measurement_color_;
std::vector<cv::Scalar> cluster_colors_;
std::default_random_engine rng_;
std::uniform_real_distribution<double> dist_;
DrawInfo draw_info_;
bool reset_img_;


static void DrawMeas(cv::Mat& img, const Meas<double>& meas, const System<tModel>& sys,  const DrawInfo& draw_info) {
    VisualizationHost::DrawMeasPolicy(img,meas,sys,draw_info);
}

static void DrawTrack(cv::Mat& img, const tModel& model, const System<tModel>& sys, const DrawInfo& draw_info) {
    VisualizationHost::DrawTrackPolicy(img,model,sys,draw_info);
} 

static cv::Scalar GetRandomColor();

static cv::Scalar ScalarHSV2BGR(const cv::Scalar& scalar) {
    cv::Mat bgr;
    cv::Mat hsv(1,1, CV_8UC3, scalar);
    cvtColor(hsv, bgr, cv::COLOR_HSV2BGR);
    return cv::Scalar(bgr.data[0], bgr.data[1], bgr.data[2]);
}

static cv::Scalar ScalarBGR2HSV(const cv::Scalar& scalar) {
    cv::Mat bgr(1,1,CV_8UC3,scalar);
    cv::Mat hsv;
    cv::cvtColor(bgr,hsv,cv::COLOR_BGR2HSV);
    return cv::Scalar(hsv.data[0], hsv.data[1], hsv.data[2]);
}

static cv::Scalar ScalarFadeColor(const cv::Scalar& scalar, double time_stamp, double curr_time, double time_window) {

    cv::Scalar faded_color = ScalarBGR2HSV(scalar);
    double scale = 0.8/time_window*(time_stamp-curr_time+time_window) + 0.2;
    faded_color[1]*=scale;
    return ScalarHSV2BGR(faded_color);


}

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                        Definitions
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename tModel, template<typename > typename tDrawMeasurementPolicy, template<typename> typename tDrawTrackPolicy>

VisualizationHost<tModel,tDrawMeasurementPolicy,tDrawTrackPolicy>::VisualizationHost(const std::vector<int>& img_dimensions, const double drawing_scale) :img_dimensions_(img_dimensions), drawing_scale_(drawing_scale) { 
    srand(time(NULL));
    dist_ = std::uniform_real_distribution<int>(0,1);
    true_track_color_ = cv::Scalar(0,255,0);
    estimated_poor_track_color_ = cv::Scalar(204,102,0);
    estimated_good_track_color_ = cv::Scalar(0,255,255);
    velocity_color_ = cv::Scalar(0,0,255);
    measurement_color_ = cv::Scalar(255,0,127);
    reset_img_ = true;

    draw_info_.img_center = cv::Point(img_dimensions_(0)/2,img_dimensions_(1)/2);



// cv::Scalar color_pos;
// cv::Scalar color_vel;
// double scale_draw_pos=1;   /** < When the object being drawn is a filled object, this value will scale the volume of the drawn object */
// double scale_draw_vel=1;   /** < This scales the size of the velocity vectors that are being drawn */
// double scale_drawing=1;    /** < Scales the position of objects being drawn by this value. It causes the image to expand or contract */
// cv::Point img_center;
// int line_thickness = 1;
// bool draw_validation_region = false;


    
}

//----------------------------------------------------------------------------------------------------------------------

template<typename tModel, template<typename > typename tDrawMeasurementPolicy, template<typename> typename tDrawTrackPolicy>
void VisualizationHost<tModel,tDrawMeasurementPolicy,tDrawTrackPolicy>::DrawSystem(const System<tModel>* sys) {


    if (reset_img_)
        img_ = cv::Mat( img_dimensions_(0), img_dimensions_(1), CV_8UC3 , cv::Scalar(255,255,255));


    // Draw measurements
    cv::Scalar color;
    cv::Scalar faded_color;
    
    for (auto cluster_iterator = sys->data_tree_.data_.begin(); cluster_iterator != sys->data_tree_.data_.end(); ++cluster_iterator) {

        if (cluster_iterator->cluster_label_ >= 0) {
            if (cluster_iterator->cluster_label_  >= cluster_colors_.size())
                cluster_color_.push_back(GetRandomColor());
            color = cluster_iterator[cluster_iterator->cluster_label_ ];
        } else {
            color = measurement_color_;
        }

        for(auto outer_iter = cluster_iterator->data_.begin(); outer_iter != cluster_iterator->data_.end(); ++ outer_iter) {
            faded_color = ScalarFadeColor(color, outer_iter->front().time_stamp, sys->current_time_, sys->params_.meas_time_window_);
            draw_info_.color_pos = faded_color;
            draw_info_.color_vel = ScalarFadeColor(velocity_color_, outer_iter->front().time_stamp, sys->current_time_, sys->params_.meas_time_window_);


            for (auto inner_iter = outer_iter->begin(); inner_iter != outer_iter->end(); ++ inner_iter) {
                DrawMeas(img_, *inner_iter, sys, draw_info_)
            }
        }
    }

    // Draw tracks

    for (auto model_iterator = sys->models_.begin(); model_iterator!= sys->models_.end(); ++ model_iterator) {
        if (model_iterator->label_ >0)
            color = estimated_good_track_color_;
        else 
            color = estimated_poor_track_color_;

        faded_color = ScalarFadeColor(color, sys->current_time_-model_iterator->missed_detection_time_, sys->current_time_, sys->params_.max_missed_detection_time_);
        draw_info_.color_pos = faded_color;
        draw_info_.color_vel = ScalarFadeColor(velocity_color_, outer_iter->front().time_stamp, sys->current_time_, sys->params_.meas_time_window_);
        DrawTrack(img_, *model_iterator,sys, draw_info_);
    }

    reset_img_ = true;

}

//----------------------------------------------------------------------------------------------------------------------

template<typename tModel, template<typename > typename tDrawMeasurementPolicy, template<typename> typename tDrawTrackPolicy>
void VisualizationHost<tModel,tDrawMeasurementPolicy,tDrawTrackPolicy>::DrawSystemAndTrueTracks(const System<tModel>* sys, const std::vector<tModel>& tracks) {

    reset_img_ = false;
    img_ = cv::Mat( img_dimensions_(0), img_dimensions_(1), CV_8UC3 , cv::Scalar(255,255,255));


    // Draw true tracks

    for (auto model_iterator = tracks.begin(); model_iterator!= tracks.end(); ++ model_iterator) {


        draw_info_.color_pos = true_track_color_;
        draw_info_.color_vel = velocity_color_;
        DrawTrack(img_, *model_iterator,sys, draw_info_);
    }


    DrawSystem(sys);

}

//----------------------------------------------------------------------------------------------------------------------

template<typename tModel, template<typename > typename tDrawMeasurementPolicy, template<typename> typename tDrawTrackPolicy>

cv::Scalar VisualizationHost<tModel,tDrawMeasurementPolicy,tDrawTrackPolicy>::GetRandomColor() { 
    cv::Scalar new_color(dist_(rand_)*360),255,255);
    
    return ScalarHSV2BGR(new_color);
}



} // namespace rransac

#endif // RRANSAC_VISUALIZATION_VISUALIZATION_HOST_INFO_H_