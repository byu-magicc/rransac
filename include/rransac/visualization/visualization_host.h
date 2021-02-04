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
#include <string>

namespace rransac
{
    
template<typename tModel, template<typename > typename tDrawMeasurementPolicy, template<typename> typename tDrawTrackPolicy>
class VisualizationHost : tDrawMeasurementPolicy<tModel>, tDrawTrackPolicy<tModel> {

public:

VisualizationHost(const std::vector<int>& img_dimensions, const double drawing_scale);

VisualizationHost(const std::vector<int>& img_dimensions, const double drawing_scale, const std::string file_path);

VisualizationHost(const std::vector<int>& img_dimensions, const double drawing_scale, const std::string file_name, const double fps);

~VisualizationHost() {

    if (record_video_)
        video_writer_.release();        

}

void DrawSystem(const System<tModel>* sys);

void DrawSystemAndTrueTracks(const System<tModel>* sys, const std::vector<tModel>& tracks);

const cv::Mat* GetImg();

private:

std::vector<int> img_dimensions_;
double drawing_scale_;
cv::Mat img_;
cv::Scalar true_track_color_;
cv::Scalar estimated_poor_track_color_;
cv::Scalar estimated_good_track_color_;
cv::Scalar velocity_color_;
cv::Scalar measurement_color_;
std::vector<cv::Scalar> cluster_colors_;
std::default_random_engine gen_;
std::uniform_real_distribution<double> dist_;
DrawInfo draw_info_;
bool reset_img_;
unsigned int img_num_;
std::string window_name_;
std::string file_path_;
cv::VideoWriter video_writer_;
bool record_video_ = false;
double fps_;


static void DrawMeas(cv::Mat& img, const Meas<double>& meas, const System<tModel>* sys,  const DrawInfo& draw_info) {
    VisualizationHost::DrawMeasPolicy(img,meas,sys,draw_info);
}

static void DrawTrack(cv::Mat& img, const tModel& model, const System<tModel>* sys, const DrawInfo& draw_info) {
    VisualizationHost::DrawTrackPolicy(img,model,sys,draw_info);
} 

cv::Scalar GetRandomColor();

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
    double scale = 0.9/time_window*(time_stamp-curr_time+time_window) + 0.1;
    faded_color[1]*=scale;
    return ScalarHSV2BGR(faded_color);


}

void SaveImgCallbackFunction(int event, int x, int y, int flags, void* userdata ) {
    cv::imwrite(file_path_+ window_name_ + std::to_string(img_num_) + ".jpg", img_);
}

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                        Definitions
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename tModel, template<typename > typename tDrawMeasurementPolicy, template<typename> typename tDrawTrackPolicy>
VisualizationHost<tModel,tDrawMeasurementPolicy,tDrawTrackPolicy>::VisualizationHost(const std::vector<int>& img_dimensions, const double drawing_scale) :img_dimensions_(img_dimensions), drawing_scale_(drawing_scale) { 
    srand(time(NULL));
    dist_ = std::uniform_real_distribution<double>(0.0,1.0);
    true_track_color_ = cv::Scalar(0,255,0);
    estimated_poor_track_color_ = cv::Scalar(204,102,0);
    estimated_good_track_color_ = cv::Scalar(0,153,153);
    velocity_color_ = cv::Scalar(0,0,255);
    measurement_color_ = cv::Scalar(255,0,127);
    reset_img_ = true;

    draw_info_.img_center = cv::Point2d(img_dimensions_[1]/2,img_dimensions_[0]/2);
    draw_info_.scale_drawing = drawing_scale;
    draw_info_.draw_validation_region = true;
    draw_info_.scale_draw_pos = 2;
    draw_info_.scale_draw_vel = 2;
    img_num_ = 0;
    window_name_ = "RRANSAC Visualization";
    cv::namedWindow(window_name_, cv::WINDOW_AUTOSIZE);
    

    
}

//----------------------------------------------------------------------------------------------------------------------

template<typename tModel, template<typename > typename tDrawMeasurementPolicy, template<typename> typename tDrawTrackPolicy>
VisualizationHost<tModel,tDrawMeasurementPolicy,tDrawTrackPolicy>::VisualizationHost(const std::vector<int>& img_dimensions, const double drawing_scale, const std::string file_path) : VisualizationHost(img_dimensions, drawing_scale) {

    file_path_ = file_path;
    cv::setMouseCallback(window_name_,SaveImgCallbackFunction,NULL);

}

//----------------------------------------------------------------------------------------------------------------------

template<typename tModel, template<typename > typename tDrawMeasurementPolicy, template<typename> typename tDrawTrackPolicy>
VisualizationHost<tModel,tDrawMeasurementPolicy,tDrawTrackPolicy>::VisualizationHost(const std::vector<int>& img_dimensions, const double drawing_scale, const std::string file_name, const double fps) : VisualizationHost(img_dimensions, drawing_scale) {

    video_writer_.open(file_name,cv::VideoWriter::fourcc('m','p','4','v'),fps,cv::Size(img_dimensions[1],img_dimensions[0]));

    fps_ = fps;
    file_path_ = file_name;

    if (!video_writer_.isOpened())
        throw std::runtime_error("Could not open the output video for write. Filename:  " + file_name);

    record_video_ = true;

} 

//----------------------------------------------------------------------------------------------------------------------

template<typename tModel, template<typename > typename tDrawMeasurementPolicy, template<typename> typename tDrawTrackPolicy>
void VisualizationHost<tModel,tDrawMeasurementPolicy,tDrawTrackPolicy>::DrawSystem(const System<tModel>* sys) {

    draw_info_.draw_validation_region = true;
    if (reset_img_)
        img_ = cv::Mat( img_dimensions_[0], img_dimensions_[1], CV_8UC3 , cv::Scalar(255,255,255));


    // Draw measurements
    cv::Scalar color(0,0,0);
    cv::Scalar faded_color(0,0,0);
    
    for (auto cluster_iterator = sys->data_tree_.data_.begin(); cluster_iterator != sys->data_tree_.data_.end(); ++cluster_iterator) {

        if (cluster_iterator->cluster_label_ >= 0) {
            while (cluster_iterator->cluster_label_  >= cluster_colors_.size())
                cluster_colors_.push_back(GetRandomColor());
            color = cluster_colors_[cluster_iterator->cluster_label_ ];
        } else {
            color = measurement_color_;
        }

        for(auto outer_iter = cluster_iterator->data_.begin(); outer_iter != cluster_iterator->data_.end(); ++ outer_iter) {
            faded_color = ScalarFadeColor(color, outer_iter->front().time_stamp, sys->current_time_, sys->params_.meas_time_window_);
            // faded_color = color;
            draw_info_.color_pos = faded_color;
            draw_info_.color_vel = ScalarFadeColor(velocity_color_, outer_iter->front().time_stamp, sys->current_time_, sys->params_.meas_time_window_);


            for (auto inner_iter = outer_iter->begin(); inner_iter != outer_iter->end(); ++ inner_iter) {
                DrawMeas(img_, *inner_iter, sys, draw_info_);
            }
        }
    }

    // Draw tracks

    for (auto model_iterator = sys->models_.begin(); model_iterator!= sys->models_.end(); ++ model_iterator) {
        if (model_iterator->model_likelihood_ >=sys->params_.track_good_model_threshold_)
            color = estimated_good_track_color_;
        else 
            color = estimated_poor_track_color_;

        faded_color = ScalarFadeColor(color, sys->current_time_-model_iterator->missed_detection_time_, sys->current_time_, sys->params_.track_max_missed_detection_time_);
        draw_info_.color_pos = faded_color;
        draw_info_.color_vel = ScalarFadeColor(velocity_color_, sys->current_time_-model_iterator->missed_detection_time_, sys->current_time_, sys->params_.track_max_missed_detection_time_);
        DrawTrack(img_, *model_iterator,sys, draw_info_);
    }


    if (record_video_) {
        video_writer_.write(img_);
        cv::imshow(window_name_, img_);
        cv::waitKey(1);
    } else {
        cv::imshow(window_name_, img_);
        cv::waitKey(0);
    }

    
    
    reset_img_ = true;

}

//----------------------------------------------------------------------------------------------------------------------

template<typename tModel, template<typename > typename tDrawMeasurementPolicy, template<typename> typename tDrawTrackPolicy>
void VisualizationHost<tModel,tDrawMeasurementPolicy,tDrawTrackPolicy>::DrawSystemAndTrueTracks(const System<tModel>* sys, const std::vector<tModel>& tracks) {

    reset_img_ = false;
    img_ = cv::Mat( img_dimensions_[0], img_dimensions_[1], CV_8UC3 , cv::Scalar(255,255,255));


    // Draw true tracks
    draw_info_.draw_validation_region = false;
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
    cv::Scalar new_color(dist_(gen_)*360,255,255);
    
    return ScalarHSV2BGR(new_color);
}



} // namespace rransac

#endif // RRANSAC_VISUALIZATION_VISUALIZATION_HOST_INFO_H_