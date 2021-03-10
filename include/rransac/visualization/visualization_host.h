#ifndef RRANSAC_VISUALIZATION_VISUALIZATION_HOST_INFO_H_
#define RRANSAC_VISUALIZATION_VISUALIZATION_HOST_INFO_H_


#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <vector>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <random>
#include <string>

#include "rransac/system.h"
#include "rransac/visualization/draw_info.h"



namespace rransac
{
    
template<typename tModel, template<typename > typename tDrawMeasurementPolicy, template<typename> typename tDrawTrackPolicy>
class VisualizationHost : tDrawMeasurementPolicy<tModel>, tDrawTrackPolicy<tModel> {

public:

/**
 * This constructor requires the user to call one the setup function before drawing any measurements. Otherwise, bad things
 * can happen.
 */ 
VisualizationHost() = default;

/**
 * Calls the setup function Setup(const std::vector<int>& img_dimensions, const DrawInfo& draw_info)
 */ 
VisualizationHost(const std::vector<int>& img_dimensions, const DrawInfo& draw_info){Setup(img_dimensions,draw_info);}

/**
 * Calls the setup function Setup(const std::vector<int>& img_dimensions, const DrawInfo& draw_info, const std::string file_path)
 */ 
VisualizationHost(const std::vector<int>& img_dimensions, const DrawInfo& draw_info, const std::string file_path)
{Setup(img_dimensions,draw_info,file_path);}

/**
 * Calls the setup function Setup(const std::vector<int>& img_dimensions, const DrawInfo& draw_info, const std::string file_name, const double fps)
 */ 
VisualizationHost(const std::vector<int>& img_dimensions, const DrawInfo& draw_info, const std::string file_name, const double fps){
   Setup(img_dimensions,draw_info,file_name,fps);}


/**
 * Sets up the visualization. The visualization is shown on an image, and the thread spins until the user
 * has the image as the current window and presses a key.
 * @param draw_info The specifications for the visualizations. 
 * @param img_dimensions The dimensions of the image used in visualization. The first element is the width and the second is the height.
 */ 
void Setup(const std::vector<int>& img_dimensions, const DrawInfo& draw_info);

/**
 * Sets up the visualization. The visualization is shown on an image, and the thread spins until the user
 * has the image as the current window and presses a key. If the user clicks on the image, then the current image will be saved to
 * the provided file path. 
 * @param draw_info The specifications for the visualizations. 
 * @param img_dimensions The dimensions of the image used in visualization. The first element is the width and the second is the height.
 * @param file_path The path to the directory where the images should be saved. If the file path doesn't exist, an error will be thrown. 
 * The class will give the image a name with a numerical counter on it. 
 */ 
void Setup(const std::vector<int>& img_dimensions, const DrawInfo& draw_info, const std::string file_path);

/**
 * Sets up the visualization. The visualization is shown on an image, and the thread spins until the user
 * has the image as the current window and presses a key. The images will be saved as a video. An error will be thrown
 * if the file name is not valid (i.e. the directory protion of the file  name does not exist). 
 * @param draw_info The specifications for the visualizations. 
 * @param img_dimensions The dimensions of the image used in visualization. The first element is the width and the second is the height.
 * @param file_name The absolute file path and name of the video.
 * @param fps The frames per second for the video recording. 
 */ 
void Setup(const std::vector<int>& img_dimensions, const DrawInfo& draw_info, const std::string file_name, const double fps);



~VisualizationHost() {

    if (record_video_) {
        video_writer_.release(); 
        std::cerr << "video released " << std::endl;
    }
               

}

/**
 * Draws new measurements on the image using the color DrawInfo::color_new_meas_pos
 * @param new_measurements The new measurements to be drawn.
 * @param sys A pointer to the system information of R-RANSAC.
 * @param clear_img If set to true, the image will be cleared.
 */ 
void DrawNewMeasurements(const std::list<Meas<double>>& new_measurements, const System<tModel>* sys, bool clear_img);

/**
 * Draws new measurements on the image using the color DrawInfo::color_new_meas_pos
 * @param true_tracks The true tracks that R-RANSAC is tyring to track.
 * @param sys A pointer to the system information of R-RANSAC.
 * @param clear_img If set to true, the image will be cleared.
 */ 
void DrawTrueTracks(const std::vector<tModel>& true_tracks, const System<tModel>* sys , bool clear_img);

/**
 * Draws the measurements in the data tree that have not been associated with a track or cluster
 * @param sys A pointer to the system information of R-RANSAC that contains the measurements.
 * @param clear_img If set to true, the image will be cleared.
 */ 
void DrawUnAssociatedMeasurements(const System<tModel>* sys, bool clear_img);

/**
 * Draws the cluster, given each cluster a random color.
 * @param sys A pointer to the system information of R-RANSAC that contains the clusters.
 * @param clear_img If set to true, the image will be cleared.
 */ 
void DrawClusters(const System<tModel>* sys, bool clear_img);

/**
 * Draws the cluster, given each cluster a random color.
 * @param sys A pointer to the system information of R-RANSAC that contains the estimated tracks.
 * @param clear_img If set to true, the image will be cleared.
 */ 
void DrawEstimatedTracks(const System<tModel>* sys, bool clear_img);


/**
 * If the object is set up to record a video, calling this function will add the current image
 * to the video and show the current image.
 */ 
void RecordImage();

/**
 * Shows the current image and spins the thread until a key is pressed while the image viewier
 * is the current window. 
 */ 
void ShowImage(); 



const cv::Mat* GetImg();

private:

std::vector<int> img_dimensions_;
cv::Mat img_;
std::vector<cv::Scalar> cluster_colors_;
std::default_random_engine gen_;
std::uniform_real_distribution<double> dist_;
DrawInfo draw_info_;
DrawInfo draw_info_original_;
unsigned int img_num_;
std::string window_name_;
std::string file_path_;
cv::VideoWriter video_writer_;
bool record_video_ = false;
double fps_;


void ClearImage();

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
void VisualizationHost<tModel,tDrawMeasurementPolicy,tDrawTrackPolicy>::Setup(const std::vector<int>& img_dimensions, const DrawInfo& draw_info) { 
    srand(time(NULL));
    dist_ = std::uniform_real_distribution<double>(0.0,1.0);
    img_dimensions_ = img_dimensions;
    draw_info_ = draw_info;

    if(draw_info_.flip_image_x_axis) {
        draw_info_.flip_y = -1;
    } else {
        draw_info_.flip_y = 1;
    }

    draw_info_.img_center = cv::Point2d(img_dimensions_[0]/2,img_dimensions_[1]/2);
    draw_info_original_ = draw_info_;
    img_num_ = 0;
    window_name_ = "RRANSAC Visualization";
    cv::namedWindow(window_name_, cv::WINDOW_AUTOSIZE);

    ClearImage();
    

    
}

//----------------------------------------------------------------------------------------------------------------------

template<typename tModel, template<typename > typename tDrawMeasurementPolicy, template<typename> typename tDrawTrackPolicy>
void VisualizationHost<tModel,tDrawMeasurementPolicy,tDrawTrackPolicy>::Setup(const std::vector<int>& img_dimensions, const DrawInfo& draw_info, const std::string file_path) {


    Setup(img_dimensions,draw_info);
    file_path_ = file_path;
    cv::setMouseCallback(window_name_,SaveImgCallbackFunction,NULL);

}

//----------------------------------------------------------------------------------------------------------------------

template<typename tModel, template<typename > typename tDrawMeasurementPolicy, template<typename> typename tDrawTrackPolicy>
void VisualizationHost<tModel,tDrawMeasurementPolicy,tDrawTrackPolicy>::Setup(const std::vector<int>& img_dimensions, const DrawInfo& draw_info, const std::string file_name, const double fps){

    Setup(img_dimensions,draw_info);
    video_writer_.open(file_name,cv::VideoWriter::fourcc('m','p','4','v'),fps,cv::Size(img_dimensions[0],img_dimensions[1]));

    fps_ = fps;
    file_path_ = file_name;

    if (!video_writer_.isOpened())
        throw std::runtime_error("Could not open the output video for write. Filename:  " + file_name);

    record_video_ = true;

} 

//----------------------------------------------------------------------------------------------------------------------

template<typename tModel, template<typename > typename tDrawMeasurementPolicy, template<typename> typename tDrawTrackPolicy>
void VisualizationHost<tModel,tDrawMeasurementPolicy,tDrawTrackPolicy>::DrawNewMeasurements(const std::list<Meas<double>>& new_measurements, const System<tModel>* sys, bool clear_img) {

    if (clear_img)
        ClearImage();

    draw_info_.color_pos = draw_info_original_.color_new_meas_pos;
    draw_info_.draw_threshold = false;
    for (auto meas_iter = new_measurements.begin(); meas_iter != new_measurements.end(); ++meas_iter) {
        DrawMeas(img_, *meas_iter, sys, draw_info_);
    }

}

//----------------------------------------------------------------------------------------------------------------------

template<typename tModel, template<typename > typename tDrawMeasurementPolicy, template<typename> typename tDrawTrackPolicy>
void VisualizationHost<tModel,tDrawMeasurementPolicy,tDrawTrackPolicy>::DrawTrueTracks(const std::vector<tModel>& true_tracks, const System<tModel>* sys , bool clear_img) {

    if (clear_img)
        ClearImage();

    draw_info_.color_pos = draw_info_.color_true_track_pos;
    bool draw_validation_region = draw_info_.draw_validation_region;
    draw_info_.draw_validation_region = false;

    for (auto model_iterator = true_tracks.begin(); model_iterator!= true_tracks.end(); ++ model_iterator) {

        DrawTrack(img_, *model_iterator,sys, draw_info_);
    }

    draw_info_.draw_validation_region = draw_validation_region;
    

}

//----------------------------------------------------------------------------------------------------------------------

template<typename tModel, template<typename > typename tDrawMeasurementPolicy, template<typename> typename tDrawTrackPolicy>
void VisualizationHost<tModel,tDrawMeasurementPolicy,tDrawTrackPolicy>::DrawUnAssociatedMeasurements(const System<tModel>* sys, bool clear_img) {

    if (clear_img)
        ClearImage();

    
    draw_info_.draw_threshold = draw_info_.draw_measurment_velocity_position_threshold;
    
    for (auto cluster_iterator = sys->data_tree_.data_.begin(); cluster_iterator != sys->data_tree_.data_.end(); ++cluster_iterator) {

        if(cluster_iterator->cluster_label_ < 0) {
            for(auto outer_iter = cluster_iterator->data_.begin(); outer_iter != cluster_iterator->data_.end(); ++ outer_iter) {
                draw_info_.color_pos = ScalarFadeColor(draw_info_original_.color_meas_pos, outer_iter->front().time_stamp, sys->current_time_, sys->params_.meas_time_window_);
                draw_info_.color_vel = ScalarFadeColor(draw_info_original_.color_vel, outer_iter->front().time_stamp, sys->current_time_, sys->params_.meas_time_window_);

                for (auto inner_iter = outer_iter->begin(); inner_iter != outer_iter->end(); ++ inner_iter) {
                    DrawMeas(img_, *inner_iter, sys, draw_info_);
                }
            }
        }


    }

    draw_info_ = draw_info_original_;


}

//----------------------------------------------------------------------------------------------------------------------

template<typename tModel, template<typename > typename tDrawMeasurementPolicy, template<typename> typename tDrawTrackPolicy>
void VisualizationHost<tModel,tDrawMeasurementPolicy,tDrawTrackPolicy>::DrawClusters(const System<tModel>* sys, bool clear_img) {

    if (clear_img)
        ClearImage();

    // Draw measurements
    cv::Scalar color(0,0,0);
    draw_info_.draw_threshold = draw_info_.draw_cluster_velocity_position_threshold;

    for (auto cluster_iterator = sys->data_tree_.data_.begin(); cluster_iterator != sys->data_tree_.data_.end(); ++cluster_iterator) {
        if(cluster_iterator->cluster_label_  >= 0) {
            while (cluster_iterator->cluster_label_  >= cluster_colors_.size()) {
                cluster_colors_.push_back(GetRandomColor());
                if (cluster_iterator->cluster_label_ > 1e3)
                    break;
            }
                
            color = cluster_colors_[cluster_iterator->cluster_label_ ];

            for(auto outer_iter = cluster_iterator->data_.begin(); outer_iter != cluster_iterator->data_.end(); ++ outer_iter) {
                draw_info_.color_pos = ScalarFadeColor(color, outer_iter->front().time_stamp, sys->current_time_, sys->params_.meas_time_window_);
                draw_info_.color_vel = ScalarFadeColor(draw_info_original_.color_vel, outer_iter->front().time_stamp, sys->current_time_, sys->params_.meas_time_window_);

                for (auto inner_iter = outer_iter->begin(); inner_iter != outer_iter->end(); ++ inner_iter) {
                    DrawMeas(img_, *inner_iter, sys, draw_info_);
                }
            }
        }

    }

    draw_info_ = draw_info_original_;

}

//----------------------------------------------------------------------------------------------------------------------

template<typename tModel, template<typename > typename tDrawMeasurementPolicy, template<typename> typename tDrawTrackPolicy>
void VisualizationHost<tModel,tDrawMeasurementPolicy,tDrawTrackPolicy>::DrawEstimatedTracks(const System<tModel>* sys, bool clear_img) {


    if (draw_info_.draw_poor_tracks) {

        for (auto& track : sys->models_) {

            if(track.model_likelihood_ >= sys->params_.track_good_model_threshold_){
                draw_info_.color_pos = ScalarFadeColor(draw_info_original_.color_estimated_good_track_pos, sys->current_time_-track.missed_detection_time_, sys->current_time_, sys->params_.meas_time_window_);
            } else {
                draw_info_.color_pos = ScalarFadeColor(draw_info_original_.color_estimated_poor_track_pos, sys->current_time_-track.missed_detection_time_, sys->current_time_, sys->params_.meas_time_window_);
            }
            

            DrawTrack(img_, track,sys, draw_info_);
        }

    } else {
        for (auto&& track : sys->good_models_) {
            draw_info_.color_pos = ScalarFadeColor(draw_info_original_.color_estimated_good_track_pos, sys->current_time_-track->missed_detection_time_, sys->current_time_, sys->params_.meas_time_window_);
            draw_info_.color_vel = ScalarFadeColor(draw_info_original_.color_vel, sys->current_time_-track->missed_detection_time_, sys->current_time_, sys->params_.meas_time_window_);
            DrawTrack(img_, *track,sys, draw_info_);
        }
    }

    draw_info_ = draw_info_original_;

}


//----------------------------------------------------------------------------------------------------------------------

template<typename tModel, template<typename > typename tDrawMeasurementPolicy, template<typename> typename tDrawTrackPolicy>
cv::Scalar VisualizationHost<tModel,tDrawMeasurementPolicy,tDrawTrackPolicy>::GetRandomColor() { 
    cv::Scalar new_color(dist_(gen_)*360,255,255);
    
    return ScalarHSV2BGR(new_color);
}

//----------------------------------------------------------------------------------------------------------------------

template<typename tModel, template<typename > typename tDrawMeasurementPolicy, template<typename> typename tDrawTrackPolicy>
void VisualizationHost<tModel,tDrawMeasurementPolicy,tDrawTrackPolicy>::ClearImage() {
    img_ = cv::Mat( img_dimensions_[1], img_dimensions_[0], CV_8UC3 , cv::Scalar(255,255,255));
}

//----------------------------------------------------------------------------------------------------------------------

template<typename tModel, template<typename > typename tDrawMeasurementPolicy, template<typename> typename tDrawTrackPolicy>
void VisualizationHost<tModel,tDrawMeasurementPolicy,tDrawTrackPolicy>::RecordImage() {
    if (record_video_) {
        video_writer_.write(img_);
        cv::imshow(window_name_, img_);
        cv::waitKey(1);
    } 
}

//----------------------------------------------------------------------------------------------------------------------

template<typename tModel, template<typename > typename tDrawMeasurementPolicy, template<typename> typename tDrawTrackPolicy>
void VisualizationHost<tModel,tDrawMeasurementPolicy,tDrawTrackPolicy>::ShowImage() {
    cv::imshow(window_name_, img_);
    cv::waitKey(0);
}

} // namespace rransac

#endif // RRANSAC_VISUALIZATION_VISUALIZATION_HOST_INFO_H_