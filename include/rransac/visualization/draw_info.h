#ifndef RRANSAC_VISUALIZATION_DRAW_INFO_H_
#define RRANSAC_VISUALIZATION_DRAW_INFO_H_


#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>

namespace rransac {

struct DrawInfo {

cv::Scalar color_meas_pos = cv::Scalar(255,0,127);                /**< The color of measurements not assigned to a cluster. */
cv::Scalar color_new_meas_pos = cv::Scalar(0,0,0);                /**< The color of new measurements not associated to a track or cluster. */
cv::Scalar color_estimated_poor_track_pos=cv::Scalar(204,102,0);  /**< The color of the objects representing the poor track's pose. */
cv::Scalar color_estimated_good_track_pos=cv::Scalar(0,153,153);  /**< The color of the objects representing the good track's pose. */
cv::Scalar color_true_track_pos=cv::Scalar(0,255,0);              /**< The color of the objects representing the true track's pose. */
cv::Scalar color_vel = cv::Scalar(0,0,255);                       /** < Every object drawn associated with velocity will be this color. */
double scale_draw_pos=2;                                          /** < When the object being drawn is a filled object, this value will scale the volume of the drawn object */
double scale_draw_vel=2;                                          /** < This scales the size of the velocity vectors that are being drawn */
double scale_drawing =1;                                          /** < Scales the position of objects being drawn by this value. It causes the image to expand or contract */
int line_thickness  =1;                                           /** < The thickness of the lines drawn for velocity, validation region, and thresholds. */
bool draw_validation_region=false;                                /** < Draws the validation region of the tracks. */
bool draw_cluster_velocity_position_threshold=false;              /** < Draw a circle around each measurement in a cluster to indicate the velocity and position threshold. The velocity threshold will be the color of color_vel 
                                                                        and the position threshold will be a randomly generated color. */
bool draw_measurment_velocity_position_threshold=false;           /** < Draw a circle around each measurement not in a cluster to indicate the velocity and position threshold. The velocity threshold will be the color of color_vel 
                                                                        and the position threshold will be the color of. */

bool draw_poor_tracks = false;                                    /**< Draw the poor tracks in addition to the good tracks. */

// Not user defined
cv::Point img_center;                                             /** < The center of the image. */
cv::Scalar color_pos;                                             /** < The color of the object representing pose that is about to be drawn. */
bool draw_threshold;

};

} // namespace rransac

#endif //RRANSAC_VISUALIZATION_DRAW_INFO_H_


