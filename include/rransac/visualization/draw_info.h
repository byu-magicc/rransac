#ifndef RRANSAC_VISUALIZATION_DRAW_INFO_H_
#define RRANSAC_VISUALIZATION_DRAW_INFO_H_


#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>

namespace rransac {

struct DrawInfo {

cv::Scalar color_pos;
cv::Scalar color_vel;
double scale_draw_pos=1;   /** < When the object being drawn is a filled object, this value will scale the volume of the drawn object */
double scale_draw_vel=1;   /** < This scales the size of the velocity vectors that are being drawn */
double scale_drawing=1;    /** < Scales the position of objects being drawn by this value. It causes the image to expand or contract */
cv::Point img_center;
int line_thickness = 1;
bool draw_validation_region = false;

};

} // namespace rransac

#endif //RRANSAC_VISUALIZATION_DRAW_INFO_H_