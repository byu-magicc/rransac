#ifndef RRANSAC_VISUALIZATION_DRAW_INFO_H_
#define RRANSAC_VISUALIZATION_DRAW_INFO_H_


#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>

namespace rransac {

struct DrawInfo {

cv::Scalar color_pos;
cv::Scalar color_vel;
double scale_pos=1;
double scale_vel=1;
cv::Point img_center;
int line_thickness = 1;

};

} // namespace rransac

#endif //RRANSAC_VISUALIZATION_DRAW_INFO_H_