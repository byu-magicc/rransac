#ifndef RRANSAC_VISUALIZATION_DRAW_MEAS_POLICIES_MEAS_R2_SE2_POS_POLICY_H_
#define RRANSAC_VISUALIZATION_DRAW_MEAS_POLICIES_MEAS_R2_SE2_POS_POLICY_H_


#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <vector>
#include <Eigen/Core>
#include "visualization/draw_info.h"
#include <cmath>
#include "system.h"

#include "common/measurement/measurement_base.h"

namespace rransac
{

template <typename tModel>    
class DrawMeasR2SE2PosPolicy {

public: 

static void DrawMeasPolicy(cv::Mat& img, const Meas<double>& meas, const System<tModel>& sys,  const DrawInfo& draw_info);

};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                        Definitions
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename tModel>    
void DrawMeasR2SE2PosPolicy<tModel>::DrawMeasPolicy(cv::Mat& img, const Meas<double>& meas, const System<tModel>& sys,  const DrawInfo& draw_info)
{

// Draw the position and velocity
double radius = 1;
cv::Point vel;
cv::Point pos(meas.pose(0),-meas.pose(1));
pos*=draw_info.scale_drawing;
pos += draw_info.img_center;

if (meas.twist.rows() >0) {
    vel = cv::Point(meas.twist(0),-meas.twist(1))*draw_info.scale_draw_vel;
}

cv::circle(img,pos,radius*draw_info.scale_draw_pos,draw_info.color_pos,-1, cv::LINE_AA);
cv::arrowedLine(img, pos,vel+pos,draw_info.color_vel,draw_info.line_thickness, cv::LINE_AA);

// Draw the position neighborhood

cv::circle(img,pos,sys.params_.cluster_position_threshold_*draw_info.scale_drawing,draw_info.color_pos,draw_info.line_thickness, cv::LINE_AA);

// Draw the velocity neighborhood
cv::circle(img,pos,sys.params_.cluster_velocity_threshold_*sys.params_.cluster_time_threshold_*draw_info.scale_drawing,draw_info.color_vel,draw_info.line_thickness, cv::LINE_AA);


}


} // namespace rransac


#endif // RRANSAC_VISUALIZATION_MEAS_R2_SE2_POS_POLICY_H_