#ifndef RRANSAC_VISUALIZATION_VISUALIZATION_HOST_INFO_H_
#define RRANSAC_VISUALIZATION_VISUALIZATION_HOST_INFO_H_

#include "system.h"
#include "draw_info.h"

namespace rransac
{
    
template<typename tModel, template<typename > typename tDrawMeasurementPolicy, template<typename> typename tDrawTrackPolicy>
class VisualizationHost : tDrawMeasurementPolicy<tModel>, tDrawTrackPolicy<tModel> {

public:



private:

static void DrawMeas(cv::Mat& img, const Meas<double>& meas, const System<tModel>& sys,  const DrawInfo& draw_info) {
    VisualizationHost::DrawMeasPolicy(img,meas,sys,draw_info);
}

static void DrawTrack(cv::Mat& img, const tModel& model, const System<tModel>& sys, const DrawInfo& draw_info) {
    VisualizationHost::DrawTrackPolicy(img,model,sys,draw_info);
} 


};


} // namespace rransac

#endif // RRANSAC_VISUALIZATION_VISUALIZATION_HOST_INFO_H_