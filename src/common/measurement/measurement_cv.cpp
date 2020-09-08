#include <math.h>

#include "rransac/common/measurement/measurement_cv.h"

namespace rransac
{

//-----------------------------------------------------------------------------

Measurement_CV::Measurement_CV() = default;

//-----------------------------------------------------------------------------

Measurement_CV::~Measurement_CV() = default;

//-----------------------------------------------------------------------------

void Measurement_CV::TransformMeasurement(Meas& meas, const Transformation& T, const double dt)
{
  meas.data = T.T * meas.data;
}

//-----------------------------------------------------------------------------

float Measurement_CV::GetDistance(const Meas& meas1, const Meas& meas2, const DistanceType& type, const Parameters& params)
{
  switch (type) {
    case kSpatial:
      return GetSpatialDistance(meas1, meas2, params);
    case kTemporal:
      return GetTemporalDistance(meas1, meas2, params);
    case kTotal:
      return GetTotalDistance(meas1, meas2, params);
  }
}

//-----------------------------------------------------------------------------

float Measurement_CV::GetSpatialDistance(const Meas& meas1, const Meas& meas2, const Parameters& params)
{
  // Euclidean/Frobenius norm distance between two points
  return (meas1.data - meas2.data).norm();
}

//-----------------------------------------------------------------------------

float Measurement_CV::GetTemporalDistance(const Meas& meas1, const Meas& meas2, const Parameters& params)
{
  return abs(meas1.time_stamp - meas2.time_stamp);
}

//-----------------------------------------------------------------------------

float Measurement_CV::GetTotalDistance(const Meas& meas1, const Meas& meas2, const Parameters& params)
{
  return sqrt(pow(GetSpatialDistance(meas1,meas2,params),2) + pow(GetTemporalDistance(meas1,meas2,params),2));
}

} // namespace rransac
