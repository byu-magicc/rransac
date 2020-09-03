#include "rransac/common/measurement/measurement_cv.h"

namespace rransac
{

//-----------------------------------------------------------------------------

Measurement_CV::Measurement() = default;

//-----------------------------------------------------------------------------

Measurement_CV::~Measurement() = default;

//-----------------------------------------------------------------------------

static void Measurement_CV::TransformMeasurement(Meas& meas, const Transformation& T, const double dt)
{

}

//-----------------------------------------------------------------------------

static float Measurement_CV::GetDistance(const Meas& meas1, const Meas& meas2, const DistanceType& type, const Parameters& params)
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

static float Measurement_CV::GetSpatialDistance(Meas& meas1, Meas& meas2, const Parameters& params);
{
  // Euclidean/Frobenius norm distance between two points
  return (meas1.data - meas2.data).norm();
}

//-----------------------------------------------------------------------------

static float Measurement_CV::GetTemporalDistance(Meas& meas1, Meas& meas2, const Parameters& params);
{
  return meas1.time_stamp - meas2.time_stamp;
}

//-----------------------------------------------------------------------------

static float Measurement_CV::GetTotalDistance(Meas& meas1, Meas& meas2, const Parameters& params);
{
  return sqrt(GetSpatialDistance(meas1,meas2,params) ^ 2 + GetTemporalDistance(meas1,meas2,params) ^ 2);
}

} // namespace rransac
