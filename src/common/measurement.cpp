#include "rransac/common/measurement/measurement.h"

namespace rransac
{

//-----------------------------------------------------------------------------
Measurement::Measurement(){};

//-----------------------------------------------------------------------------
Measurement::~Measurement(){};


//-----------------------------------------------------------------------------
static void Measurement::TransformMeasurement(Meas& meas, const Transformation& T, const double dt)
{
  //TODO:: Research transforming a measurement given a transformation matrix
}


//-----------------------------------------------------------------------------
static float Measurement::GetDistance(const Meas& meas1, const Meas& meas2, const DistanceType& type, const Parameters& params)
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
static float Measurement::GetSpatialDistance(Meas& meas1, Meas& meas2, const Parameters& params);
{

}

//-----------------------------------------------------------------------------
static float Measurement::GetTemporalDistance(Meas& meas1, Meas& meas2, const Parameters& params);
{

}

//-----------------------------------------------------------------------------
static float Measurement::GetTotalDistance(Meas& meas1, Meas& meas2, const Parameters& params);
{
  
}

} // namespace rransac
