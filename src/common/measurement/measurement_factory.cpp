#include "rransac/common/measurement/measurement_factory.h"
#include "rransac/common/measurement/measurement_cv.h"

namespace rransac
{

//-----------------------------------------------------------------------------

MeasurementFactory::MeasurementFactory() = default;

//-----------------------------------------------------------------------------

MeasurementFactory::~MeasurementFactory() = default;

//-----------------------------------------------------------------------------

std::unique_ptr<MeasurementBase> MeasurementFactory::CreateMeasurementClass(const DerivedMeasurement type)
{
    switch (type) {
        case kMeasurementCV:
            return std::unique_ptr<MeasurementBase>{new Measurement_CV()};
    }
}


} // namespace rransac