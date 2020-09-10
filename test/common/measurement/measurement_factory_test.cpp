#include <gtest/gtest.h>

#include "common/measurement/measurement_factory.h"
#include "common/measurement/measurement_cv.h"

namespace rransac 
{

//---------------------------------------------------------------------------------------

TEST(Measurement_Factory_Test, uniq_ptr_test) 
{
    DerivedMeasurement meas_type = kMeasurementCV;
    MeasurementFactory factory_obj;

    std::unique_ptr<MeasurementBase> meas_obj = factory_obj.CreateMeasurementClass(meas_type);

    EXPECT_EQ(meas_obj->measurement_type_,"measurement_cv") << "Should have created ptr to measurement_cv class";
}


} // namespace rransac