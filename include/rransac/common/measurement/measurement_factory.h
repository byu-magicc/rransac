#ifindef RRANSAC_COMMON_MEASUREMENT_MEASUREMENT_FACTORY_H_
#define RRANSAC_COMMON_MEASUREMENT_MEASUREMENT_FACTORY_H_

#include "rransac/common/measurement/measurement.h

namespace rransac
{
  /** \class MeasurementFactory
   * Creates a new instant of a derived measurement class.
  */
  class MeasurementFactory
  {
  public:

      /**
       * Returns a unique pointer to a derived measurement class.
       * @param[in] type The type of derived measurement class to use.
       * @see DerivedMeasurement
       */
      std::unique_ptr<MeasurementBase> CreateMeasurementClass(const DerivedMeasurement type);
  };
} // namespace rransac

#endif // RRANSAC_COMMON_MEASUREMENT_MEASUREMENT_FACTORY_H_
