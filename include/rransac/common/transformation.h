#ifndef RRANSAC_COMMON_TRANSFORMATION_H_
#define RRANSAC_COMMON_TRANSFORMATION_H_

#include <Eigen/Core>

namespace rransac
{
/** \struct Transformation
 * When the global frame changes, the measurements and models need to be transformed into the current global frame.
 * This struct provides the necessary data in order to transform the measurements and models. It is provided by the
 * user and used by R-RANSAC.
*/

struct Transformation
{
    Eigen::MatrixXd T;
};
} // namespace rransac

#endif // RRANSAC_COMMON_TRANSFORMATION_H_
