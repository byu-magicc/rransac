#ifndef RRANSAC_SYSTEM_H_
#define RRANSAC_SYSTEM_H_

#include <vector>
#include <list>

#include "parameters.h"
#include "common/measurement/measurement_base.h"

// // Sources
// #include "common/sources/source_base.h"
// #include "common/sources/source_RN.h"
// #include "common/sources/source_SEN_pos_vel.h"
// #include "common/sources/source_SEN_pose_twist.h"

// // Models
// #include "common/models/model_base.h"
// #include "common/models/model_RN.h"
// #include "common/models/model_SEN_pos_vel.h"
// #include "common/models/model_SEN_pose_twist.h"

// // Transformations
// #include "common/transformations/transformation_base.h"
// #include "common/transformations/transformation_null.h"
// #include "common/transformations/trans_homography.h"

namespace rransac {


/**
 * \class System
 * This class holds all of the system information including parameters, sources, and models. 
 * 
 */ 
template <typename tModel> 
class System {

typedef tModel Model;

public:

    Parameters params_;                                        /** < System parameters */
    std::vector<typename tModel::Source> sources_;             /** < Contains all of the instantiated sources. */
    std::vector<std::vector<Meas>> new_meas_;                  /** < Contains all of the new measurements. */
    typename tModel::Transformation transformaion_;            /** < The transformation for the measurements and tracks */
    std::list<tModel> models_;                                 /** < The models created by rransac */
    std::list<std::unique_ptr<tModel>> good_models_;           /** < A list of pointers to the good models */

    // TODO add clusters and data tree

private:

};



} // rransac

#endif // RRANSAC_SYSTEM_H_
