#include "common/sources/source_base.h"
#include <typeinfo>
#include "common/sources/source_base.h"

namespace rransac {

template<class S>
class SourceRN : public SourceBase<lie_groups::R2_r2> {

public:

/** Initializes the measurement source. This function must set the parameters.  */
// void Init(const SourceParameters& params) override;      

// /** Returns the jacobian of the observation function w.r.t. the states */
// Eigen::MatrixXd GetLinObsMatState(const S& state) override {return H_;}                              

// /** Returns the jacobian of the observation function w.r.t. the sensor noise */
// Eigen::MatrixXd GetLinObsMatSensorNoise(const S& state) override {return V_;}                         

// /** Computes the estimated measurement given a state */
// Eigen::MatrixXd GetEstMeas(const S& state) override {return state.g_.data_;} /** Returns an estimated measurement according to the state. */

void Test(){
    H_ = H_*2;
}

};


} // namesapce rransac